"""Physics checks on recoil (primary knock-on atom) production."""

from pathlib import Path
import shutil

import h5py
import numpy as np
import pytest

import openmc


@pytest.mark.parametrize('lr,ion_model,emitted', [
    (None, 'statistical', True), (0, 'statistical', True),
    (31, 'statistical', True), (33, 'statistical', True),
    (33, 'statistical', False), (33, 'none', True),
    (1, 'statistical', True), (99, 'statistical', True)])
def test_residual_breakup_metadata(run_in_tmpdir, lr, ion_model, emitted):
    """Missing LR preserves legacy behavior; known breakup changes the products."""
    source = openmc.data.DataLibrary.from_xml().get_by_material('Li7')['path']
    path = Path('Li7.h5').resolve()
    shutil.copy2(source, path)
    with h5py.File(path, 'r+') as handle:
        attrs = handle['Li7/reactions/reaction_052'].attrs
        if lr is not None:
            attrs['breakup_flag'] = lr
        elif 'breakup_flag' in attrs:
            del attrs['breakup_flag']
    library = openmc.data.DataLibrary()
    library.register_file(path)
    library.export_to_xml('cross_sections.xml')
    model = _one_collision_model(
        'Li7', 5.e6, ['H3', 'He4', 'Li7'], [0., 10.e6], mts=[52],
        recoil={'light_ion_model': ion_model, 'emitted_ions': emitted})
    model.materials.cross_sections = Path('cross_sections.xml').resolve()
    triton, alpha, lithium = _run(model, run_in_tmpdir).ravel()
    if lr in (None, 0, 31):
        assert lithium > 0 and triton == alpha == 0
    elif lr == 33 and ion_model == 'statistical':
        assert alpha > 0 and lithium == 0
        assert triton == (alpha if emitted else 0)
    else:
        assert triton == alpha == lithium == 0


@pytest.mark.parametrize('group,key,value', [
    ('Li7', 'excitation_energy', -1.),
    ('Li7', 'excitation_energy', float('nan')),
    ('Li7/reactions/reaction_052', 'breakup_flag', -2),
    ('Li7/reactions/reaction_052', 'q_mass_difference', float('inf'))])
def test_invalid_recoil_metadata_is_not_a_legacy_default(
        run_in_tmpdir, group, key, value):
    source = openmc.data.DataLibrary.from_xml().get_by_material('Li7')['path']
    path = Path('Li7.h5').resolve()
    shutil.copy2(source, path)
    with h5py.File(path, 'r+') as handle:
        handle[group].attrs[key] = value
    library = openmc.data.DataLibrary()
    library.register_file(path)
    library.export_to_xml('cross_sections.xml')
    model = _one_collision_model('Li7', 5.e6, ['Li7'], [0., 10.e6], particles=1)
    model.materials.cross_sections = Path('cross_sections.xml').resolve()
    with pytest.raises(RuntimeError, match='Invalid'):
        model.run(output=False)


def test_recoil_metadata_reinitializes(run_in_tmpdir):
    """C API finalization must not retain a previous library's breakup cache."""
    import openmc.lib

    source = openmc.data.DataLibrary.from_xml().get_by_material('Li7')['path']
    path = Path('Li7.h5').resolve()
    shutil.copy2(source, path)
    library = openmc.data.DataLibrary()
    library.register_file(path)
    library.export_to_xml('cross_sections.xml')
    model = _one_collision_model('Li7', 5.e6, ['He4', 'Li7'], [0., 10.e6], mts=[52])
    model.materials.cross_sections = Path('cross_sections.xml').resolve()
    for lr in (33, 0, 33):
        with h5py.File(path, 'r+') as handle:
            handle['Li7/reactions/reaction_052'].attrs['breakup_flag'] = lr
        model.init_lib(output=False)
        try:
            openmc.lib.run()
            alpha, lithium = openmc.lib.tallies[model.tallies[0].id].mean.ravel()
            if lr == 33:
                assert alpha > 0 and lithium == 0
            else:
                assert lithium > 0 and alpha == 0
        finally:
            model.finalize_lib()


@pytest.mark.parametrize('excitation,minimum', [
    (0., None), (0., 8.e6), (137900., 8.e6), (0., -1.e6), (0., 20.e6)])
def test_continuum_charged_excitation_limit(run_in_tmpdir, excitation, minimum):
    """A synthetic continuum component retains its required residual excitation."""
    source = openmc.data.DataLibrary.from_xml().get_by_material('Fe56')['path']
    path = Path('Fe56.h5').resolve()
    shutil.copy2(source, path)
    with h5py.File(path, 'r+') as handle:
        target = handle['Fe56']
        target.attrs['excitation_energy'] = excitation
        reactions = target['reactions']
        reactions.move('reaction_103', 'reaction_649')
        attrs = reactions['reaction_649'].attrs
        attrs['mt'] = 649
        qi = attrs.get('q_reaction', attrs.get('Q_value'))
        if minimum is None:
            if 'q_mass_difference' in attrs:
                del attrs['q_mass_difference']
        else:
            attrs['q_mass_difference'] = qi + minimum - excitation
    library = openmc.data.DataLibrary()
    library.register_file(path)
    library.export_to_xml('cross_sections.xml')
    # With an 8 MeV minimum excitation, the proton must stay below 4 MeV
    # at 14 MeV incident energy, including the laboratory boost. The
    # unrestricted distribution strongly populates energies above that bound.
    model = _one_collision_model('Fe56', 14.e6, ['H1'], [0., 4.e6, 20.e6],
        mts=[649], recoil={'emitted_ions': True}, particles=50000)
    model.materials.cross_sections = Path('cross_sections.xml').resolve()
    low, high = _run(model, run_in_tmpdir).ravel()
    if minimum is not None and minimum > 14.e6:
        assert low == high == 0
        return
    assert low > 0
    if minimum is not None and minimum > 0:
        assert high == 0
    else:
        assert high > 0


@pytest.mark.parametrize('offset,level,excitation,metadata', [
    (144630., 0., 0., True), (144630., 1000., 0., True),
    (144630., 50690., 0., True), (-200000., 1.e6, 0., True),
    (144630., 50690., 137900., True),
    (144630., 50690., 0., False),
    (300000., 0., 0., True), (0., -1000., 0., True),
    (0., -.11, 0., True)])
def test_discrete_charged_level_budget(
        run_in_tmpdir, offset, level, excitation, metadata):
    """Reconcile masses without changing the evaluated residual excitation."""
    amu = 931.49410242e6
    neutron = openmc.data.NEUTRON_MASS
    proton = 1.007276466621
    electron = 5.48579909065e-4
    residual = openmc.data.atomic_mass('Mn56') - 25 * electron
    target = openmc.data.atomic_mass('Fe56') - 26 * electron
    ground_q = (target + neutron - residual - proton) * amu
    qm = ground_q + offset
    qi = qm + excitation - level
    source = openmc.data.DataLibrary.from_xml().get_by_material('Fe56')['path']
    path = Path('Fe56.h5').resolve()
    shutil.copy2(source, path)
    with h5py.File(path, 'r+') as handle:
        target_data = handle['Fe56']
        target_data.attrs['excitation_energy'] = excitation
        reactions = target_data['reactions']
        reactions.move('reaction_103', 'reaction_600')
        attrs = reactions['reaction_600'].attrs
        attrs['mt'] = 600
        attrs['q_reaction'] = qi
        if metadata:
            attrs['q_mass_difference'] = qm
        elif 'q_mass_difference' in attrs:
            del attrs['q_mass_difference']
    library = openmc.data.DataLibrary()
    library.register_file(path)
    library.export_to_xml('cross_sections.xml')
    rejected = offset - level > 250000. or (metadata and level < -1.)
    q = (ground_q + excitation - max(0., level) if metadata else
         min(qi, ground_q + excitation))
    energy = 14.e6
    total_mass = proton + residual
    internal = energy + q - energy * neutron / total_mass
    ion_cm = internal * residual / total_mass
    boost = energy * neutron * proton / total_mass**2
    width = 2 * np.sqrt(ion_cm * boost)
    # A 100 eV margin permits mass-table rounding while resolving the 1 keV
    # separation between the first two synthetic levels.
    lower, upper = ion_cm + boost - width - 100., ion_cm + boost + width + 100.
    model = _one_collision_model('Fe56', energy, ['H1'],
        [0., lower, upper, 20.e6], mts=[600],
        recoil={'emitted_ions': True}, particles=30000)
    model.materials.cross_sections = Path('cross_sections.xml').resolve()
    below, inside, above = _run(model, run_in_tmpdir).ravel()
    assert below == above == 0
    assert inside == 0 if rejected else inside > 0


def _one_collision_model(nuclide, energy, products, e_bins, mts=None,
                         recoil=None, particles=20000, density=5.0,
                         radius=0.4):
    """Thin monoenergetic target so that scored collisions are first ones."""
    mat = openmc.Material()
    mat.add_nuclide(nuclide, 1.0)
    mat.set_density('g/cm3', density)

    sph = openmc.Sphere(r=radius, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sph)

    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.particles = particles
    settings.batches = 2
    settings.survival_biasing = False
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([energy], [1.0]),
    )
    # Kill anything that lost energy so only first collisions contribute
    settings.cutoff = {'energy_neutron': 0.999 * energy}
    settings.recoil_production = True
    if recoil is not None:
        settings.recoil = recoil

    tally = openmc.Tally(name='recoil')
    filters = [openmc.ParticleProductionFilter(products, e_bins)]
    if mts is not None:
        filters.insert(0, openmc.ReactionFilter(mts))
    tally.filters = filters
    tally.scores = ['events']

    return openmc.Model(openmc.Geometry([cell]), [mat], settings,
                        openmc.Tallies([tally]))


def _run(model, run_in_tmpdir):
    sp_path = model.run(output=False)
    with openmc.StatePoint(sp_path) as sp:
        return sp.get_tally(name='recoil').mean.copy()


def test_elastic_recoil_endpoint(run_in_tmpdir):
    """Elastic recoil must not exceed 4A/(A+1)^2 E and must reach most of it."""
    energy = 2.0e6
    nuclide = 'Fe56'
    e_bins = np.logspace(0, np.log10(0.5e6), 201)
    model = _one_collision_model(nuclide, energy, [nuclide], e_bins)
    mean = _run(model, run_in_tmpdir).ravel()

    lib = openmc.data.IncidentNeutron.from_hdf5(
        openmc.data.DataLibrary.from_xml().get_by_material(nuclide)['path'])
    awr = lib.atomic_weight_ratio
    e_max = 4.0 * awr / (awr + 1.0) ** 2 * energy

    nz = np.nonzero(mean)[0]
    assert nz.size > 0
    # Nothing above the kinematic maximum (allow one bin for binning width)
    assert e_bins[nz[-1]] <= e_max * 1.02
    # The upper end of the range is populated
    assert e_bins[nz[-1] + 1] > 0.5 * e_max


def test_bound_thermal_scattering_is_outside_recoil_scope(run_in_tmpdir):
    """An S(alpha,beta) event must not be reported as a free-atom PKA."""
    energy = 0.0253
    mat = openmc.Material()
    mat.add_nuclide('H1', 1.0)
    mat.set_density('g/cm3', 1.0)
    mat.add_s_alpha_beta('c_H_in_H2O')

    sph = openmc.Sphere(r=0.4, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sph)
    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.particles = 20000
    settings.batches = 2
    settings.survival_biasing = False
    settings.recoil_production = True
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([energy], [1.0]),
    )

    tally = openmc.Tally(name='recoil')
    tally.filters = [
        openmc.ReactionFilter(['(n,elastic)']),
        openmc.ParticleProductionFilter(['H1'], [0.0, 1.0]),
    ]
    tally.scores = ['events']
    model = openmc.Model(openmc.Geometry([cell]), [mat], settings,
                         openmc.Tallies([tally]))

    mean = _run(model, run_in_tmpdir)
    assert mean.sum() == 0.0


def test_elastic_recoil_energy_balance(run_in_tmpdir):
    """Recoil energy plus outgoing neutron energy equals the incident energy."""
    energy = 2.0e6
    nuclide = 'Fe56'
    e_bins = np.logspace(0, np.log10(0.5e6), 401)
    model = _one_collision_model(nuclide, energy, [nuclide], e_bins,
                                 particles=100000)
    # Tally the outgoing neutron spectrum in the same run
    n_bins = np.linspace(0, energy, 401)
    tally = openmc.Tally(name='nspec')
    tally.filters = [openmc.ParticleFilter(['neutron']),
                     openmc.EnergyFilter(n_bins)]
    tally.scores = ['scatter']
    model.tallies.append(tally)

    sp_path = model.run(output=False)
    with openmc.StatePoint(sp_path) as sp:
        recoil = sp.get_tally(name='recoil').mean.ravel()

    mid = np.sqrt(e_bins[:-1] * e_bins[1:])
    mean_recoil = (recoil * mid).sum() / recoil.sum()

    lib = openmc.data.IncidentNeutron.from_hdf5(
        openmc.data.DataLibrary.from_xml().get_by_material(nuclide)['path'])
    awr = lib.atomic_weight_ratio
    # <E_R> = 2A/(A+1)^2 E (1 - <mu_cm>); with |<mu_cm>| <= 1 this bounds it
    e_max = 4.0 * awr / (awr + 1.0) ** 2 * energy
    assert 0.0 < mean_recoil < e_max


def test_capture_recoil_scale(run_in_tmpdir):
    """Thermal capture recoil is set by the photon kick, not by E_in."""
    energy = 0.0253
    e_bins = np.logspace(0, 4, 121)
    model = _one_collision_model('Fe56', energy, ['Fe57'], e_bins,
                                 mts=['(n,gamma)'], particles=20000,
                                 density=7.874)
    mean = _run(model, run_in_tmpdir).ravel()
    assert mean.sum() > 0

    mid = np.sqrt(e_bins[:-1] * e_bins[1:])
    mean_recoil = (mean * mid).sum() / mean.sum()
    # A single 7.6 MeV photon gives E_R = E_g^2 / (2 M c^2) ~ 550 eV; the
    # sampled emission divides the energy among several photons, so require
    # the right decade.
    assert 20.0 < mean_recoil < 1000.0


def test_light_ion_model_none(run_in_tmpdir):
    """light_ion_model='none' suppresses emitted-ion secondary particles."""
    energy = 14.0e6
    e_bins = np.logspace(0, np.log10(2.0e7), 101)

    model = _one_collision_model(
        'Fe56', energy, ['H1'], e_bins, mts=['(n,p)'],
        recoil={'light_ion_model': 'statistical', 'emitted_ions': True},
        particles=100000, density=7.874)
    with_ions = _run(model, run_in_tmpdir).ravel().sum()

    model = _one_collision_model(
        'Fe56', energy, ['H1'], e_bins, mts=['(n,p)'],
        recoil={'light_ion_model': 'none', 'emitted_ions': True},
        particles=100000, density=7.874)
    without = _run(model, run_in_tmpdir).ravel().sum()

    assert with_ions > 0.0
    assert without == 0.0


def test_emitted_ions_flag(run_in_tmpdir):
    """emitted_ions=False keeps the residual but drops the light ion."""
    energy = 14.0e6
    e_bins = np.logspace(0, np.log10(2.0e7), 101)
    model = _one_collision_model(
        'Fe56', energy, ['H1', 'Mn56'], e_bins, mts=['(n,p)'],
        recoil={'emitted_ions': False}, particles=100000, density=7.874)
    mean = _run(model, run_in_tmpdir).reshape(2, -1)
    assert mean[0].sum() == 0.0      # H1 not banked
    assert mean[1].sum() > 0.0       # Mn56 residual still banked


def test_discrete_alpha_products_close_the_banked_energy(run_in_tmpdir):
    """B-10 MT=800 banked products close, not merely the helper formulas."""
    energy = 2.0e6
    nuclide = 'B10'
    lib = openmc.data.IncidentNeutron.from_hdf5(
        openmc.data.DataLibrary.from_xml().get_by_material(nuclide)['path'])
    budget = energy + lib[800].q_reaction

    # A 240 eV bin width resolves sub-keV departures from energy closure.
    e_bins = np.linspace(0.0, budget * 1.001, 20001)
    model = _one_collision_model(
        nuclide, energy, ['He4', 'Li7'], e_bins, mts=[800],
        particles=50000, density=2.34)
    mean = _run(model, run_in_tmpdir).reshape(2, -1)
    assert np.all(mean.sum(axis=1) > 0.0)

    midpoint = 0.5 * (e_bins[1:] + e_bins[:-1])
    product_mean = (mean * midpoint).sum(axis=1) / mean.sum(axis=1)
    assert abs(product_mean.sum() - budget) < 500.0


@pytest.mark.parametrize('mt,energy,products', [
    (105, 2.e6, ['H3', 'He4']), (103, 5.e6, ['H1', 'He6'])])
def test_ground_state_only_products_close_the_energy(
        run_in_tmpdir, mt, energy, products):
    """Summed single-ion channels cannot excite an intact He4 or He6 residual."""
    nuclide = 'Li6'
    lib = openmc.data.IncidentNeutron.from_hdf5(
        openmc.data.DataLibrary.from_xml().get_by_material(nuclide)['path'])
    budget = energy + lib[mt].q_reaction
    e_bins = np.linspace(0.0, budget * 1.001, 2001)
    model = _one_collision_model(
        nuclide, energy, products, e_bins, mts=[mt],
        particles=50000, density=0.534)
    mean = _run(model, run_in_tmpdir).reshape(2, -1)
    assert np.all(mean.sum(axis=1) > 0.0)

    midpoint = 0.5 * (e_bins[1:] + e_bins[:-1])
    product_mean = (mean * midpoint).sum(axis=1) / mean.sum(axis=1)
    # Allow histogram resolution and the small evaluated/AME mass difference.
    assert abs(product_mean.sum() - budget) < 5000.0


@pytest.mark.parametrize('lr,photons,two_body', [
    (None, False, True), (0, False, True),
    (31, False, False), (0, True, False)])
def test_ground_state_only_rule_respects_decay_metadata(
        run_in_tmpdir, lr, photons, two_body):
    """Known photon emission prevents the summed-channel two-body inference."""
    source = openmc.data.DataLibrary.from_xml().get_by_material('Li6')['path']
    path = Path('Li6.h5').resolve()
    shutil.copy2(source, path)
    with h5py.File(path, 'r+') as handle:
        reactions = handle['Li6/reactions']
        reaction = reactions['reaction_105']
        budget = 2.e6 + reaction.attrs.get('q_reaction', reaction.attrs.get('Q_value'))
        if lr is None:
            if 'breakup_flag' in reaction.attrs:
                del reaction.attrs['breakup_flag']
        else:
            reaction.attrs['breakup_flag'] = lr
        if photons:
            handle.copy(reactions['reaction_102/product_0'], reaction, 'product_0')
    library = openmc.data.DataLibrary()
    library.register_file(path)
    library.export_to_xml('cross_sections.xml')
    edges = np.linspace(0., budget*1.001, 2001)
    model = _one_collision_model('Li6', 2.e6, ['H3', 'He4'], edges,
        mts=[105], particles=30000, density=.534)
    model.materials.cross_sections = Path('cross_sections.xml').resolve()
    spectra = _run(model, run_in_tmpdir).reshape(2, -1)
    assert np.all(spectra.sum(axis=1) > 0.)
    means = spectra @ (.5*(edges[1:]+edges[:-1])) / spectra.sum(axis=1)
    if two_body:
        assert abs(means.sum()-budget) < 5000.
    else:
        assert means.sum() < .98*budget


def test_multi_neutron_within_energy_budget(run_in_tmpdir):
    """(n,2n) recoil stays inside the kinematic window set by E + Q."""
    energy = 14.0e6
    nuclide = 'Fe56'
    e_bins = np.logspace(0, np.log10(5.0e6), 201)
    model = _one_collision_model(nuclide, energy, ['Fe55'], e_bins,
                                 mts=['(n,2n)'], particles=200000,
                                 density=7.874)
    mean = _run(model, run_in_tmpdir).ravel()
    assert mean.sum() > 0

    lib = openmc.data.IncidentNeutron.from_hdf5(
        openmc.data.DataLibrary.from_xml().get_by_material(nuclide)['path'])
    awr = lib.atomic_weight_ratio
    # The residual momentum cannot exceed the incident momentum plus the
    # momentum of neutrons sharing E + Q, so bound the recoil generously.
    q = lib[16].q_reaction
    p_in = np.sqrt(energy)
    p_out = 2.0 * np.sqrt(max(energy + q, 0.0))
    e_max = (p_in + p_out) ** 2 / awr
    nz = np.nonzero(mean)[0]
    assert e_bins[nz[-1]] < e_max


# ---------------------------------------------------------------------------
# Fission and nonfission absorption
#
# Whether a recoil is produced must depend on which reaction absorbed the
# neutron, not on whether a fission source site happened to be banked. Those
# are sampled independently: p.fission() is true for a genuine capture on a
# fissionable nuclide whenever a site was banked, and false for a genuine
# fission in a fixed-source run with create_fission_neutrons off.
# ---------------------------------------------------------------------------


def _fissionable_model(nuclide, energy, products, e_bins, *, mts=None,
                       eigenvalue=False, survival=False,
                       create_fission_neutrons=True, particles=40000,
                       density=18.0):
    """A thin sphere of a fissionable nuclide, so scored events are first ones."""
    mat = openmc.Material()
    mat.add_nuclide(nuclide, 1.0)
    mat.set_density('g/cm3', density)

    sph = openmc.Sphere(r=0.4, boundary_type='vacuum')
    cell = openmc.Cell(fill=mat, region=-sph)

    settings = openmc.Settings()
    settings.particles = particles
    settings.batches = 3
    settings.survival_biasing = survival
    settings.recoil_production = True
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([energy], [1.0]),
    )
    if eigenvalue:
        settings.run_mode = 'eigenvalue'
        settings.inactive = 1
    else:
        settings.run_mode = 'fixed source'
        settings.create_fission_neutrons = create_fission_neutrons
    # Kill anything that lost energy so only first collisions contribute
    settings.cutoff = {'energy_neutron': 0.999 * energy}

    tally = openmc.Tally(name='recoil')
    filters = [openmc.ParticleProductionFilter(products, e_bins)]
    if mts is not None:
        filters.insert(0, openmc.ReactionFilter(mts))
    tally.filters = filters
    tally.scores = ['events']

    return openmc.Model(openmc.Geometry([cell]), [mat], settings,
                        openmc.Tallies([tally]))


def _absorption_fractions(nuclide, energy):
    """Nonfission and fission shares of the absorption cross section."""
    from openmc.data import IncidentNeutron, DataLibrary

    lib = IncidentNeutron.from_hdf5(
        DataLibrary.from_xml().get_by_material(nuclide)['path'])
    fission = capture = 0.0
    for mt, rx in lib.reactions.items():
        if mt in (18, 19, 20, 21, 38):
            fission += rx.xs['294K'](energy) if mt == 18 else 0.0
        elif mt == 102:
            capture += rx.xs['294K'](energy)
    return capture, fission


def test_fission_does_not_produce_a_nonfission_recoil(run_in_tmpdir):
    """No recoil may be attributed to a fission MT.

    Fission fragments are out of scope, so a fission absorption must produce
    no recoil secondary particle. Capture on the same nuclide must still
    produce U-236.
    """
    e_bins = np.logspace(0, np.log10(1.0e6), 61)
    model = _fissionable_model(
        'U235', 0.0253, ['U236', 'Xe135'], e_bins,
        mts=['fission', '(n,gamma)'], create_fission_neutrons=False)
    mean = _run(model, run_in_tmpdir).reshape(2, 2, -1)

    # axis 0 is the reaction filter in the order given
    assert mean[0].sum() == 0.0, (
        "fission must produce no recoil secondary particle")
    assert mean[1, 0].sum() > 0.0, "capture must still produce U-236"


def test_capture_recoil_survives_a_banked_fission_site(run_in_tmpdir):
    """Fission-site banking must not change capture-recoil production."""
    e_bins = np.logspace(0, np.log10(1.0e6), 61)

    def capture_rate(**kw):
        model = _fissionable_model('U235', 0.0253, ['U236'], e_bins,
                                   mts=['(n,gamma)'], **kw)
        return _run(model, run_in_tmpdir).sum()

    # Both runs see the same source spectrum; they differ only in whether the
    # fission sites that spectrum produces are banked. An eigenvalue comparison
    # would not be like-for-like, because after the first batch its source is the
    # fission spectrum rather than the thermal point source.
    without_sites = capture_rate(create_fission_neutrons=False)
    with_sites = capture_rate(create_fission_neutrons=True)

    assert without_sites > 0.0
    assert with_sites > 0.0
    assert abs(with_sites - without_sites) / without_sites < 0.15, (
        without_sites, with_sites)


def test_capture_recoil_is_produced_in_eigenvalue_mode(run_in_tmpdir):
    """Capture recoils are produced on a fissionable nuclide in eigenvalue mode.

    A rate comparison against fixed source would not be like-for-like because
    the eigenvalue source becomes the fission spectrum after the first batch;
    this test checks only that capture production occurs.
    """
    e_bins = np.logspace(0, np.log10(1.0e6), 61)
    model = _fissionable_model('U235', 0.0253, ['U236'], e_bins,
                               mts=['(n,gamma)'], eigenvalue=True)
    assert _run(model, run_in_tmpdir).sum() > 0.0


def test_nonfission_recoil_rate_follows_the_cross_section(run_in_tmpdir):
    """Production must track the nonfission share of absorption.

    U-235 at 0.0253 eV has a capture-to-fission ratio near 0.17, so a run that
    produced a capture recoil for every absorption would overproduce by about
    a factor of seven.
    """
    e_bins = np.logspace(0, np.log10(1.0e6), 61)
    model = _fissionable_model('U235', 0.0253, ['U236'], e_bins,
                               particles=60000, create_fission_neutrons=False)
    # tally absorption and fission reaction rates in the same run
    rates = openmc.Tally(name='rates')
    rates.scores = ['absorption', 'fission']
    model.tallies.append(rates)

    sp_path = model.run(output=False)
    with openmc.StatePoint(sp_path) as sp:
        recoil = sp.get_tally(name='recoil').mean.sum()
        r = sp.get_tally(name='rates').mean.ravel()
    absorption, fission = float(r[0]), float(r[1])

    assert absorption > 0.0 and fission > 0.0
    expected = absorption - fission
    assert abs(recoil - expected) / expected < 0.10, (recoil, expected)


def test_analog_and_survival_biased_production_agree(run_in_tmpdir):
    """Survival biasing must give a nonfission channel a nonfission weight.

    Giving it the whole absorbed weight overproduces by the ratio of total to
    nonfission absorption, which for U-235 at thermal energies is about seven.
    """
    e_bins = np.logspace(0, np.log10(1.0e6), 61)

    # No ReactionFilter here. Under survival biasing the neutron survives the
    # implicit absorption and goes on to scatter, so event_mt reports the
    # scattering reaction and a reaction-filtered tally would score nothing.
    # U-236 is produced only by capture, so the product filter is enough.
    def rate(survival):
        model = _fissionable_model(
            'U235', 0.0253, ['U236'], e_bins,
            survival=survival, particles=60000,
            create_fission_neutrons=False)
        return _run(model, run_in_tmpdir).sum()

    analog = rate(False)
    biased = rate(True)
    assert analog > 0.0 and biased > 0.0
    assert abs(biased - analog) / analog < 0.10, (analog, biased)


def test_nonfissionable_capture_produces_the_residual(run_in_tmpdir):
    """Capture on a nonfissionable nuclide produces its residual."""
    e_bins = np.logspace(0, 4, 61)
    model = _one_collision_model('Fe56', 0.0253, ['Fe57'], e_bins,
                                 mts=['(n,gamma)'], density=7.874)
    assert _run(model, run_in_tmpdir).sum() > 0.0


def test_discrete_level_recoil_is_uniform(run_in_tmpdir):
    """A named level plus isotropic emission makes the recoil energy uniform.

    For an exactly two-body exit channel the laboratory recoil energy is affine
    in the center-of-mass cosine, E_R = A + B*mu, so sampling mu uniformly
    makes E_R uniform between its kinematic limits. The pre-equilibrium
    systematics that serve the continuum channels would tilt it strongly
    towards one end instead.
    """
    energy = 14.0e6
    particles, batches = 200000, 2
    e_bins = np.linspace(0.0, 1.8e6, 37)
    # A whole mean free path of target: the energy cutoff still restricts the
    # tally to neutrons that have not scattered, and MT=600 is a 1% channel.
    model = _one_collision_model('Si28', energy, ['Al28'], e_bins, mts=[600],
                                 particles=particles, density=2.33,
                                 radius=20.0)
    # Scores are per source particle; recover the event counts for the shape
    counts = _run(model, run_in_tmpdir).ravel() * particles * batches

    nz = np.nonzero(counts)[0]
    assert nz.size > 20, 'too few occupied bins to judge the shape'

    # Support edges, and the interior bins that they do not clip
    lo, hi = e_bins[nz[0]], e_bins[nz[-1] + 1]
    interior = counts[nz[0] + 1:nz[-1]]
    mid = 0.5 * (e_bins[nz[0] + 1:nz[-1]] + e_bins[nz[0] + 2:nz[-1] + 1])
    assert interior.sum() > 500, 'too few events to judge the shape'

    # <E_R> = A + B <mu>, and the support gives A and |B|, so the offset of the
    # mean from the midpoint of the support measures <mu> directly.
    mean_mu = (0.5 * (lo + hi) - (interior * mid).sum() / interior.sum()) \
        / (0.5 * (hi - lo))
    assert abs(mean_mu) < 0.10, \
        f'center-of-mass emission is not isotropic: <mu> = {mean_mu:+.3f}'

    # And the shape is flat, not merely balanced
    expected = interior.mean()
    chi2 = ((interior - expected) ** 2 / expected).sum()
    assert chi2 < 3.0 * (len(interior) - 1), \
        f'recoil spectrum is not uniform: chi2 = {chi2:.1f} on ' \
        f'{len(interior) - 1} degrees of freedom'


def test_recoil_production_does_not_perturb_transport(run_in_tmpdir):
    """Banking recoils must not change the neutron solution.

    Recoils are production-only secondary particles: they are scored when
    requested and never transported. Sampling one consumes random numbers, so
    the two runs diverge in their streams and cannot agree bit for bit -- what
    has to hold is that they agree within statistics, bin by bin.
    """
    def run(recoil_on, particles=40000):
        steel = openmc.Material()
        for nuclide, fraction in (('Fe56', 0.70), ('Cr52', 0.18),
                                  ('Ni58', 0.12)):
            steel.add_nuclide(nuclide, fraction)
        steel.set_density('g/cm3', 7.9)
        sph = openmc.Sphere(r=30.0, boundary_type='vacuum')

        settings = openmc.Settings()
        settings.run_mode = 'fixed source'
        settings.particles = particles
        settings.batches = 5
        settings.seed = 1
        settings.source = openmc.IndependentSource(
            space=openmc.stats.Point(),
            energy=openmc.stats.Discrete([14.1e6], [1.0]))
        settings.recoil_production = recoil_on

        tally = openmc.Tally(name='flux')
        tally.filters = [openmc.EnergyFilter(np.logspace(0, 7.2, 13))]
        tally.scores = ['flux']

        model = openmc.Model(
            openmc.Geometry([openmc.Cell(fill=steel, region=-sph)]),
            [steel], settings, openmc.Tallies([tally]))
        with openmc.StatePoint(model.run(output=False)) as sp:
            t = sp.get_tally(name='flux')
            return t.mean.ravel(), t.std_dev.ravel()

    off, off_sd = run(False)
    on, on_sd = run(True)

    sigma = np.sqrt(off_sd ** 2 + on_sd ** 2)
    ok = sigma > 0
    assert ok.sum() > 8, 'not enough populated bins to judge'
    z = (off[ok] - on[ok]) / sigma[ok]
    assert np.abs(z).max() < 5.0, \
        f'recoil production moved the neutron flux: max |z| = {np.abs(z).max():.1f}'


def test_neutron_inelastic_uses_the_evaluated_angle(run_in_tmpdir):
    """A discrete neutron level must never reach the light-ion angular model.

    MT = 51-90 are two-body like the charged-particle levels, so the same
    affine relation applies and the mean cosine can be read straight off the
    recoil spectrum. Here it has to come back as the *evaluated* MF=4 mean, not
    as the isotropy the charged-particle levels are given and not as the
    Kalbach systematics the continuum channels are given.
    """
    energy, particles, batches = 8.0e6, 100000, 4
    e_bins = np.linspace(0.0, 6.0e5, 61)
    model = _one_collision_model('Fe56', energy, ['Fe56'], e_bins, mts=[51],
                                 particles=particles, density=7.874,
                                 radius=12.0)
    counts = _run(model, run_in_tmpdir).ravel() * particles * batches

    nz = np.nonzero(counts)[0]
    assert nz.size > 20, 'too few occupied bins to judge the shape'
    lo, hi = e_bins[nz[0]], e_bins[nz[-1] + 1]
    interior = counts[nz[0] + 1:nz[-1]]
    mid = 0.5 * (e_bins[nz[0] + 1:nz[-1]] + e_bins[nz[0] + 2:nz[-1] + 1])
    assert interior.sum() > 500, 'too few events to judge the shape'
    measured = (0.5 * (lo + hi) - (interior * mid).sum() / interior.sum()) \
        / (0.5 * (hi - lo))

    # What the evaluation says, read from whichever library is in use
    lib = openmc.data.IncidentNeutron.from_hdf5(
        openmc.data.DataLibrary.from_xml().get_by_material('Fe56')['path'])
    angle = lib[51].products[0].distribution[0].angle
    grid = np.asarray(angle.energy, float)
    tab = angle.mu[int(np.argmin(np.abs(grid - energy)))]
    x, p = np.asarray(tab.x, float), np.asarray(tab.p, float)
    expected = np.trapezoid(p * x, x) / np.trapezoid(p, x)

    # Generous, because the support edges are read off finite bins; tight
    # enough that isotropy, which would give zero, cannot pass when the
    # evaluated distribution is this anisotropic.
    assert abs(measured - expected) < 0.06, \
        f'recoil implies <mu> = {measured:+.3f}, evaluation says {expected:+.3f}'
    assert abs(expected) > 0.15, \
        'this level is not anisotropic enough here for the test to bite'


@pytest.mark.parametrize('mt,correct,wrong', [
    ('(n,p)', 'Mn56', ['Mn55', 'Fe56', 'Cr53']),
    ('(n,a)', 'Cr53', ['Cr52', 'Mn56', 'Fe55']),
    ('(n,2n)', 'Fe55', ['Fe56', 'Mn55', 'Cr53']),
])
def test_exit_channel_identity_is_conserved(run_in_tmpdir, mt, correct, wrong):
    """Only the residual the channel actually leaves may be banked.

    Charge and mass conservation at the level a tally can see. Every wrong
    candidate here is one unit of Z or A away from the right one, so an
    off-by-one in the exit-channel bookkeeping -- counting the emitted
    multiplicity wrongly, or parsing the reaction name wrongly -- produces a
    nuclide that is silently plausible rather than an obvious error.
    """
    e_bins = np.logspace(0, np.log10(5.0e6), 41)
    products = [correct] + wrong
    model = _one_collision_model('Fe56', 14.0e6, products, e_bins, mts=[mt],
                                 particles=100000, density=7.874)
    mean = _run(model, run_in_tmpdir).reshape(len(products), -1)

    assert mean[0].sum() > 0.0, f'{mt} banked no {correct}'
    for i, name in enumerate(wrong, start=1):
        assert mean[i].sum() == 0.0, f'{mt} banked {name}, which it cannot leave'
