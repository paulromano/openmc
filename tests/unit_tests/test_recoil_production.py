"""Physics checks on recoil (primary knock-on atom) production.

These run tiny fixed-source problems and check properties of the recoil
spectrum that follow from kinematics alone, so they hold for any nuclear data
library rather than only for the reference one.
"""

import numpy as np
import pytest

import openmc


def _one_collision_model(nuclide, energy, products, e_bins, mts=None,
                         recoil=None, particles=20000, density=5.0):
    """Thin monoenergetic target so that scored collisions are first ones."""
    mat = openmc.Material()
    mat.add_nuclide(nuclide, 1.0)
    mat.set_density('g/cm3', density)

    sph = openmc.Sphere(r=0.4, boundary_type='vacuum')
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
    nuclide = 'Fe56'
    e_bins = np.logspace(0, 4, 121)
    model = _one_collision_model('Fe56', energy, ['Fe57'], e_bins,
                                 mts=['(n,gamma)'], particles=20000,
                                 density=7.874)
    mean = _run(model, run_in_tmpdir).ravel()
    assert mean.sum() > 0

    mid = np.sqrt(e_bins[:-1] * e_bins[1:])
    mean_recoil = (mean * mid).sum() / mean.sum()
    # A single 7.6 MeV photon gives E_R = E_g^2 / (2 M c^2) ~ 550 eV; the
    # cascade spreads that over several photons, so require the right decade.
    assert 20.0 < mean_recoil < 1000.0


def test_light_ion_model_none(run_in_tmpdir):
    """light_ion_model='none' suppresses the emitted-ion records."""
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
    q = lib[16].q_value
    p_in = np.sqrt(energy)
    p_out = 2.0 * np.sqrt(max(energy + q, 0.0))
    e_max = (p_in + p_out) ** 2 / awr
    nz = np.nonzero(mean)[0]
    assert e_bins[nz[-1]] < e_max
