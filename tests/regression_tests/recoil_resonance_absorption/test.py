import numpy as np
import openmc
import pytest
from openmc.utility_funcs import change_directory

from tests.testing_harness import PyAPITestHarness


FISSION_MTS = [18, 19, 20, 21, 38]


def make_model(energy, multipole, survival):
    model = openmc.Model()

    uranium = openmc.Material()
    uranium.add_nuclide('U235', 1.0)
    uranium.set_density('g/cm3', 18.0)
    model.materials.append(uranium)

    sphere = openmc.Sphere(r=0.4, boundary_type='vacuum')
    model.geometry = openmc.Geometry([
        openmc.Cell(fill=uranium, region=-sphere)
    ])

    model.settings.run_mode = 'fixed source'
    model.settings.batches = 10
    model.settings.particles = 10000
    model.settings.create_fission_neutrons = False
    model.settings.survival_biasing = survival
    model.settings.recoil_production = True
    model.settings.ptables = True
    model.settings.temperature = {'multipole': multipole}
    model.settings.cutoff = {'energy_neutron': 0.999 * energy}
    model.settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([energy], [1.0]),
    )

    rates = openmc.Tally(name='rates')
    rates.nuclides = ['U235']
    rates.scores = ['(n,gamma)', 'absorption', 'fission']

    production = openmc.Tally(name='production')
    production.filters = [
        openmc.ParticleProductionFilter(['U236'], [0.0, 1.0e7])
    ]
    if not survival:
        production.filters.insert(
            0, openmc.ReactionFilter([102, *FISSION_MTS]))
    production.scores = ['events']

    model.tallies = [rates, production]
    return model


class ResonanceAbsorptionHarness(PyAPITestHarness):
    def __init__(self, model, survival):
        super().__init__('statepoint.10.h5', model)
        self._survival = survival

    def _get_results(self):
        results = super()._get_results()
        with openmc.StatePoint(self._sp_name) as sp:
            rates = sp.get_tally(name='rates')
            ngamma, absorption, fission = rates.mean.ravel()
            ngamma_sd, _, _ = rates.std_dev.ravel()

            production = sp.get_tally(name='production')
            produced = production.mean.ravel()
            produced_sd = production.std_dev.ravel()

        assert ngamma > 0.0
        assert fission > 0.0
        assert ngamma == pytest.approx(absorption - fission, rel=1.0e-12)

        if self._survival:
            capture_production = produced.sum()
            capture_production_sd = np.sqrt((produced_sd**2).sum())
        else:
            # MT=102 is first in the ReactionFilter. A fission absorption must
            # not be relabelled as capture and produce U236.
            capture_production = produced[0]
            capture_production_sd = produced_sd[0]
            assert np.all(produced[1:] == 0.0)

        difference = abs(capture_production - ngamma)
        combined_sd = np.hypot(capture_production_sd, ngamma_sd)
        assert difference <= max(5.0 * combined_sd, 0.02 * ngamma), (
            capture_production, ngamma, combined_sd)
        return results


@pytest.mark.parametrize('energy,multipole,survival,subdir', [
    (1.0e3, True, False, 'wmp_analog'),
    (1.0e3, True, True, 'wmp_survival'),
    (1.0e4, False, False, 'urr_analog'),
    (1.0e4, False, True, 'urr_survival'),
    (1.0e4, False, False, 'pointwise_control'),
])
def test_recoil_resonance_absorption(energy, multipole, survival, subdir):
    with change_directory(subdir):
        openmc.reset_auto_ids()
        model = make_model(energy, multipole, survival)
        if subdir == 'pointwise_control':
            model.settings.ptables = False
        harness = ResonanceAbsorptionHarness(model, survival)
        harness.main()
