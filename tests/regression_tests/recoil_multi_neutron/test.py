import numpy as np
import openmc

from tests.testing_harness import PyAPITestHarness


def make_model():
    material = openmc.Material()
    material.add_nuclide('Fe56', 1.0)
    material.set_density('g/cm3', 7.874)

    sphere = openmc.Sphere(r=1.0, boundary_type='vacuum')
    geometry = openmc.Geometry([
        openmc.Cell(fill=material, region=-sphere)
    ])

    settings = openmc.Settings()
    settings.run_mode = 'fixed source'
    settings.batches = 10
    settings.particles = 100000
    settings.recoil_production = True
    settings.cutoff = {'energy_neutron': 0.999 * 14.0e6}
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([14.0e6], [1.0]),
    )

    rate = openmc.Tally(name='rate')
    rate.nuclides = ['Fe56']
    rate.scores = ['(n,2n)']

    production = openmc.Tally(name='production')
    production.filters = [
        openmc.ReactionFilter([16]),
        openmc.ParticleProductionFilter(
            ['Fe55'], np.linspace(0.0, 1.0e6, 41)),
    ]
    production.scores = ['events']

    return openmc.Model(geometry, [material], settings,
                        openmc.Tallies([rate, production]))


class MultiNeutronHarness(PyAPITestHarness):
    def _get_results(self):
        results = super()._get_results()
        with openmc.StatePoint(self._sp_name) as sp:
            rate_tally = sp.get_tally(name='rate')
            rate = rate_tally.mean.sum()
            rate_sd = np.sqrt((rate_tally.std_dev**2).sum())

            production = sp.get_tally(name='production')
            produced = production.mean.sum()
            produced_sd = np.sqrt((production.std_dev**2).sum())

        assert rate > 0.0
        assert produced > 0.0
        loss = rate - produced
        combined_sd = np.hypot(rate_sd, produced_sd)
        assert loss >= -5.0 * combined_sd
        # Some evaluated first-neutron marginal draws leave no feasible
        # remainder. They cannot be repaired without changing neutron
        # transport, so pin the measured loss instead of normalizing it away.
        assert loss <= max(5.0 * combined_sd, 0.02 * rate), (
            rate, produced, combined_sd)
        return results


def test_recoil_multi_neutron():
    openmc.reset_auto_ids()
    harness = MultiNeutronHarness('statepoint.10.h5', make_model())
    harness.main()
