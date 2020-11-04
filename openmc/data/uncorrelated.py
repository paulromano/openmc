import numpy as np

import openmc.checkvalue as cv
from openmc.stats.univariate import interpolate_tabular, Tabular
from .angle_energy import AngleEnergy
from .energy_distribution import EnergyDistribution
from .angle_distribution import AngleDistribution
from .correlated import CorrelatedAngleEnergy


class UncorrelatedAngleEnergy(AngleEnergy):
    """Uncorrelated angle-energy distribution

    Parameters
    ----------
    angle : openmc.data.AngleDistribution
        Distribution of outgoing angles represented as scattering cosines
    energy : openmc.data.EnergyDistribution
        Distribution of outgoing energies

    Attributes
    ----------
    angle : openmc.data.AngleDistribution
        Distribution of outgoing angles represented as scattering cosines
    energy : openmc.data.EnergyDistribution
        Distribution of outgoing energies

    """

    def __init__(self, angle=None, energy=None):
        self._angle = None
        self._energy = None

        if angle is not None:
            self.angle = angle
        if energy is not None:
            self.energy = energy

    @property
    def angle(self):
        return self._angle

    @property
    def energy(self):
        return self._energy

    @angle.setter
    def angle(self, angle):
        cv.check_type('uncorrelated angle distribution', angle,
                      AngleDistribution)
        self._angle = angle

    @energy.setter
    def energy(self, energy):
        cv.check_type('uncorrelated energy distribution', energy,
                      EnergyDistribution)
        self._energy = energy

    def to_hdf5(self, group):
        """Write distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """
        group.attrs['type'] = np.string_('uncorrelated')
        if self.angle is not None:
            angle_group = group.create_group('angle')
            self.angle.to_hdf5(angle_group)

        if self.energy is not None:
            energy_group = group.create_group('energy')
            self.energy.to_hdf5(energy_group)

    @classmethod
    def from_hdf5(cls, group):
        """Generate uncorrelated angle-energy distribution from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.UncorrelatedAngleEnergy
            Uncorrelated angle-energy distribution

        """
        dist = cls()
        if 'angle' in group:
            dist.angle = AngleDistribution.from_hdf5(group['angle'])
        if 'energy' in group:
            dist.energy = EnergyDistribution.from_hdf5(group['energy'])
        return dist

    def to_correlated(self):
        # Need breakpoints, interpolation, energy, energy_out, and mu
        energy_dist = self.energy.to_continuous_tabular()

        energy = energy_dist.energy
        energy_out = energy_dist.energy_out

        # Get angle distribution in tabular form
        mu = []
        for ein_i, eout_i in zip(energy, energy_out):
            if self.angle is not None:
                # Determine correct mu distribution to use
                ein = self.angle.energy
                if ein_i >= ein[-1]:
                    idx = len(ein) - 2
                else:
                    idx = np.searchsorted(ein, ein_i, 'right') - 1

                # Interpolate tabular mu distributions
                f = (ein_i - ein[idx]) / (ein[idx + 1] - ein[idx])
                f = max(0.0, min(1.0, f))
                mu_i = interpolate_tabular(
                    self.angle.mu[idx], self.angle.mu[idx + 1], f)

            else:
                # If not angle distribution specified, it is isotropic
                mu_i = Tabular([-1., 1.], [0.5, 0.5])

            mu.append([mu_i]*len(eout_i))

        breakpoints = [len(energy)]
        interpolation = [2]
        return CorrelatedAngleEnergy(breakpoints, interpolation, energy, energy_out, mu)
