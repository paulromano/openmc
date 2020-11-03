from abc import ABC, abstractmethod
from collections.abc import Iterable
from functools import partial
from math import pi, sqrt, erf, exp, sinh
from numbers import Integral, Real
from warnings import warn

import numpy as np

import openmc.checkvalue as cv
from openmc.mixin import EqualityMixin
import openmc.stats.univariate as uni
from .data import EV_PER_MEV
from .endf import get_tab1_record, get_tab2_record
from .function import Tabulated1D, INTERPOLATION_SCHEME
from .grid import linearize


def _sample_with_max(sample_func, max_value, size):
    # Initial sampling
    rvs = sample_func(size)

    while True:
        # Check for values exceeding maximum energy
        idx, = (rvs > max_value).nonzero()
        if idx.size == 0:
            break

        # Replace values that were out of range
        rvs[idx] = sample_func(idx.size)
    return rvs


class EnergyDistribution(EqualityMixin, ABC):
    """Abstract superclass for all energy distributions."""
    def __init__(self):
        pass

    @abstractmethod
    def to_hdf5(self, group):
        pass

    @staticmethod
    def from_hdf5(group):
        """Generate energy distribution from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.EnergyDistribution
            Energy distribution

        """
        energy_type = group.attrs['type'].decode()
        if energy_type == 'maxwell':
            return MaxwellEnergy.from_hdf5(group)
        elif energy_type == 'evaporation':
            return Evaporation.from_hdf5(group)
        elif energy_type == 'watt':
            return WattEnergy.from_hdf5(group)
        elif energy_type == 'madland-nix':
            return MadlandNix.from_hdf5(group)
        elif energy_type == 'discrete_photon':
            return DiscretePhoton.from_hdf5(group)
        elif energy_type == 'level':
            return LevelInelastic.from_hdf5(group)
        elif energy_type == 'continuous':
            return ContinuousTabular.from_hdf5(group)
        else:
            raise ValueError("Unknown energy distribution type: {}"
                             .format(energy_type))

    @staticmethod
    def from_endf(file_obj, params):
        """Generate energy distribution from an ENDF evaluation

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of a section for an energy
            distribution.
        params : list
            List of parameters at the start of the energy distribution that
            includes the LF value indicating what type of energy distribution is
            present.

        Returns
        -------
        openmc.data.EnergyDistribution
            A sub-class of :class:`openmc.data.EnergyDistribution`

        """
        lf = params[3]
        if lf == 1:
            return ArbitraryTabulated.from_endf(file_obj, params)
        elif lf == 5:
            return GeneralEvaporation.from_endf(file_obj, params)
        elif lf == 7:
            return MaxwellEnergy.from_endf(file_obj, params)
        elif lf == 9:
            return Evaporation.from_endf(file_obj, params)
        elif lf == 11:
            return WattEnergy.from_endf(file_obj, params)
        elif lf == 12:
            return MadlandNix.from_endf(file_obj, params)

    @abstractmethod
    def sample(self, E, size=None):
        """Draw random samples for energy distribution

        Parameters
        ----------
        E : float
            Incident energy in [eV]
        size : int or tuple or ints, optional
            Output shape. If the given shape is, e.g., `(m, n)`, then `m*n`
            samples are drawn. If size is `None` (default), a single value is
            returned.

        Returns
        -------
        numpy.ndarray or scalar
            Sampled energies in [eV]

        """


class ArbitraryTabulated(EnergyDistribution):
    r"""Arbitrary tabulated function given in ENDF MF=5, LF=1 represented as

    .. math::
         f(E \rightarrow E') = g(E \rightarrow E')

    Parameters
    ----------
    energy : numpy.ndarray
        Array of incident neutron energies
    pdf : list of openmc.data.Tabulated1D
        Tabulated outgoing energy distribution probability density functions

    Attributes
    ----------
    energy : numpy.ndarray
        Array of incident neutron energies
    pdf : list of openmc.data.Tabulated1D
        Tabulated outgoing energy distribution probability density functions

    """

    def __init__(self, energy, pdf):
        super().__init__()
        self.energy = energy
        self.pdf = pdf

    def to_hdf5(self, group):
        raise NotImplementedError

    @classmethod
    def from_endf(cls, file_obj, params):
        """Generate arbitrary tabulated distribution from an ENDF evaluation

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of a section for an energy
            distribution.
        params : list
            List of parameters at the start of the energy distribution that
            includes the LF value indicating what type of energy distribution is
            present.

        Returns
        -------
        openmc.data.ArbitraryTabulated
            Arbitrary tabulated distribution

        """
        params, tab2 = get_tab2_record(file_obj)
        n_energies = params[5]

        energy = np.zeros(n_energies)
        pdf = []
        for j in range(n_energies):
            params, func = get_tab1_record(file_obj)
            energy[j] = params[1]
            pdf.append(func)
        return cls(energy, pdf)

    def sample(self, E, size=None):
        raise NotImplementedError


class GeneralEvaporation(EnergyDistribution):
    r"""General evaporation spectrum given in ENDF MF=5, LF=5 represented as

    .. math::
        f(E \rightarrow E') = g(E'/\theta(E))

    Parameters
    ----------
    theta : openmc.data.Tabulated1D
        Tabulated function of incident neutron energy :math:`E`
    g : openmc.data.Tabulated1D
        Tabulated function of :math:`x = E'/\theta(E)`
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    Attributes
    ----------
    theta : openmc.data.Tabulated1D
        Tabulated function of incident neutron energy :math:`E`
    g : openmc.data.Tabulated1D
        Tabulated function of :math:`x = E'/\theta(E)`
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    """

    def __init__(self, theta, g, u):
        super().__init__()
        self.theta = theta
        self.g = g
        self.u = u

    def to_hdf5(self, group):
        raise NotImplementedError

    @classmethod
    def from_ace(cls, ace, idx=0):
        raise NotImplementedError

    @classmethod
    def from_endf(cls, file_obj, params):
        """Generate general evaporation spectrum from an ENDF evaluation

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of a section for an energy
            distribution.
        params : list
            List of parameters at the start of the energy distribution that
            includes the LF value indicating what type of energy distribution is
            present.

        Returns
        -------
        openmc.data.GeneralEvaporation
            General evaporation spectrum

        """
        u = params[0]
        params, theta = get_tab1_record(file_obj)
        params, g = get_tab1_record(file_obj)
        return cls(theta, g, u)

    def sample(self, E, size=None):
        raise NotImplementedError


class MaxwellEnergy(EnergyDistribution):
    r"""Simple Maxwellian fission spectrum represented as

    .. math::
        f(E \rightarrow E') = \frac{\sqrt{E'}}{I} e^{-E'/\theta(E)}

    Parameters
    ----------
    theta : openmc.data.Tabulated1D
        Tabulated function of incident neutron energy
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    Attributes
    ----------
    theta : openmc.data.Tabulated1D
        Tabulated function of incident neutron energy
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    """

    def __init__(self, theta, u):
        super().__init__()
        self.theta = theta
        self.u = u

    @property
    def theta(self):
        return self._theta

    @property
    def u(self):
        return self._u

    @theta.setter
    def theta(self, theta):
        cv.check_type('Maxwell theta', theta, Tabulated1D)
        self._theta = theta

    @u.setter
    def u(self, u):
        cv.check_type('Maxwell restriction energy', u, Real)
        self._u = u

    def to_hdf5(self, group):
        """Write distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        group.attrs['type'] = np.string_('maxwell')
        group.attrs['u'] = self.u
        self.theta.to_hdf5(group, 'theta')

    @classmethod
    def from_hdf5(cls, group):
        """Generate Maxwell distribution from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.MaxwellEnergy
            Maxwell distribution

        """
        theta = Tabulated1D.from_hdf5(group['theta'])
        u = group.attrs['u']
        return cls(theta, u)

    @classmethod
    def from_ace(cls, ace, idx=0):
        """Create a Maxwell distribution from an ACE table

        Parameters
        ----------
        ace : openmc.data.ace.Table
            An ACE table
        idx : int
            Offset to read from in XSS array (default of zero)

        Returns
        -------
        openmc.data.MaxwellEnergy
            Maxwell distribution

        """
        # Read nuclear temperature -- since units are MeV, convert to eV
        theta = Tabulated1D.from_ace(ace, idx)
        theta.y *= EV_PER_MEV

        # Restriction energy
        nr = int(ace.xss[idx])
        ne = int(ace.xss[idx + 1 + 2*nr])
        u = ace.xss[idx + 2 + 2*nr + 2*ne]*EV_PER_MEV

        return cls(theta, u)

    @classmethod
    def from_endf(cls, file_obj, params):
        """Generate Maxwell distribution from an ENDF evaluation

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of a section for an energy
            distribution.
        params : list
            List of parameters at the start of the energy distribution that
            includes the LF value indicating what type of energy distribution is
            present.

        Returns
        -------
        openmc.data.MaxwellEnergy
            Maxwell distribution

        """
        u = params[0]
        params, theta = get_tab1_record(file_obj)
        return cls(theta, u)

    def sample(self, E, size=None):
        # Determine theta as function of incident energy
        theta = self.theta(E)

        # Sample outgoing energies
        dist = uni.Maxwell(theta)
        return _sample_with_max(dist.sample, E - self.u, size)

    def to_continuous_tabular(self):
        energy = self.theta.x
        breakpoints = [energy.size]
        interpolation = [2]

        # Define function for probability density
        def pdf(E):
            return sqrt(E)/norm * exp(-E / theta)

        energy_out = []
        u = self.u
        for E, theta in zip(self.theta.x, self.theta.y):
            # Determine normalization constant (I from ENDF-102 Section 5.1.1.3)
            v = (E - u)/theta
            norm = theta**(3/2)*(sqrt(pi)/2*erf(sqrt(v)) - sqrt(v)*exp(-v))

            # Convert density to tabulated form
            dist = uni.Tabular(*linearize([1.0e-5, E - u], pdf))
            energy_out.append(dist)

        return ContinuousTabular(breakpoints, interpolation, energy, energy_out)



class Evaporation(EnergyDistribution):
    r"""Evaporation spectrum represented as

    .. math::
        f(E \rightarrow E') = \frac{E'}{I} e^{-E'/\theta(E)}

    Parameters
    ----------
    theta : openmc.data.Tabulated1D
        Tabulated function of incident neutron energy
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    Attributes
    ----------
    theta : openmc.data.Tabulated1D
        Tabulated function of incident neutron energy
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    """

    def __init__(self, theta, u):
        super().__init__()
        self.theta = theta
        self.u = u

    @property
    def theta(self):
        return self._theta

    @property
    def u(self):
        return self._u

    @theta.setter
    def theta(self, theta):
        cv.check_type('Evaporation theta', theta, Tabulated1D)
        self._theta = theta

    @u.setter
    def u(self, u):
        cv.check_type('Evaporation restriction energy', u, Real)
        self._u = u

    def to_hdf5(self, group):
        """Write distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        group.attrs['type'] = np.string_('evaporation')
        group.attrs['u'] = self.u
        self.theta.to_hdf5(group, 'theta')

    @classmethod
    def from_hdf5(cls, group):
        """Generate evaporation spectrum from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.Evaporation
            Evaporation spectrum

        """
        theta = Tabulated1D.from_hdf5(group['theta'])
        u = group.attrs['u']
        return cls(theta, u)

    @classmethod
    def from_ace(cls, ace, idx=0):
        """Create an evaporation spectrum from an ACE table

        Parameters
        ----------
        ace : openmc.data.ace.Table
            An ACE table
        idx : int
            Offset to read from in XSS array (default of zero)

        Returns
        -------
        openmc.data.Evaporation
            Evaporation spectrum

        """
        # Read nuclear temperature -- since units are MeV, convert to eV
        theta = Tabulated1D.from_ace(ace, idx)
        theta.y *= EV_PER_MEV

        # Restriction energy
        nr = int(ace.xss[idx])
        ne = int(ace.xss[idx + 1 + 2*nr])
        u = ace.xss[idx + 2 + 2*nr + 2*ne]*EV_PER_MEV

        return cls(theta, u)

    @classmethod
    def from_endf(cls, file_obj, params):
        """Generate evaporation spectrum from an ENDF evaluation

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of a section for an energy
            distribution.
        params : list
            List of parameters at the start of the energy distribution that
            includes the LF value indicating what type of energy distribution is
            present.

        Returns
        -------
        openmc.data.Evaporation
            Evaporation spectrum

        """
        u = params[0]
        params, theta = get_tab1_record(file_obj)
        return cls(theta, u)

    def sample(self, E, size=None):
        # Determine theta as function of incident energy
        theta = self.theta(E)

        # Sample outgoing energies
        rg = np.random.default_rng()
        sample_func = partial(rg.gamma, 2, scale=theta)
        return _sample_with_max(sample_func, E - self.u, size)

    def to_continuous_tabular(self):
        energy = self.theta.x
        breakpoints = [energy.size]
        interpolation = [2]

        # Define function for probability density
        def pdf(E):
            return E/norm * exp(-E / theta)

        energy_out = []
        u = self.u
        for E, theta in zip(self.theta.x, self.theta.y):
            # If incident energy matches restriction energy, outgoing energy is zero
            if E == u:
                dist = uni.Discrete([0.0], [1.0])
            else:
                # Determine normalization constant (I from ENDF-102 Section 5.1.1.4)
                v = (E - u)/theta
                norm = theta**2*(1 - exp(-v)*(1 + v))

                # Convert density to tabulated form
                dist = uni.Tabular(*linearize([1.0e-5, E - u], pdf))

            energy_out.append(dist)

        return ContinuousTabular(breakpoints, interpolation, energy, energy_out)


class WattEnergy(EnergyDistribution):
    r"""Energy-dependent Watt spectrum represented as

    .. math::
        f(E \rightarrow E') = \frac{e^{-E'/a}}{I} \sinh \left ( \sqrt{bE'}
        \right )

    Parameters
    ----------
    a, b : openmc.data.Tabulated1D
        Energy-dependent parameters tabulated as function of incident neutron
        energy
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    Attributes
    ----------
    a, b : openmc.data.Tabulated1D
        Energy-dependent parameters tabulated as function of incident neutron
        energy
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    """

    def __init__(self, a, b, u):
        super().__init__()
        self.a = a
        self.b = b
        self.u = u

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

    @property
    def u(self):
        return self._u

    @a.setter
    def a(self, a):
        cv.check_type('Watt a', a, Tabulated1D)
        self._a = a

    @b.setter
    def b(self, b):
        cv.check_type('Watt b', b, Tabulated1D)
        self._b = b

    @u.setter
    def u(self, u):
        cv.check_type('Watt restriction energy', u, Real)
        self._u = u

    def to_hdf5(self, group):
        """Write distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        group.attrs['type'] = np.string_('watt')
        group.attrs['u'] = self.u
        self.a.to_hdf5(group, 'a')
        self.b.to_hdf5(group, 'b')

    @classmethod
    def from_hdf5(cls, group):
        """Generate Watt fission spectrum from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.WattEnergy
            Watt fission spectrum

        """
        a = Tabulated1D.from_hdf5(group['a'])
        b = Tabulated1D.from_hdf5(group['b'])
        u = group.attrs['u']
        return cls(a, b, u)

    @classmethod
    def from_ace(cls, ace, idx):
        """Create a Watt fission spectrum from an ACE table

        Parameters
        ----------
        ace : openmc.data.ace.Table
            An ACE table
        idx : int
            Offset to read from in XSS array (default of zero)

        Returns
        -------
        openmc.data.WattEnergy
            Watt fission spectrum

        """
        # Energy-dependent a parameter -- units are MeV, convert to eV
        a = Tabulated1D.from_ace(ace, idx)
        a.y *= EV_PER_MEV

        # Advance index
        nr = int(ace.xss[idx])
        ne = int(ace.xss[idx + 1 + 2*nr])
        idx += 2 + 2*nr + 2*ne

        # Energy-dependent b parameter -- units are MeV^-1
        b = Tabulated1D.from_ace(ace, idx)
        b.y /= EV_PER_MEV

        # Advance index
        nr = int(ace.xss[idx])
        ne = int(ace.xss[idx + 1 + 2*nr])
        idx += 2 + 2*nr + 2*ne

        # Restriction energy
        u = ace.xss[idx]*EV_PER_MEV

        return cls(a, b, u)

    @classmethod
    def from_endf(cls, file_obj, params):
        """Generate Watt fission spectrum from an ENDF evaluation

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of a section for an energy
            distribution.
        params : list
            List of parameters at the start of the energy distribution that
            includes the LF value indicating what type of energy distribution is
            present.

        Returns
        -------
        openmc.data.WattEnergy
            Watt fission spectrum

        """
        u = params[0]
        params, a = get_tab1_record(file_obj)
        params, b = get_tab1_record(file_obj)
        return cls(a, b, u)

    def sample(self, E, size=None):
        # Determine theta as function of incident energy
        a = self.a(E)
        b = self.b(E)

        # Sample outgoing energies
        dist = uni.Watt(a, b)
        return _sample_with_max(dist.sample, E - self.u, size)

    def to_continuous_tabular(self):
        energy = self.a.x
        breakpoints = [energy.size]
        interpolation = [2]

        # Define function for probability density
        def pdf(E):
            return exp(-E/a)/norm * sinh(sqrt(b*E))

        energy_out = []
        u = self.u
        for E, a, b in zip(self.a.x, self.a.y, self.b.y):
            # Determine normalization constant (I from ENDF-102 Section 5.1.1.5)
            v = (E - u)/a
            norm = 1/2*sqrt(pi*a**3*b/4)*exp(a*b/4)*(erf(sqrt(v) - sqrt(a*b/4)) \
                + erf(sqrt(v) + sqrt(a*b/4))) - a*exp(-v)*sinh(sqrt(b*(E - u)))

            # Convert density to tabulated form
            dist = uni.Tabular(*linearize([1.0e-5, E - u], pdf))
            energy_out.append(dist)

        return ContinuousTabular(breakpoints, interpolation, energy, energy_out)


class MadlandNix(EnergyDistribution):
    r"""Energy-dependent fission neutron spectrum (Madland and Nix) given in
    ENDF MF=5, LF=12 represented as

    .. math::
        f(E \rightarrow E') = \frac{1}{2} [ g(E', E_F(L)) + g(E', E_F(H))]

    where

    .. math::
        g(E',E_F) = \frac{1}{3\sqrt{E_F T_M}} \left [ u_2^{3/2} E_1 (u_2) -
        u_1^{3/2} E_1 (u_1) + \gamma \left ( \frac{3}{2}, u_2 \right ) - \gamma
        \left ( \frac{3}{2}, u_1 \right ) \right ] \\ u_1 = \left ( \sqrt{E'} -
        \sqrt{E_F} \right )^2 / T_M \\ u_2 = \left ( \sqrt{E'} + \sqrt{E_F}
        \right )^2 / T_M.

    Parameters
    ----------
    efl, efh : float
        Constants which represent the average kinetic energy per nucleon of the
        fission fragment (efl = light, efh = heavy)
    tm : openmc.data.Tabulated1D
        Parameter tabulated as a function of incident neutron energy

    Attributes
    ----------
    efl, efh : float
        Constants which represent the average kinetic energy per nucleon of the
        fission fragment (efl = light, efh = heavy)
    tm : openmc.data.Tabulated1D
        Parameter tabulated as a function of incident neutron energy

    """

    def __init__(self, efl, efh, tm):
        super().__init__()
        self.efl = efl
        self.efh = efh
        self.tm = tm

    @property
    def efl(self):
        return self._efl

    @property
    def efh(self):
        return self._efh

    @property
    def tm(self):
        return self._tm

    @efl.setter
    def efl(self, efl):
        name = 'Madland-Nix light fragment energy'
        cv.check_type(name, efl, Real)
        cv.check_greater_than(name, efl, 0.)
        self._efl = efl

    @efh.setter
    def efh(self, efh):
        name = 'Madland-Nix heavy fragment energy'
        cv.check_type(name, efh, Real)
        cv.check_greater_than(name, efh, 0.)
        self._efh = efh

    @tm.setter
    def tm(self, tm):
        cv.check_type('Madland-Nix maximum temperature', tm, Tabulated1D)
        self._tm = tm

    def to_hdf5(self, group):
        """Write distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        group.attrs['type'] = np.string_('madland-nix')
        group.attrs['efl'] = self.efl
        group.attrs['efh'] = self.efh
        self.tm.to_hdf5(group)

    @classmethod
    def from_hdf5(cls, group):
        """Generate Madland-Nix fission spectrum from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.MadlandNix
            Madland-Nix fission spectrum

        """
        efl = group.attrs['efl']
        efh = group.attrs['efh']
        tm = Tabulated1D.from_hdf5(group['tm'])
        return cls(efl, efh, tm)

    @classmethod
    def from_endf(cls, file_obj, params):
        """Generate Madland-Nix fission spectrum from an ENDF evaluation

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of a section for an energy
            distribution.
        params : list
            List of parameters at the start of the energy distribution that
            includes the LF value indicating what type of energy distribution is
            present.

        Returns
        -------
        openmc.data.MadlandNix
            Madland-Nix fission spectrum

        """
        params, tm = get_tab1_record(file_obj)
        efl, efh = params[0:2]
        return cls(efl, efh, tm)



class DiscretePhoton(EnergyDistribution):
    """Discrete photon energy distribution

    Parameters
    ----------
    primary_flag : int
        Indicator of whether the photon is a primary or non-primary photon.
    energy : float
        Photon energy (if lp==0 or lp==1) or binding energy (if lp==2).
    atomic_weight_ratio : float
        Atomic weight ratio of the target nuclide responsible for the emitted
        particle

    Attributes
    ----------
    primary_flag : int
        Indicator of whether the photon is a primary or non-primary photon.
    energy : float
        Photon energy (if lp==0 or lp==1) or binding energy (if lp==2).
    atomic_weight_ratio : float
        Atomic weight ratio of the target nuclide responsible for the emitted
        particle

    """

    def __init__(self, primary_flag, energy, atomic_weight_ratio):
        super().__init__()
        self.primary_flag = primary_flag
        self.energy = energy
        self.atomic_weight_ratio = atomic_weight_ratio

    @property
    def primary_flag(self):
        return self._primary_flag

    @property
    def energy(self):
        return self._energy

    @property
    def atomic_weight_ratio(self):
        return self._atomic_weight_ratio

    @primary_flag.setter
    def primary_flag(self, primary_flag):
        cv.check_type('discrete photon primary_flag', primary_flag, Integral)
        self._primary_flag = primary_flag

    @energy.setter
    def energy(self, energy):
        cv.check_type('discrete photon energy', energy, Real)
        self._energy = energy

    @atomic_weight_ratio.setter
    def atomic_weight_ratio(self, atomic_weight_ratio):
        cv.check_type('atomic weight ratio', atomic_weight_ratio, Real)
        self._atomic_weight_ratio = atomic_weight_ratio

    def to_hdf5(self, group):
        """Write distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        group.attrs['type'] = np.string_('discrete_photon')
        group.attrs['primary_flag'] = self.primary_flag
        group.attrs['energy'] = self.energy
        group.attrs['atomic_weight_ratio'] = self.atomic_weight_ratio

    @classmethod
    def from_hdf5(cls, group):
        """Generate discrete photon energy distribution from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.DiscretePhoton
            Discrete photon energy distribution

        """
        primary_flag = group.attrs['primary_flag']
        energy = group.attrs['energy']
        awr = group.attrs['atomic_weight_ratio']
        return cls(primary_flag, energy, awr)

    @classmethod
    def from_ace(cls, ace, idx):
        """Generate discrete photon energy distribution from an ACE table

        Parameters
        ----------
        ace : openmc.data.ace.Table
            An ACE table
        idx : int
            Offset to read from in XSS array (default of zero)

        Returns
        -------
        openmc.data.DiscretePhoton
            Discrete photon energy distribution

        """
        primary_flag = int(ace.xss[idx])
        energy = ace.xss[idx + 1]*EV_PER_MEV
        return cls(primary_flag, energy, ace.atomic_weight_ratio)

    def sample(self, E, size=None):
        if self.primary_flag == 2:
            A = self.atomic_weight_ratio
            rvs = np.full(size, self.energy + (A + 1)/A * E)
        else:
            rvs = np.full(size, self.energy)
        return rvs.item() if rvs.shape == () else rvs

    def to_continuous_tabular(self):
        energy = [0.0, 150e6]
        breakpoints = [2]
        interpolation = [2]

        energy_out = []
        for E in energy:
            # Get outgoing energy
            E_out = self.sample(E)

            # Create discrete distribution
            energy_out.append(uni.Discrete([E_out], [1.0]))

        return ContinuousTabular(breakpoints, interpolation, energy, energy_out)


class LevelInelastic(EnergyDistribution):
    r"""Level inelastic scattering

    Parameters
    ----------
    threshold : float
        Energy threshold in the laboratory system, :math:`(A + 1)/A * |Q|`
    mass_ratio : float
        :math:`(A/(A + 1))^2`

    Attributes
    ----------
    threshold : float
        Energy threshold in the laboratory system, :math:`(A + 1)/A * |Q|`
    mass_ratio : float
        :math:`(A/(A + 1))^2`

    """

    def __init__(self, threshold, mass_ratio):
        super().__init__()
        self.threshold = threshold
        self.mass_ratio = mass_ratio

    @property
    def threshold(self):
        return self._threshold

    @property
    def mass_ratio(self):
        return self._mass_ratio

    @threshold.setter
    def threshold(self, threshold):
        cv.check_type('level inelastic threhsold', threshold, Real)
        self._threshold = threshold

    @mass_ratio.setter
    def mass_ratio(self, mass_ratio):
        cv.check_type('level inelastic mass ratio', mass_ratio, Real)
        self._mass_ratio = mass_ratio

    def to_hdf5(self, group):
        """Write distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        group.attrs['type'] = np.string_('level')
        group.attrs['threshold'] = self.threshold
        group.attrs['mass_ratio'] = self.mass_ratio

    @classmethod
    def from_hdf5(cls, group):
        """Generate level inelastic distribution from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.LevelInelastic
            Level inelastic scattering distribution

        """
        threshold = group.attrs['threshold']
        mass_ratio = group.attrs['mass_ratio']
        return cls(threshold, mass_ratio)

    @classmethod
    def from_ace(cls, ace, idx):
        """Generate level inelastic distribution from an ACE table

        Parameters
        ----------
        ace : openmc.data.ace.Table
            An ACE table
        idx : int
            Offset to read from in XSS array (default of zero)

        Returns
        -------
        openmc.data.LevelInelastic
            Level inelastic scattering distribution

        """
        threshold = ace.xss[idx]*EV_PER_MEV
        mass_ratio = ace.xss[idx + 1]
        return cls(threshold, mass_ratio)

    def sample(self, E, size=None):
        E_out = self.mass_ratio*(E - self.threshold)
        rvs = np.full(size, E_out)
        return rvs.item() if rvs.shape == () else rvs

    to_continuous_tabular = DiscretePhoton.to_continuous_tabular


class ContinuousTabular(EnergyDistribution):
    """Continuous tabular distribution

    Parameters
    ----------
    breakpoints : Iterable of int
        Breakpoints defining interpolation regions
    interpolation : Iterable of int
        Interpolation codes
    energy : Iterable of float
        Incoming energies at which distributions exist
    energy_out : Iterable of openmc.stats.Univariate
        Distribution of outgoing energies corresponding to each incoming energy

    Attributes
    ----------
    breakpoints : Iterable of int
        Breakpoints defining interpolation regions
    interpolation : Iterable of int
        Interpolation codes
    energy : Iterable of float
        Incoming energies at which distributions exist
    energy_out : Iterable of openmc.stats.Univariate
        Distribution of outgoing energies corresponding to each incoming energy

    """

    def __init__(self, breakpoints, interpolation, energy, energy_out):
        super().__init__()
        self.breakpoints = breakpoints
        self.interpolation = interpolation
        self.energy = energy
        self.energy_out = energy_out

    @property
    def breakpoints(self):
        return self._breakpoints

    @property
    def interpolation(self):
        return self._interpolation

    @property
    def energy(self):
        return self._energy

    @property
    def energy_out(self):
        return self._energy_out

    @breakpoints.setter
    def breakpoints(self, breakpoints):
        cv.check_type('continuous tabular breakpoints', breakpoints,
                      Iterable, Integral)
        self._breakpoints = breakpoints

    @interpolation.setter
    def interpolation(self, interpolation):
        cv.check_type('continuous tabular interpolation', interpolation,
                      Iterable, Integral)
        self._interpolation = interpolation

    @energy.setter
    def energy(self, energy):
        cv.check_type('continuous tabular incoming energy', energy,
                      Iterable, Real)
        self._energy = energy

    @energy_out.setter
    def energy_out(self, energy_out):
        cv.check_type('continuous tabular outgoing energy', energy_out,
                      Iterable, uni.Univariate)
        self._energy_out = energy_out

    def to_hdf5(self, group):
        """Write distribution to an HDF5 group

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to write to

        """

        group.attrs['type'] = np.string_('continuous')

        dset = group.create_dataset('energy', data=self.energy)
        dset.attrs['interpolation'] = np.vstack((self.breakpoints,
                                                 self.interpolation))

        # Determine total number of (E,p) pairs and create array
        n_pairs = sum(len(d) for d in self.energy_out)
        pairs = np.empty((3, n_pairs))

        # Create array for offsets
        offsets = np.empty(len(self.energy_out), dtype=int)
        interpolation = np.empty(len(self.energy_out), dtype=int)
        n_discrete_lines = np.empty(len(self.energy_out), dtype=int)
        j = 0

        # Populate offsets and pairs array
        for i, eout in enumerate(self.energy_out):
            n = len(eout)
            offsets[i] = j

            if isinstance(eout, uni.Mixture):
                discrete, continuous = eout.distribution
                n_discrete_lines[i] = m = len(discrete)
                interpolation[i] = 1 if continuous.interpolation == 'histogram' else 2
                pairs[0, j:j+m] = discrete.x
                pairs[1, j:j+m] = discrete.p
                pairs[2, j:j+m] = discrete.c
                pairs[0, j+m:j+n] = continuous.x
                pairs[1, j+m:j+n] = continuous.p
                pairs[2, j+m:j+n] = continuous.c
            else:
                if isinstance(eout, uni.Tabular):
                    n_discrete_lines[i] = 0
                    interpolation[i] = 1 if eout.interpolation == 'histogram' else 2
                elif isinstance(eout, uni.Discrete):
                    n_discrete_lines[i] = n
                    interpolation[i] = 1
                pairs[0, j:j+n] = eout.x
                pairs[1, j:j+n] = eout.p
                pairs[2, j:j+n] = eout.c
            j += n

        # Create dataset for distributions
        dset = group.create_dataset('distribution', data=pairs)

        # Write interpolation as attribute
        dset.attrs['offsets'] = offsets
        dset.attrs['interpolation'] = interpolation
        dset.attrs['n_discrete_lines'] = n_discrete_lines

    @classmethod
    def from_hdf5(cls, group):
        """Generate continuous tabular distribution from HDF5 data

        Parameters
        ----------
        group : h5py.Group
            HDF5 group to read from

        Returns
        -------
        openmc.data.ContinuousTabular
            Continuous tabular energy distribution

        """
        interp_data = group['energy'].attrs['interpolation']
        energy_breakpoints = interp_data[0, :]
        energy_interpolation = interp_data[1, :]
        energy = group['energy'][()]

        data = group['distribution']
        offsets = data.attrs['offsets']
        interpolation = data.attrs['interpolation']
        n_discrete_lines = data.attrs['n_discrete_lines']

        energy_out = []
        n_energy = len(energy)
        for i in range(n_energy):
            # Determine length of outgoing energy distribution and number of
            # discrete lines
            j = offsets[i]
            if i < n_energy - 1:
                n = offsets[i+1] - j
            else:
                n = data.shape[1] - j
            m = n_discrete_lines[i]

            # Create discrete distribution if lines are present
            if m > 0:
                eout_discrete = uni.Discrete(data[0, j:j+m], data[1, j:j+m])
                eout_discrete.c = data[2, j:j+m]
                p_discrete = eout_discrete.c[-1]

            # Create continuous distribution
            if m < n:
                interp = INTERPOLATION_SCHEME[interpolation[i]]
                eout_continuous = uni.Tabular(data[0, j+m:j+n], data[1, j+m:j+n], interp)
                eout_continuous.c = data[2, j+m:j+n]

            # If both continuous and discrete are present, create a mixture
            # distribution
            if m == 0:
                eout_i = eout_continuous
            elif m == n:
                eout_i = eout_discrete
            else:
                eout_i = uni.Mixture([p_discrete, 1. - p_discrete],
                                     [eout_discrete, eout_continuous])
            energy_out.append(eout_i)

        return cls(energy_breakpoints, energy_interpolation,
                   energy, energy_out)

    @classmethod
    def from_ace(cls, ace, idx, ldis):
        """Generate continuous tabular energy distribution from ACE data

        Parameters
        ----------
        ace : openmc.data.ace.Table
            ACE table to read from
        idx : int
            Index in XSS array of the start of the energy distribution data
            (LDIS + LOCC - 1)
        ldis : int
            Index in XSS array of the start of the energy distribution block
            (e.g. JXS[11])

        Returns
        -------
        openmc.data.ContinuousTabular
            Continuous tabular energy distribution

        """
        # Read number of interpolation regions and incoming energies
        n_regions = int(ace.xss[idx])
        n_energy_in = int(ace.xss[idx + 1 + 2*n_regions])

        # Get interpolation information
        idx += 1
        if n_regions > 0:
            breakpoints = ace.xss[idx:idx + n_regions].astype(int)
            interpolation = ace.xss[idx + n_regions:idx + 2*n_regions].astype(int)
        else:
            breakpoints = np.array([n_energy_in])
            interpolation = np.array([2])

        # Incoming energies at which distributions exist
        idx += 2*n_regions + 1
        energy = ace.xss[idx:idx + n_energy_in]*EV_PER_MEV

        # Location of distributions
        idx += n_energy_in
        loc_dist = ace.xss[idx:idx + n_energy_in].astype(int)

        # Initialize variables
        energy_out = []

        # Read each outgoing energy distribution
        for i in range(n_energy_in):
            idx = ldis + loc_dist[i] - 1

            # intt = interpolation scheme (1=hist, 2=lin-lin)
            INTTp = int(ace.xss[idx])
            intt = INTTp % 10
            n_discrete_lines = (INTTp - intt)//10
            if intt not in (1, 2):
                warn("Interpolation scheme for continuous tabular distribution "
                     "is not histogram or linear-linear.")
                intt = 2

            n_energy_out = int(ace.xss[idx + 1])
            data = ace.xss[idx + 2:idx + 2 + 3*n_energy_out].copy()
            data.shape = (3, n_energy_out)
            data[0,:] *= EV_PER_MEV

            # Create continuous distribution
            eout_continuous = uni.Tabular(data[0][n_discrete_lines:],
                                      data[1][n_discrete_lines:]/EV_PER_MEV,
                                      INTERPOLATION_SCHEME[intt])
            eout_continuous.c = data[2][n_discrete_lines:]

            # If discrete lines are present, create a mixture distribution
            if n_discrete_lines > 0:
                eout_discrete = uni.Discrete(data[0][:n_discrete_lines],
                                         data[1][:n_discrete_lines])
                eout_discrete.c = data[2][:n_discrete_lines]
                if n_discrete_lines == n_energy_out:
                    eout_i = eout_discrete
                else:
                    p_discrete = min(sum(eout_discrete.p), 1.0)
                    eout_i = uni.Mixture([p_discrete, 1. - p_discrete],
                                     [eout_discrete, eout_continuous])
            else:
                eout_i = eout_continuous

            energy_out.append(eout_i)

        return cls(breakpoints, interpolation, energy, energy_out)

    def sample(self, E, size=None):
        i = np.searchsorted(self.energy, E, side='right') - 1
        r = (E - self.energy[i]) / (self.energy[i + 1] - self.energy[i])

        histogram_interp = (self.interpolation[0] == 1)
        rg = np.random.default_rng()
        if histogram_interp:
            l = i
        else:
            l = i + 1 if r > rg.uniform() else i

        # TODO: Account for discrete distributions
        E_i_1 = self.energy_out[i].x[0]
        E_i_K = self.energy_out[i].x[-1]
        E_i1_1 = self.energy_out[i+1].x[0]
        E_i1_K = self.energy_out[i+1].x[-1]
        E_1 = E_i_1 + r * (E_i1_1 - E_i_1)
        E_K = E_i_K + r * (E_i1_K - E_i_K)

        E_out = self.energy_out[l].sample(size=size)

        if not histogram_interp:
            if l == i:
                return E_1 + (E_out - E_i_1)*(E_K - E_1)/(E_i_K - E_i_1)
            else:
                return E_1 + (E_out - E_i1_1)*(E_K - E_1)/(E_i1_K - E_i1_1)
        else:
            return E_out
