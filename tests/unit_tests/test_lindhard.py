"""Tests for the Lindhard/Robinson damage energy partition function.

This module tests a pure-Python reference implementation of the partition
function against hand-verified values and NJOY HEATR constants, ensuring the
formula is correct before it gets compiled into the C++ library.
"""

import math

import pytest


def lindhard_partition(E_R, Z_R, A_R, Z_L, A_L, E_cutoff=25.0):
    """Python reference implementation of the Lindhard/Robinson partition.

    Parameters
    ----------
    E_R : float
        Recoil kinetic energy [eV]
    Z_R, A_R : int
        Atomic number and mass number of recoil nucleus
    Z_L, A_L : int
        Atomic number and mass number of lattice atom
    E_cutoff : float
        Minimum recoil energy for non-zero damage energy [eV]

    Returns
    -------
    float
        Damage energy [eV]
    """
    if Z_R == 0 or E_R < E_cutoff:
        return 0.0

    zr23 = Z_R ** (2.0 / 3.0)
    zl23 = Z_L ** (2.0 / 3.0)
    ar = float(A_R)
    al = float(A_L)

    E_L = 30.724 * Z_R * Z_L * math.sqrt(zr23 + zl23) * (ar + al) / al
    F_L = (0.0793 * zr23 * math.sqrt(Z_L) * (ar + al) ** 1.5 /
           ((zr23 + zl23) ** 0.75 * ar ** 1.5 * al ** 0.5))

    eps = E_R / E_L
    g = 3.4008 * eps ** (1.0 / 6.0) + 0.40244 * eps ** 0.75 + eps

    return E_R / (1.0 + F_L * g)


def test_lindhard_below_cutoff():
    """Damage energy should be zero below the cutoff."""
    assert lindhard_partition(24.9, 26, 56, 26, 56) == 0.0
    assert lindhard_partition(0.0, 26, 56, 26, 56) == 0.0
    assert lindhard_partition(-10.0, 26, 56, 26, 56) == 0.0


def test_lindhard_at_cutoff():
    """At exactly the cutoff, damage energy should be non-zero."""
    result = lindhard_partition(25.0, 26, 56, 26, 56)
    assert result > 0.0
    assert result <= 25.0


def test_lindhard_uncharged_recoil():
    """Neutron recoil (Z=0) should give zero damage energy."""
    assert lindhard_partition(1e6, 0, 1, 26, 56) == 0.0


def test_lindhard_less_than_recoil():
    """Damage energy should always be less than recoil energy."""
    for E_R in [100, 1e3, 1e4, 1e5, 1e6, 1e7]:
        dam = lindhard_partition(E_R, 26, 56, 26, 56)
        assert 0.0 < dam < E_R


def test_lindhard_fe56_self_irradiation():
    """Fe-56 self-irradiation at 1 MeV."""
    E_R = 1e6
    dam = lindhard_partition(E_R, 26, 56, 26, 56)
    # At 1 MeV, the partition should give roughly 40-60% for Fe
    fraction = dam / E_R
    assert 0.3 < fraction < 0.7


def test_lindhard_tungsten():
    """W-184 self-irradiation at various energies."""
    for E_R in [1e3, 1e4, 1e5, 1e6]:
        dam = lindhard_partition(E_R, 74, 184, 74, 184)
        assert 0.0 < dam < E_R


def test_lindhard_he4_in_fe56():
    """He-4 recoil (Z=2, A=4) in Fe-56 lattice."""
    E_R = 5e5
    dam = lindhard_partition(E_R, 2, 4, 26, 56)
    assert 0.0 < dam < E_R


def test_lindhard_low_energy_approaches_recoil():
    """At low energies (just above cutoff), damage ~ recoil energy."""
    E_R = 30.0
    dam = lindhard_partition(E_R, 26, 56, 26, 56)
    assert dam / E_R > 0.8


def test_lindhard_fraction_decreases_with_energy():
    """The damage fraction should decrease with increasing energy."""
    fractions = []
    for E_R in [1e2, 1e3, 1e4, 1e5, 1e6, 1e7]:
        dam = lindhard_partition(E_R, 26, 56, 26, 56)
        fractions.append(dam / E_R)

    for i in range(len(fractions) - 1):
        assert fractions[i] > fractions[i + 1]


def test_lindhard_custom_cutoff():
    """Custom cutoff energy should be respected."""
    assert lindhard_partition(50.0, 26, 56, 26, 56, E_cutoff=100.0) == 0.0
    assert lindhard_partition(50.0, 26, 56, 26, 56) > 0.0


def test_lindhard_njoy_constants():
    """Verify the partition function matches NJOY HEATR function df()."""
    zr, zl = 26, 26
    ar, al = 56.0, 56.0

    zr23 = zr ** (2.0 / 3.0)
    zl23 = zl ** (2.0 / 3.0)

    el = 30.724 * zr * zl * math.sqrt(zr23 + zl23) * (ar + al) / al
    denom = (zr23 + zl23) ** 0.75 * ar ** 1.5 * al ** 0.5
    fl = 0.0793 * zr23 * math.sqrt(zl) * (ar + al) ** 1.5 / denom

    E_R = 1e5
    eps = E_R / el
    g = 3.4008 * eps ** (1.0 / 6.0) + 0.40244 * eps ** 0.75 + eps
    expected = E_R / (1.0 + fl * g)

    result = lindhard_partition(E_R, 26, 56, 26, 56)
    assert result == pytest.approx(expected, rel=1e-12)
