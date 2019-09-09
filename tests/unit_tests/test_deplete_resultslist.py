"""Tests the ResultsList class"""

from pathlib import Path

import numpy as np
import pytest
import openmc.deplete


@pytest.fixture
def res():
    """Load the reference results"""
    filename = (Path(__file__).parents[1] / 'regression_tests' / 'deplete'
                / 'test_reference.h5')
    return openmc.deplete.ResultsList.from_hdf5(filename)


def test_get_atoms(res):
    """Tests evaluating single nuclide concentration."""
    t, n = res.get_atoms("1", "Xe135")

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    n_ref = [6.67473282e+08, 3.46601139e+14, 3.76391716e+14, 3.54026315e+14]

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(n, n_ref)


def test_get_reaction_rate(res):
    """Tests evaluating reaction rate."""
    t, r = res.get_reaction_rate("1", "Xe135", "(n,gamma)")

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    n_ref = [6.67473282e+08, 3.46601139e+14, 3.76391716e+14, 3.54026315e+14]
    xs_ref = [3.60702482e-05, 5.05757116e-05, 3.48167553e-05, 2.15382652e-05]

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(r, np.array(n_ref) * xs_ref)


def test_get_eigenvalue(res):
    """Tests evaluating eigenvalue."""
    t, k = res.get_eigenvalue()

    t_ref = [0.0, 1296000.0, 2592000.0, 3888000.0]
    k_ref = [1.19344045, 1.16489694, 1.18740893, 1.1635849]
    u_ref = [0.0348618196, 0.0247216694, 0.0364101555, 0.0293198358]

    np.testing.assert_allclose(t, t_ref)
    np.testing.assert_allclose(k[:, 0], k_ref)
    np.testing.assert_allclose(k[:, 1], u_ref)
