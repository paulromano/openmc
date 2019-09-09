import random
from math import sqrt

import numpy as np
import openmc
import openmc.model
import pytest

from tests.testing_harness import PyAPITestHarness


@pytest.fixture
def triso_model():
    model = openmc.model.Model()

    # Define TRISO matrials
    fuel = openmc.Material()
    fuel.set_density('g/cm3', 10.5)
    fuel.add_nuclide('U235', 1.0)

    porous_carbon = openmc.Material()
    porous_carbon.set_density('g/cm3', 0.1)
    porous_carbon.add_nuclide('C0', 1.0)
    porous_carbon.add_s_alpha_beta('c_Graphite')

    sic = openmc.Material()
    sic.set_density('g/cm3', 3.20)
    sic.add_nuclide('C0', 1.0)
    sic.add_element('Si', 1.0)

    ipyc = porous_carbon
    opyc = porous_carbon
    graphite = porous_carbon

    model.materials.extend([fuel, porous_carbon, sic])

    # Create TRISO particles
    spheres = [openmc.Sphere(r=r*1e-4)
               for r in [212.5, 312.5, 347.5, 382.5]]
    c1 = openmc.Cell(fill=fuel, region=-spheres[0])
    c2 = openmc.Cell(fill=porous_carbon, region=+spheres[0] & -spheres[1])
    c3 = openmc.Cell(fill=ipyc, region=+spheres[1] & -spheres[2])
    c4 = openmc.Cell(fill=sic, region=+spheres[2] & -spheres[3])
    c5 = openmc.Cell(fill=opyc, region=+spheres[3])
    inner_univ = openmc.Universe(cells=[c1, c2, c3, c4, c5])

    # Define box to contain TRISO particles in
    w = 0.25
    min_x = openmc.XPlane(x0=-w/2, boundary_type='reflective')
    max_x = openmc.XPlane(x0=w/2, boundary_type='reflective')
    min_y = openmc.YPlane(y0=-w/2, boundary_type='reflective')
    max_y = openmc.YPlane(y0=w/2, boundary_type='reflective')
    min_z = openmc.ZPlane(z0=-w/2, boundary_type='reflective')
    max_z = openmc.ZPlane(z0=w/2, boundary_type='reflective')
    box_region = +min_x & -max_x & +min_y & -max_y & +min_z & -max_z

    outer_radius = 422.5*1e-4
    centers = openmc.model.pack_spheres(radius=outer_radius, region=box_region,
                                        num_spheres=10)
    trisos = [openmc.model.TRISO(outer_radius, inner_univ, c)
              for c in centers]

    for t in trisos:
        box_region &= ~t.region
    box = openmc.Cell(fill=graphite, region=box_region)

    model.geometry = openmc.Geometry(trisos + [box])

    model.settings.batches = 4
    model.settings.inactive = 0
    model.settings.particles = 100
    model.settings.source = openmc.Source(
        space=openmc.stats.Point(),
        energy=openmc.stats.Discrete([0.0253], [1.0])
    )

    return model


def test_triso(triso_model):
    harness = PyAPITestHarness('statepoint.4.h5', triso_model)
    harness.main()
