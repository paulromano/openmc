#!/usr/bin/env python

import os
import sys

sys.path.insert(0, os.pardir)
from testing_harness import PyAPITestHarness
import openmc


class SensitivityIFPTest(PyAPITestHarness):
    def __init__(self, statepoint_name):
        super().__init__(statepoint_name)

        self._model.settings.inactive = 30
        self._model.settings.batches = 100

        ll, ur = self._model.geometry.bounding_box

        mesh = openmc.Mesh(1)
        mesh.lower_left = ll
        mesh.upper_right = ur
        mesh.dimension = (1, 1, 1)

        s = openmc.Sensitivity(1)
        s.mesh = mesh
        s.energies = [0.0, 20.0e6]
        s.nuclides = ['U235', 'Zr90', 'H1']
        s.scores = ['total', 'scatter', '(n,2n)', 'fission', 'nubar',
                    'chi', 'capture']

        self._model.sensitivities = openmc.Sensitivities([s])
        self._model.sensitivities.adjoint_method = 'ifp'
        self._model.sensitivities.ifp_block_length = 10

    def _build_inputs(self):
        super()._build_inputs()
        self._model.sensitivities.export_to_xml()


if __name__ == '__main__':
    harness = SensitivityIFPTest('statepoint.100.h5')
    harness.main()
