from collections import Iterable
from numbers import Integral, Real
from xml.etree import ElementTree as ET

from six import string_types
import numpy as np
import pandas as pd

import openmc
import openmc.checkvalue as cv
from openmc.clean_xml import clean_xml_indentation


class Sensitivity(object):
    def __init__(self, sensitivity_id):
        self.id = sensitivity_id
        self._mesh = None
        self._importance_mesh = None
        self._energies = None
        self._nuclides = None
        self._scores = None

    @property
    def id(self):
        return self._id

    @property
    def mesh(self):
        return self._mesh

    @property
    def importance_mesh(self):
        return self._importance_mesh

    @property
    def energies(self):
        return self._energies

    @property
    def nuclides(self):
        return self._nuclides

    @property
    def scores(self):
        return self._scores

    @id.setter
    def id(self, sid):
        cv.check_type('sensitivity ID', sid, Integral)
        self._id = sid

    @mesh.setter
    def mesh(self, mesh):
        cv.check_type('sensitivity mesh', mesh, openmc.Mesh)
        self._mesh = mesh

    @importance_mesh.setter
    def importance_mesh(self, importance_mesh):
        cv.check_type('sensitivity importance mesh', importance_mesh,
                      openmc.Mesh)
        self._importance_mesh = importance_mesh

    @energies.setter
    def energies(self, energies):
        cv.check_type('sensitivity energies', energies, Iterable, Real)
        self._energies = np.asarray(energies)

    @nuclides.setter
    def nuclides(self, nuclides):
        cv.check_type('sensitivity nuclides', nuclides, Iterable, string_types)
        self._nuclides = nuclides

    @scores.setter
    def scores(self, scores):
        cv.check_type('sensitivity scores', scores, Iterable, string_types)
        self._scores = scores

    @classmethod
    def from_hdf5(cls, group):
        sid = int(group.name.split('/')[-1].lstrip('sensitivity '))
        data = cls(sid)
        data.energies = group['energies'].value
        data.nuclides = [i.decode() for i in group['nuclides'].value]
        data.scores = [i.decode() for i in group['scores'].value]
        data.n_realizations = group['n_realizations'].value
        mesh_id = group['mesh_id'].value
        data.mesh = openmc.Mesh.from_hdf5(group.parent['mesh {}'.format(mesh_id)])
        data._results = group['results'].value
        return data

    def get_pandas_dataframe(self):
        columns = ['Energy low [eV]', 'Energy high [eV]', 'Mesh', 'Score', 'Nuclide',
                   'mean', 'std. dev.']

        n = self.n_realizations
        mean = self._results[..., 0] / n
        stdev = np.sqrt((self._results[..., 1]/n - mean*mean)/(n - 1))

        n_mesh = self._results.shape[1]
        records = []
        for i, (e_lo, e_hi) in enumerate(zip(self.energies[:-1], self.energies[1:])):
            for j in range(n_mesh):
                for k, score in enumerate(self.scores):
                    for m, nuclide in enumerate(self.nuclides):
                        records.append((e_lo, e_hi, j + 1, score, nuclide,
                                        mean[i, j, k, m], stdev[i, j, k, m]))

        return pd.DataFrame.from_records(records, columns=columns)

    def to_xml_element(self):
        element = ET.Element("sensitivity")
        element.set("id", str(self._id))

        subelement = ET.SubElement(element, "mesh")
        subelement.text = str(self.mesh.id)
        subelement = ET.SubElement(element, "impmesh")
        subelement.text = str(self.importance_mesh.id)
        subelement = ET.SubElement(element, "energies")
        subelement.text = ' '.join(str(E) for E in self.energies)
        subelement = ET.SubElement(element, "nuclides")
        subelement.text = ' '.join(self.nuclides)
        subelement = ET.SubElement(element, "scores")
        subelement.text = ' '.join(self.scores)

        return element


class Sensitivities(cv.CheckedList):
    valid_methods = ['ifp', 'clutch-ifp', 'clutch-fm', 'gpt-ifp',
                     'gpt-clutch-ifp', 'gpt-clutch-fm']

    def __init__(self, items=None):
        super(Sensitivities, self).__init__(Sensitivity,
                                            'sensitivity collection')
        if items is not None:
            self += items

    @property
    def adjoint_method(self):
        return self._adjoint_method

    @property
    def ifp_block_length(self):
        return self._ifp_block_length

    @adjoint_method.setter
    def adjoint_method(self, method):
        cv.check_value('adjoint method', method, Sensitivities.valid_methods)
        self._adjoint_method = method

    @ifp_block_length.setter
    def ifp_block_length(self, length):
        cv.check_type('IFP block length', length, Integral)
        cv.check_greater_than('IFP block length', length, 1)
        self._ifp_block_length = length

    def export_to_xml(self, path='sensitivities.xml'):
        root_element = ET.Element("sensitivities")

        # Write adjoint method and IFP block length
        subelement = ET.SubElement(root_element, 'adjoint_method')
        subelement.text = self.adjoint_method
        if 'fm' not in self.adjoint_method:
            subelement = ET.SubElement(root_element, 'ifp_block_length')
            subelement.text = str(self.ifp_block_length)

        # Write mesh/sensitivity elements
        already_written = set()
        for s in self:
            if s.mesh is not None:
                if s.mesh not in already_written:
                    root_element.append(s.mesh.to_xml_element())
                    already_written.add(s.mesh)
            if s.importance_mesh is not None:
                if s.importance_mesh not in already_written:
                    root_element.append(s.importance_mesh.to_xml_element())
                    already_written.add(s.importance_mesh)
            root_element.append(s.to_xml_element())

        # Clean the indentation in the file to be user-readable
        clean_xml_indentation(root_element)

        # Write the XML Tree to the sensitivities.xml file
        tree = ET.ElementTree(root_element)
        tree.write(path, xml_declaration=True, encoding='utf-8', method="xml")
