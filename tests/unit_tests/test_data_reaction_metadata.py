"""Preservation of ENDF fields unavailable in ACE and older HDF5 files."""

from types import SimpleNamespace

import h5py
import numpy as np
import pytest

import openmc.data
from openmc.data.endf import Evaluation


@pytest.fixture
def evaluation():
    # A minimal MF=3 section with distinct QM and QI and a breakup flag.
    def record(*fields):
        return ''.join(f'{x:11.4E}' if isinstance(x, float) else f'{x:11d}'
                       for x in fields) + '\n'

    ev = Evaluation.__new__(Evaluation)
    ev.target = dict(atomic_number=95, mass_number=242, isomeric_state=2,
                     mass=240.0, temperature=0.0, excitation_energy=48600.0)
    ev.section = {(3, 51): (
        record(95242.0, 240.0, 0, 0, 0, 0)
        + record(0.0, 48600.0, 0, 33, 1, 2)
        + record(2, 2, 0, 0, 0, 0)
        + record(1.0, 0.0, 2.0e7, 1.0, 0.0, 0.0))}
    # An empty resonance section avoids invoking unrelated resonance models.
    ev.section[2, 151] = record(95242.0, 240.0, 0, 0, 0, 0)
    ev.reaction_list = [(3, 51, 4, 0)]
    return ev


def test_endf_metadata(evaluation):
    data = openmc.data.IncidentNeutron.from_endf(evaluation)
    assert data.excitation_energy == 48600.0
    assert data[51].q_mass_difference == 0.0
    assert data[51].q_reaction == 48600.0
    assert data[51].breakup_flag == 33


@pytest.mark.parametrize('values', [(None, None, None), (0.0, 0.0, 0),
                                   (48600.0, -1.2e6, 33)])
def test_metadata_hdf5(tmp_path, values):
    data = openmc.data.IncidentNeutron('Am242_m1', 95, 242, 1, 240.0, [0.0])
    data.energy['0K'] = np.array([1.0, 2.0e7])
    rx = openmc.data.Reaction(51)
    rx.q_reaction = 48600.0
    rx.xs['0K'] = openmc.data.Tabulated1D(data.energy['0K'], [0.0, 1.0])
    data.reactions[51] = rx
    data.excitation_energy, rx.q_mass_difference, rx.breakup_flag = values
    path = tmp_path / 'data.h5'
    data.export_to_hdf5(path)
    with h5py.File(path, 'r+') as handle:
        attrs = handle[f'{data.name}/reactions/reaction_051'].attrs
        assert attrs['q_reaction'] == rx.q_reaction
        assert 'Q_value' not in attrs
        if values[0] is None:
            # A legacy file has neither the new attributes nor a minor bump.
            handle.attrs['version'] = [3, 0]
            assert 'excitation_energy' not in handle[data.name].attrs
            attrs['Q_value'] = attrs['q_reaction']
            del attrs['q_reaction']
            assert 'q_mass_difference' not in attrs
            assert 'breakup_flag' not in attrs
    restored = openmc.data.IncidentNeutron.from_hdf5(path)
    assert (restored.excitation_energy, restored[51].q_mass_difference,
            restored[51].breakup_flag) == values
    assert restored[51].q_reaction == rx.q_reaction


def test_njoy_metadata(monkeypatch, evaluation):
    import openmc.data.neutron as neutron

    # Isolate the ENDF enrichment step from the external NJOY executable.
    data = openmc.data.IncidentNeutron('Am242_m1', 95, 242, 1, 240.0, [0.0])
    data.energy['0K'] = np.array([1.0, 2.0e7])
    data.reactions = {mt: openmc.data.Reaction(mt) for mt in (51, 444)}
    data[51].q_reaction = 48000.0  # Preserve the processed ACE reaction Q value.
    monkeypatch.setattr(neutron, 'make_ace', lambda *args, **kwargs: None)
    monkeypatch.setattr(neutron, 'Library',
                        lambda path: SimpleNamespace(tables=[object()]))
    monkeypatch.setattr(neutron.IncidentNeutron, 'from_ace', lambda table: data)
    result = neutron.IncidentNeutron.from_njoy(
        'unused.endf', evaluation=evaluation, heatr=False)
    assert result.name == 'Am242_m2'
    assert result.metastable == 2
    assert result.excitation_energy == 48600.0
    assert result[51].q_mass_difference == 0.0
    assert result[51].breakup_flag == 33
    assert result[51].q_reaction == 48000.0
    assert result[444].q_mass_difference is None
    assert result[444].breakup_flag is None


def test_q_value_alias():
    reaction = openmc.data.Reaction(51)
    reaction.q_value = 48600.0
    assert reaction.q_reaction == 48600.0
    reaction.q_reaction = -1.0e6
    assert reaction.q_value == -1.0e6
    with pytest.raises(TypeError):
        reaction.q_value = 'invalid'


def test_q_reaction_precedence(tmp_path):
    with h5py.File(tmp_path / 'reaction.h5', 'w') as handle:
        reaction = openmc.data.Reaction(51)
        reaction.q_reaction = 48600.0
        reaction.to_hdf5(handle)
        handle.attrs['Q_value'] = -1.0e6
        restored = openmc.data.Reaction.from_hdf5(handle, {})
        assert restored.q_reaction == 48600.0
