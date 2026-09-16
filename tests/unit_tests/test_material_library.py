from importlib import resources
import json

import lxml.etree as ET
import pytest

import openmc


@pytest.fixture
def cross_sections(tmp_path):
    """Cross section index covering the natural elements used in tests."""
    root = ET.Element('cross_sections')
    for nuclide in (
        'O16', 'O17', 'Na23', 'Si28', 'Si29', 'Si30', 'Cl35', 'Cl37',
        'Y89', 'Ce136', 'Ce138', 'Ce140', 'Ce142', 'Cs133', 'Lu175',
        'Lu176'
    ):
        ET.SubElement(
            root,
            'library',
            materials=nuclide,
            path=f'{nuclide}.h5',
            type='neutron',
        )
    path = tmp_path / 'cross_sections.xml'
    ET.ElementTree(root).write(path)
    return path


def test_natural_elements(cross_sections):
    """Natural elements use nuclides present in the cross section library."""
    with openmc.config.patch('cross_sections', cross_sections):
        material = openmc.Material.from_library('Sodium Oxide')

    assert material.name == 'Sodium Oxide'
    assert material.density == pytest.approx(2.27)
    assert material.density_units == 'g/cm3'

    nuclides = {nuc.name: nuc.percent for nuc in material.nuclides}
    assert set(nuclides) == {'O16', 'O17', 'Na23'}
    assert nuclides['O16'] + nuclides['O17'] == pytest.approx(0.333333)
    assert nuclides['Na23'] == pytest.approx(0.666667)


def test_isotopic_components(cross_sections):
    """Non-natural isotopic compositions are retained explicitly."""
    he3 = openmc.Material.from_library('He-3 Proportional Gas')
    assert he3.nuclides == [('He3', 1.0, 'ao')]

    with openmc.config.patch('cross_sections', cross_sections):
        clyc = openmc.Material.from_library(
            'Cesium Lithium Yttrium Chloride (CLYC) with 95% Li6 Enrichment'
        )
    lithium = {
        nuc.name: nuc.percent
        for nuc in clyc.nuclides
        if nuc.name.startswith('Li')
    }
    assert lithium == pytest.approx({'Li6': 0.095682, 'Li7': 0.004318})

    leu = openmc.Material.from_library('Uranium, Low Enriched (LEU)')
    uranium = {nuc.name: nuc.percent for nuc in leu.nuclides}
    assert uranium == pytest.approx({
        'U234': 0.000271,
        'U235': 0.030372,
        'U236': 0.000139,
        'U238': 0.969217,
    })

    aged_pu = openmc.Material.from_library(
        'Plutonium, Aged WGPu (A: 4-7% Pu240)'
    )
    assert 'Am241' in aged_pu.get_nuclides()
    assert sum(nuc.percent for nuc in aged_pu.nuclides) == pytest.approx(1.0)


def test_lyso_uses_pdf_composition(cross_sections):
    """LYSO includes the cerium dopant specified in the Rev. 2 PDF."""
    with openmc.config.patch('cross_sections', cross_sections):
        lyso = openmc.Material.from_library(
            'Lutetium Yttrium OxyorthoSilicate: 0.5 atom% Cerium (LYSO)'
        )

    fractions = {nuc.name: nuc.percent for nuc in lyso.nuclides}
    assert lyso.density == pytest.approx(7.25)
    assert fractions['O16'] + fractions['O17'] == pytest.approx(0.621875)
    assert sum(
        fraction for name, fraction in fractions.items()
        if name.startswith('Si')
    ) == pytest.approx(0.124375)
    assert fractions['Y89'] == pytest.approx(0.012438)
    assert fractions['Lu175'] + fractions['Lu176'] == pytest.approx(0.236313)
    assert sum(
        fraction for name, fraction in fractions.items()
        if name.startswith('Ce')
    ) == pytest.approx(0.005)


def test_independent_materials():
    """Each lookup creates an independent material with a fresh ID."""
    first = openmc.Material.from_library('He-3 Proportional Gas')
    second = openmc.Material.from_library('He-3 Proportional Gas')

    assert first is not second
    assert first.id != second.id
    first.add_nuclide('H1', 1.0)
    assert second.get_nuclides() == ['He3']


def test_constructor_kwargs():
    """Keyword arguments are forwarded to the Material constructor."""
    material = openmc.Material.from_library(
        'He-3 Proportional Gas',
        material_id=987654,
        name='Helium-3 detector gas',
        temperature=293.15,
        volume=10.0,
        depletable=True,
    )

    assert material.id == 987654
    assert material.name == 'Helium-3 detector gas'
    assert material.temperature == 293.15
    assert material.volume == 10.0
    assert material.depletable
    assert material.density == pytest.approx(0.0001252645124733361)
    assert material.nuclides == [('He3', 1.0, 'ao')]


def test_unknown_library_and_material():
    with pytest.raises(ValueError, match="Unknown material library 'unknown'"):
        openmc.Material.from_library('Sodium Oxide', library='unknown')

    with pytest.raises(ValueError, match="Material 'unknown' not found"):
        openmc.Material.from_library('unknown')


def test_pnnl_data():
    """Check invariants of the complete bundled PNNL material library."""
    path = resources.files('openmc.data').joinpath(
        'material_libraries/pnnl_v2.json'
    )
    data = json.loads(path.read_text(encoding='utf-8'))

    assert data['schema_version'] == 1
    assert data['density_units'] == 'g/cm3'
    assert data['percent_type'] == 'ao'
    assert data['source']['data_origin'] == (
        'PNNL Materials Compendium downloadable JSON'
    )
    assert 'LYSO' in data['source']['note']
    assert data['source']['report_sha256'] == (
        '72b26dba2c3b5583b86fe5d5fe27a43d2890331d0515ce787f1c18fd7321cee6'
    )
    assert data['source']['data_sha256'] == (
        '5db6f9ca58793659e73c66cfd624736ffa59846239cd423c11eb6dbe5a56607b'
    )
    lyso = data['materials'][
        'Lutetium Yttrium OxyorthoSilicate: 0.5 atom% Cerium (LYSO)'
    ]
    assert lyso['density'] == 7.25
    assert lyso['elements'] == {
        'O': 0.621875,
        'Si': 0.124375,
        'Y': 0.012438,
        'Lu': 0.236313,
        'Ce': 0.005,
    }
    assert 'PNNL-15870 Rev. 2' in lyso['source_note']
    assert len(data['materials']) == 411
    assert sum('nuclides' in mat for mat in data['materials'].values()) == 45

    element_symbols = set(openmc.data.ATOMIC_SYMBOL.values())
    for name, material in data['materials'].items():
        assert name
        assert material['density'] > 0.0

        elements = material.get('elements', {})
        nuclides = material.get('nuclides', {})
        assert elements or nuclides
        assert set(elements) <= element_symbols
        for nuclide in nuclides:
            openmc.data.zam(nuclide)

        fractions = [*elements.values(), *nuclides.values()]
        assert all(fraction > 0.0 for fraction in fractions)
        assert sum(fractions) == pytest.approx(1.0, abs=5.0e-6)
