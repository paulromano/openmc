#!/usr/bin/env python3
"""Generate OpenMC's compact PNNL material library from the PNNL JSON data."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


EXPECTED_MATERIALS = 411
FRACTION_TOLERANCE = 5.0e-6
REPORT_URL = (
    "https://www.pnnl.gov/main/publications/external/technical_reports/"
    "PNNL-15870Rev2.pdf"
)
REPORT_SHA256 = (
    "72b26dba2c3b5583b86fe5d5fe27a43d2890331d0515ce787f1c18fd7321cee6"
)
SOURCE_URL = "https://compendium.cwmd.pnnl.gov/"
SOURCE_MIRROR_URL = (
    "https://raw.githubusercontent.com/pyne/materials-compendium/"
    "ae2f6c43f422737d659a8b19081586fb0bad09aa/"
    "src/materials_compendium/MaterialsCompendium.json"
)

# The PNNL downloadable JSON is the source of record for the bundled library;
# the pinned PyNE copy is used as a stable mirror of that data. The Rev. 2 PDF
# is a companion publication and has known differences from the downloadable
# data, most noticeably the density and elemental composition of LYSO.
DATA_ORIGIN = "PNNL Materials Compendium downloadable JSON"
REPORT_NOTE = (
    "Values are taken from the PNNL downloadable JSON, which is authoritative "
    "for this library. The Rev. 2 PDF has known differences, notably for LYSO."
)

# Names in the downloadable JSON that differ from the headings in the report.
NAME_OVERRIDES = {
    122: "Explosive Compound, AN",
    123: "Explosive Compound, EGDN",
    124: "Explosive Compound, HMX",
    125: "Explosive Compound, NC",
    126: "Explosive Compound, NG",
    127: "Explosive Compound, PETN",
    128: "Explosive Compound, RDX",
    129: "Explosive Compound, TNT",
    163: "He-3 Proportional Gas",
    164: "He-4 Gas Detector",
    203: "Lutetium Iodide",
    206: "Lutetium Yttrium OxyorthoSilicate: 0.5 atom% Cerium (LYSO)",
    321: "Sodium Iodide with 8 wt% Lithium - 0.10 wt% Thalium doped",
}


def _positive_float(value, description):
    """Convert a value to a finite, positive float."""
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{description} must be finite and positive")
    return value


def _nonnegative_float(value, description):
    """Convert a value to a finite, nonnegative float."""
    value = float(value)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"{description} must be finite and nonnegative")
    return value


def convert(source: Path) -> dict:
    """Convert the PNNL material data to OpenMC's compact representation."""
    source_bytes = source.read_bytes()
    source_data = json.loads(source_bytes)
    pnnl_materials = source_data["data"]

    if len(pnnl_materials) != EXPECTED_MATERIALS:
        raise ValueError(
            f"Expected {EXPECTED_MATERIALS} materials; found "
            f"{len(pnnl_materials)}"
        )

    materials = {}
    for material_data in pnnl_materials:
        material_number = material_data["MatNum"]
        name = NAME_OVERRIDES.get(
            material_number, " ".join(material_data["Name"].split())
        )
        if name in materials:
            raise ValueError(f"Duplicate material name: {name}")

        density = _positive_float(
            material_data["Density"], f"Density for {name}"
        )
        elements = {}
        nuclides = {}

        for element_data in material_data["Elements"]:
            element = element_data["Element"]
            element_fraction = _nonnegative_float(
                element_data["AtomFraction"],
                f"Atom fraction for {element} in {name}",
            )

            if element_data.get("NonIsotopic") is True:
                if element_fraction == 0.0:
                    continue
                if element in elements:
                    raise ValueError(f"Duplicate element {element} in {name}")
                elements[element] = element_fraction
            else:
                isotope_total = 0.0
                for isotope_data in element_data["Isotopes"]:
                    nuclide = isotope_data["Isotope"].replace("-", "")
                    fraction = _nonnegative_float(
                        isotope_data["AtomFraction"],
                        f"Atom fraction for {nuclide} in {name}",
                    )
                    isotope_total += fraction
                    if fraction == 0.0:
                        continue
                    if nuclide in nuclides:
                        raise ValueError(f"Duplicate nuclide {nuclide} in {name}")
                    nuclides[nuclide] = fraction

                if not math.isclose(
                    isotope_total,
                    element_fraction,
                    rel_tol=0.0,
                    abs_tol=FRACTION_TOLERANCE,
                ):
                    raise ValueError(
                        f"Isotope fractions for {element} in {name} sum to "
                        f"{isotope_total}, not {element_fraction}"
                    )

        total = sum(elements.values()) + sum(nuclides.values())
        if not math.isclose(
            total, 1.0, rel_tol=0.0, abs_tol=FRACTION_TOLERANCE
        ):
            raise ValueError(f"Atom fractions for {name} sum to {total}")

        material = {"density": density}
        if elements:
            material["elements"] = elements
        if nuclides:
            material["nuclides"] = nuclides
        materials[name] = material

    return {
        "schema_version": 1,
        "source": {
            "title": (
                "Compendium of Material Composition Data for Radiation "
                "Transport Modeling"
            ),
            "report": "PNNL-15870 Rev. 2",
            "report_url": REPORT_URL,
            "report_sha256": REPORT_SHA256,
            "data_origin": DATA_ORIGIN,
            "data_url": SOURCE_URL,
            "data_mirror_url": SOURCE_MIRROR_URL,
            "site_version": source_data["siteVersion"],
            "data_sha256": hashlib.sha256(source_bytes).hexdigest(),
            "note": REPORT_NOTE,
        },
        "density_units": "g/cm3",
        "percent_type": "ao",
        "materials": dict(sorted(materials.items())),
    }


def main():
    """Parse command-line arguments and write the compact material library."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="PNNL MaterialsCompendium.json")
    parser.add_argument("output", type=Path, help="Output JSON file")
    args = parser.parse_args()

    library = convert(args.source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(library, indent=2) + "\n")


if __name__ == "__main__":
    main()
