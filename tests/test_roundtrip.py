from __future__ import annotations

import pytest

pytest.importorskip("fastapi")
pytest.importorskip("rdkit")

from fastapi.testclient import TestClient

from app.api import app
from app.adapters.input_rdkit import parse_smiles
from app.adapters.output_rdkit import to_smiles, to_svg_bondline


client = TestClient(app)


def test_smiles_roundtrip_canonical() -> None:
    core = parse_smiles("C[C@H](O)C(=O)O")
    canonical_1 = to_smiles(core)
    canonical_2 = to_smiles(parse_smiles(canonical_1))
    assert canonical_1 == canonical_2


def test_svg_contains_bond_paths_for_stereo_molecule() -> None:
    core = parse_smiles("F[C@H](Cl)Br")
    svg = to_svg_bondline(core, {"stereo_annotation": True})
    assert "<svg" in svg
    assert "class='bond-" in svg or 'class="bond-' in svg


def test_convert_meta_stereo_summary_schema() -> None:
    resp = client.post(
        "/convert",
        json={
            "input_type": "smiles",
            "input_payload": "F/C=C/F",
            "output_type": "smiles",
            "output_options": {},
        },
    )
    assert resp.status_code == 200
    data = resp.json()
    assert "stereo_summary" in data["meta"]
    assert isinstance(data["meta"]["stereo_summary"]["double_bonds"], list)
