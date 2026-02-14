"""RDKit input adapters: external format -> ChemCore."""
from __future__ import annotations

from rdkit import Chem

from app.core.chemcore import mol_to_chemcore
from app.models import ChemCore


def parse_smiles(smiles: str) -> ChemCore:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES payload")
    return mol_to_chemcore(mol)


def parse_inchi(inchi: str) -> ChemCore:
    mol = Chem.MolFromInchi(inchi)
    if mol is None:
        raise ValueError("Invalid InChI payload")
    return mol_to_chemcore(mol)


def parse_molfile(mol_block: str) -> ChemCore:
    mol = Chem.MolFromMolBlock(mol_block, sanitize=True, removeHs=False)
    if mol is None:
        raise ValueError("Invalid MOL/SDF payload")
    return mol_to_chemcore(mol)


def parse_fischer_json(fischer_json: dict) -> ChemCore:
    """TODO(plugin): Fischer JSON -> ChemCore."""
    raise NotImplementedError("FischerJSON parser plugin is not implemented yet")
