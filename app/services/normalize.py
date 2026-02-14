"""Normalization service for ChemCore."""
from __future__ import annotations

from rdkit import Chem

from app.core.chemcore import chemcore_to_mol, mol_to_chemcore
from app.models import ChemCore


def normalize(core: ChemCore) -> ChemCore:
    mol = chemcore_to_mol(core)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    normalized = mol_to_chemcore(mol)
    normalized.warnings = [*core.warnings, *normalized.warnings]
    return normalized
