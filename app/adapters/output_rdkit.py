"""RDKit output adapters: ChemCore -> target formats."""
from __future__ import annotations

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from app.core.chemcore import chemcore_to_mol
from app.models import ChemCore


def to_svg_bondline(core: ChemCore, options: dict | None = None) -> str:
    options = options or {}
    mol = chemcore_to_mol(core)

    if mol.GetNumConformers() == 0:
        AllChem.Compute2DCoords(mol)

    rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=True, addChiralHs=True, wedgeBonds=True)

    width = int(options.get("width", 450))
    height = int(options.get("height", 320))
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    draw_options = drawer.drawOptions()
    draw_options.addStereoAnnotation = bool(options.get("stereo_annotation", True))

    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def to_smiles(core: ChemCore) -> str:
    return Chem.MolToSmiles(chemcore_to_mol(core), isomericSmiles=True, canonical=True)


def to_inchi(core: ChemCore) -> str:
    return Chem.MolToInchi(chemcore_to_mol(core))


def to_mol(core: ChemCore) -> str:
    return Chem.MolToMolBlock(chemcore_to_mol(core))


def to_fischer(core: ChemCore, target_spec: dict | None = None) -> dict:
    """TODO(plugin): ChemCore -> Fischer JSON."""
    raise NotImplementedError("Fischer output plugin is not implemented yet")
