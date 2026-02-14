"""ChemCore conversion helpers: RDKit Mol <-> ChemCore."""
from __future__ import annotations

from rdkit import Chem

from app.models import (
    AtomModel,
    BondModel,
    CIPLabel,
    ChiralCenterSummary,
    ChemCore,
    DoubleBondSummary,
    StereoSummary,
)


def _extract_stereo_summary(mol: Chem.Mol, warnings: list[str] | None = None) -> StereoSummary:
    summary = StereoSummary()
    warnings = warnings if warnings is not None else []

    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    for atom_idx, label in Chem.FindMolChiralCenters(
        mol, includeUnassigned=True, useLegacyImplementation=False
    ):
        cip = CIPLabel.UNKNOWN
        if label == "R":
            cip = CIPLabel.R
        elif label == "S":
            cip = CIPLabel.S
        else:
            warnings.append(f"Chiral center at atom {atom_idx} is unassigned.")
        summary.chiral_centers.append(ChiralCenterSummary(atom_idx=atom_idx, cip=cip))

    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        if bond.GetIsAromatic():
            continue

        stereo = bond.GetStereo()
        conf = CIPLabel.UNKNOWN
        if stereo == Chem.rdchem.BondStereo.STEREOE:
            conf = CIPLabel.E
        elif stereo == Chem.rdchem.BondStereo.STEREOZ:
            conf = CIPLabel.Z
        else:
            warnings.append(f"Double bond {bond.GetIdx()} has unknown E/Z stereochemistry.")
        summary.double_bonds.append(DoubleBondSummary(bond_idx=bond.GetIdx(), configuration=conf))

    return summary


def mol_to_chemcore(mol: Chem.Mol) -> ChemCore:
    mol = Chem.Mol(mol)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    warnings: list[str] = []
    stereo = _extract_stereo_summary(mol, warnings)

    atoms = [
        AtomModel(
            idx=atom.GetIdx(),
            atomic_num=atom.GetAtomicNum(),
            symbol=atom.GetSymbol(),
            formal_charge=atom.GetFormalCharge(),
            isotope=atom.GetIsotope() or None,
            explicit_hs=atom.GetNumExplicitHs(),
            aromatic=atom.GetIsAromatic(),
        )
        for atom in mol.GetAtoms()
    ]

    bonds = [
        BondModel(
            idx=bond.GetIdx(),
            begin_atom_idx=bond.GetBeginAtomIdx(),
            end_atom_idx=bond.GetEndAtomIdx(),
            order=float(bond.GetBondTypeAsDouble()),
            aromatic=bond.GetIsAromatic(),
        )
        for bond in mol.GetBonds()
    ]

    return ChemCore(
        atoms=atoms,
        bonds=bonds,
        stereo=stereo,
        warnings=warnings,
        source_molblock=Chem.MolToMolBlock(mol),
    )


def chemcore_to_mol(core: ChemCore) -> Chem.Mol:
    if core.source_molblock:
        mol = Chem.MolFromMolBlock(core.source_molblock, sanitize=True, removeHs=False)
        if mol is not None:
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            return mol

    rw_mol = Chem.RWMol()
    atom_map: dict[int, int] = {}

    for atom in core.atoms:
        rd_atom = Chem.Atom(atom.atomic_num)
        rd_atom.SetFormalCharge(atom.formal_charge)
        rd_atom.SetNumExplicitHs(atom.explicit_hs)
        rd_atom.SetNoImplicit(True)
        if atom.isotope:
            rd_atom.SetIsotope(atom.isotope)
        rd_atom.SetIsAromatic(atom.aromatic)
        atom_map[atom.idx] = rw_mol.AddAtom(rd_atom)

    def _bond_type(order: float) -> Chem.BondType:
        if order == 1.0:
            return Chem.BondType.SINGLE
        if order == 2.0:
            return Chem.BondType.DOUBLE
        if order == 3.0:
            return Chem.BondType.TRIPLE
        if order == 1.5:
            return Chem.BondType.AROMATIC
        return Chem.BondType.UNSPECIFIED

    for bond in core.bonds:
        rw_mol.AddBond(atom_map[bond.begin_atom_idx], atom_map[bond.end_atom_idx], _bond_type(bond.order))
        rd_bond = rw_mol.GetBondBetweenAtoms(atom_map[bond.begin_atom_idx], atom_map[bond.end_atom_idx])
        if rd_bond is not None:
            rd_bond.SetIsAromatic(bond.aromatic)

    mol = rw_mol.GetMol()
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol
