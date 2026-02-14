"""Pydantic schemas for ChemConvert MVP."""
from __future__ import annotations

from enum import Enum
from typing import Any, Literal

from pydantic import BaseModel, Field


class InputType(str, Enum):
    SMILES = "smiles"
    INCHI = "inchi"
    MOL = "mol"
    SDF = "sdf"
    FISCHER_JSON = "fischer_json"


class OutputType(str, Enum):
    SVG = "svg"
    SMILES = "smiles"
    INCHI = "inchi"
    MOL = "mol"
    FISCHER = "fischer"


class CIPLabel(str, Enum):
    R = "R"
    S = "S"
    E = "E"
    Z = "Z"
    UNKNOWN = "unknown"


class AtomModel(BaseModel):
    idx: int
    atomic_num: int
    symbol: str
    formal_charge: int = 0
    isotope: int | None = None
    explicit_hs: int = 0
    aromatic: bool = False


class BondModel(BaseModel):
    idx: int
    begin_atom_idx: int
    end_atom_idx: int
    order: float
    aromatic: bool = False


class ChiralCenterSummary(BaseModel):
    atom_idx: int
    cip: CIPLabel | Literal["unknown"]


class DoubleBondSummary(BaseModel):
    bond_idx: int
    configuration: CIPLabel | Literal["unknown"]


class StereoSummary(BaseModel):
    chiral_centers: list[ChiralCenterSummary] = Field(default_factory=list)
    double_bonds: list[DoubleBondSummary] = Field(default_factory=list)


class ChemCore(BaseModel):
    atoms: list[AtomModel]
    bonds: list[BondModel]
    stereo: StereoSummary = Field(default_factory=StereoSummary)
    warnings: list[str] = Field(default_factory=list)
    source_molblock: str | None = None


class FischerTargetSpec(BaseModel):
    anchor_atom_idx: int | None = None
    main_chain: list[int] = Field(default_factory=list)
    orientation: Literal["vertical", "horizontal"] = "vertical"


class NewmanTargetSpec(BaseModel):
    front_atom_idx: int | None = None
    back_atom_idx: int | None = None


class HaworthTargetSpec(BaseModel):
    ring_atoms: list[int] = Field(default_factory=list)


class ChairTargetSpec(BaseModel):
    ring_atoms: list[int] = Field(default_factory=list)


class TargetSpec(BaseModel):
    fischer: FischerTargetSpec | None = None
    newman: NewmanTargetSpec | None = None
    haworth: HaworthTargetSpec | None = None
    chair: ChairTargetSpec | None = None


class ConvertRequest(BaseModel):
    input_type: InputType
    input_payload: str | dict[str, Any]
    output_type: OutputType
    output_options: dict[str, Any] = Field(default_factory=dict)
    target_spec: TargetSpec | None = None


class ConvertMeta(BaseModel):
    canonical_smiles: str
    stereo_summary: StereoSummary


class ConvertResponse(BaseModel):
    output_payload: str | dict[str, Any]
    warnings: list[str]
    meta: ConvertMeta
