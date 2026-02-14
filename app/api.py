"""FastAPI app and convert endpoint."""
from __future__ import annotations

from fastapi import FastAPI, HTTPException

from app.adapters.input_rdkit import (
    parse_fischer_json,
    parse_inchi,
    parse_molfile,
    parse_smiles,
)
from app.adapters.output_rdkit import to_fischer, to_inchi, to_mol, to_smiles, to_svg_bondline
from app.models import ConvertMeta, ConvertRequest, ConvertResponse, InputType, OutputType
from app.services.normalize import normalize

app = FastAPI(title="ChemConvert", version="0.2.0")


@app.post("/convert", response_model=ConvertResponse)
def convert(req: ConvertRequest) -> ConvertResponse:
    try:
        if req.input_type == InputType.SMILES:
            core = parse_smiles(str(req.input_payload))
        elif req.input_type == InputType.INCHI:
            core = parse_inchi(str(req.input_payload))
        elif req.input_type in {InputType.MOL, InputType.SDF}:
            core = parse_molfile(str(req.input_payload))
        elif req.input_type == InputType.FISCHER_JSON:
            if not isinstance(req.input_payload, dict):
                raise ValueError("fischer_json input payload must be object")
            core = parse_fischer_json(req.input_payload)
        else:
            raise ValueError(f"Unsupported input type: {req.input_type}")

        core = normalize(core)

        if req.output_type == OutputType.SVG:
            payload = to_svg_bondline(core, req.output_options)
        elif req.output_type == OutputType.SMILES:
            payload = to_smiles(core)
        elif req.output_type == OutputType.INCHI:
            payload = to_inchi(core)
        elif req.output_type == OutputType.MOL:
            payload = to_mol(core)
        elif req.output_type == OutputType.FISCHER:
            payload = to_fischer(core, req.target_spec.model_dump() if req.target_spec else None)
        else:
            raise ValueError(f"Unsupported output type: {req.output_type}")

        return ConvertResponse(
            output_payload=payload,
            warnings=core.warnings,
            meta=ConvertMeta(canonical_smiles=to_smiles(core), stereo_summary=core.stereo),
        )
    except NotImplementedError as exc:
        raise HTTPException(status_code=501, detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
