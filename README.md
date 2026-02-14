# chemconvert

基于 **Input -> ChemCore -> Output** 的结构式自由转换 MVP，后续可扩展为插件。

## 最小目录结构

```text
chemconvert/
  app/
    api.py
    models.py
    core/chemcore.py
    services/normalize.py
    adapters/input_rdkit.py
    adapters/output_rdkit.py
    adapters/projections_fischer.py    # TODO
    tests/test_roundtrip.py
  pyproject.toml
  README.md
```

## 已支持（阶段1）

- 输入：SMILES / InChI / MOL(SDF)
- 输出：bond-line SVG / isomeric SMILES / InChI / MOL
- 统一中间层：ChemCore（分子图 + 立体摘要 + warnings）

## 插件位（阶段2）

- `parse_fischer_json(...)`（TODO）
- `to_fischer(...)`（TODO）
- `TargetSpec` 中预留 Newman/Haworth/Chair 参数结构

## 运行

```bash
uvicorn app.main:app --reload
```

## API

`POST /convert`

请求示例：

```json
{
  "input_type": "smiles",
  "input_payload": "F/C=C/F",
  "output_type": "svg",
  "output_options": {"width": 500, "height": 360, "stereo_annotation": true}
}
```

响应中的 `meta.stereo_summary`：

- `chiral_centers`: `[{atom_idx, cip: R|S|unknown}]`
- `double_bonds`: `[{bond_idx, configuration: E|Z|unknown}]`

当无法判定时，标记 `unknown` 并写入 `warnings`。
