#!/usr/bin/env python3
"""
RDKit_SDF_AddH.py

读取 SDF 文件，对每个分子进行加氢（AddHs），并输出到新的 SDF 文件。

用法:
    python RDKit_SDF_AddH.py -i input.sdf -o output_AddH.sdf
"""

import argparse
import sys

from rdkit import Chem
from rdkit.Chem import AllChem


def parse_args():
    parser = argparse.ArgumentParser(
        description="使用 RDKit 对 SDF 文件中的分子加氢 (AddHs)。"
    )
    parser.add_argument(
        "-i", "--input", required=True, help="输入 SDF 文件路径"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="输出加氢后 SDF 文件路径"
    )
    parser.add_argument(
        "--no-coords",
        action="store_true",
        help="加氢后不做坐标优化（默认会做一次 MMFF/UFF 优化以生成合理氢坐标）",
    )
    return parser.parse_args()


def add_hydrogens(mol, optimize_coords=True):
    """
    对分子加氢。

    参数:
        mol: RDKit Mol 对象
        optimize_coords: 是否在加氢后进行坐标优化

    返回:
        加氢后的 Mol 对象，失败返回 None
    """
    if mol is None:
        return None

    try:
        # 加氢（显式氢）
        mol_h = Chem.AddHs(mol, addCoords=True)

        # 可选：优化坐标，使氢原子位置更合理
        if optimize_coords and mol_h.GetNumConformers() > 0:
            try:
                # 优先使用 MMFF94
                if AllChem.MMFFHasAllMoleculeParams(mol_h):
                    AllChem.MMFFOptimizeMolecule(mol_h, maxIters=200)
                else:
                    AllChem.UFFOptimizeMolecule(mol_h, maxIters=200)
            except Exception as e:
                print(f"[警告] 坐标优化失败: {e}", file=sys.stderr)

        return mol_h
    except Exception as e:
        print(f"[错误] 加氢失败: {e}", file=sys.stderr)
        return None


def main():
    args = parse_args()

    # 使用 SDMolSupplier 读取，保留所有分子属性（如 score）
    supplier = Chem.SDMolSupplier(args.input, removeHs=False, sanitize=True)

    writer = Chem.SDWriter(args.output)

    n_total = 0
    n_written = 0
    n_failed = 0

    for idx, mol in enumerate(supplier):
        n_total += 1

        if mol is None:
            print(f"[警告] 第 {idx} 个分子读取失败，跳过。", file=sys.stderr)
            n_failed += 1
            continue

        mol_h = add_hydrogens(mol, optimize_coords=not args.no_coords)

        if mol_h is None:
            print(f"[警告] 第 {idx} 个分子加氢失败，跳过。", file=sys.stderr)
            n_failed += 1
            continue

        # 保留原分子名称（_Name）和属性
        if mol.HasProp("_Name"):
            mol_h.SetProp("_Name", mol.GetProp("_Name"))

        writer.write(mol_h)
        n_written += 1

    writer.close()

    print(f"完成: 共读取 {n_total} 个分子，成功写出 {n_written} 个，失败 {n_failed} 个。")
    print(f"输出文件: {args.output}")


if __name__ == "__main__":
    main()
