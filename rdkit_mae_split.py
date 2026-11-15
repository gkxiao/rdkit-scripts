#!/usr/bin/env python
"""利用 RDKit 实现类似薛定谔 pv_convert.py 的受体与配体分离提取功能。

用法:
    python rd_mae_split.py <input.maegz> --mode ligand   # 仅提取并保存配体 (SDF)
    python rd_mae_split.py <input.maegz> --mode receptor # 仅提取并保存受体 (PDB/MAE)
    python rd_mae_split.py <input.maegz> --mode all      # 同时导出受体和配体
"""
import argparse
import gzip
import os
from rdkit import Chem
from rdkit.Chem import Descriptors, GetMolFrags
from rdkit.Chem import rdmolfiles


def open_mae_stream(path):
    if path.endswith((".maegz", ".mae.gz")):
        return gzip.open(path, "rb")
    return open(path, "rb")


def split_receptor_ligand(infile, output_prefix, mode):
    try:
        with open_mae_stream(infile) as stream:
            suppl = rdmolfiles.MaeMolSupplier(stream, removeHs=False)
            st = next(suppl)
    except Exception as e:
        print(f"ERROR: 读取文件失败: {e}")
        return

    title = st.GetProp("_Name") if st.HasProp("_Name") else "complex"
    print(f"正在处理构象: {title}")

    # 将复合物按非共价键断开，切分为独立碎片
    frags = GetMolFrags(st, asMols=True, sanitizeFrags=False)

    protein_mols = []
    ligand_mols = []
    other_mols = []

    for f in frags:
        num_atoms = f.GetNumAtoms()
        mw = Descriptors.MolWt(f)
        
        # 1. 判断是否为金属离子或水（跳过或归为受体环境视需求而定）
        if num_atoms == 1 or (num_atoms <= 3 and mw < 30):
            continue  # 忽略水和金属离子，或放入受体中
            
        # 2. 蛋白质片段判定（原子数多于 50）
        if num_atoms > 50:
            protein_mols.append(f)
        # 3. 配体判定（分子量和原子数适中的有机小分子）
        elif mw > 100 and num_atoms > 5:
            ligand_mols.append(f)
        else:
            other_mols.append(f)

    print(f" -> 识别到蛋白片段数: {len(protein_mols)}")
    print(f" -> 识别到配体分子数: {len(ligand_mols)}")

    # 导出配体模式
    if mode in ["ligand", "all"] and ligand_mols:
        ligand_out = f"{output_prefix}_ligand.sdf"
        writer = Chem.SDWriter(ligand_out)
        for lig in ligand_mols:
            writer.write(lig)
        writer.close()
        print(f" [成功] 配体已保存至: {ligand_out}")

    # 导出受体模式（将所有蛋白片段合并为一个大分子或组合输出）
    if mode in ["receptor", "all"] and protein_mols:
        receptor_out = f"{output_prefix}_protein.sdf" # 或者保存为 PDB
        # RDKit 可以通过 CombineMols 将多条链合并为一个受体对象
        receptor_mol = protein_mols[0]
        for p in protein_mols[1:]:
            receptor_mol = Chem.CombineMols(receptor_mol, p)
            
        # 写入文件 (支持保存为 SDF 或利用 Chem.MolToPDBFile 保存为 PDB)
        pdb_out = f"{output_prefix}_receptor.pdb"
        Chem.MolToPDBFile(receptor_mol, pdb_out)
        print(f" [成功] 受体已保存至: {pdb_out}")


def main():
    parser = argparse.ArgumentParser(description="使用 RDKit 提取 MAE/MAEGZ 中的受体与配体")
    parser.add_argument("input", help="输入的 .mae 或 .maegz 文件路径")
    parser.add_argument("--mode", choices=["ligand", "receptor", "all"], default="all", help="提取模式: ligand, receptor, 或 all")
    parser.add_argument("--prefix", default="output", help="输出文件的前缀")
    
    args = parser.parse_args()
    split_receptor_ligand(args.input, args.prefix, args.mode)


if __name__ == "__main__":
    main()
