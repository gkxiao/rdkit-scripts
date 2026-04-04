#!/usr/bin/env python
import sys
import argparse
from rdkit import Chem
from rdkit.Chem import rdFreeSASA

def get_residue_delta_sasa(pdb_file, chain_id, res_num):
    """
    计算指定残基的 Delta SASA (A^2)
    """
    # 1. 加载分子 (必须保留氢原子以获得准确的物理表面)
    # RDKit 的 PDB 解析器会自动处理 PDBResidueInfo
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    if not mol:
        print(f"Error: 无法读取 PDB 文件: {pdb_file}")
        return None

    # 内部辅助函数：计算特定状态下的残基 SASA 贡献
    def calculate_sasa(input_mol, target_chain, target_res):
        ptable = Chem.GetPeriodicTable()
        # 使用 RDKit 周期表内置的 Van der Waals 半径
        radii = [ptable.GetRvdw(atom.GetAtomicNum()) for atom in input_mol.GetAtoms()]

        # 执行 SASA 计算 (默认使用 1.4A 探针半径)
        rdFreeSASA.CalcSASA(input_mol, radii)

        res_sum = 0.0
        found = False
        for atom in input_mol.GetAtoms():
            info = atom.GetPDBResidueInfo()
            if info:
                # 匹配 Chain ID 和 Residue Number (修正后的方法名)
                if info.GetChainId().strip() == target_chain and info.GetResidueNumber() == target_res:
                    # 从原子属性中提取该原子的 SASA 贡献值
                    res_sum += float(atom.GetProp('SASA'))
                    found = True
        return res_sum if found else None

    # 2. 计算复合物 (Complex) 状态下的 SASA
    sasa_complex = calculate_sasa(mol, chain_id, res_num)

    # 3. 构建单体 (Monomer) 并计算
    # 逻辑：从分子中删除所有不属于目标链的原子
    edit_mol = Chem.EditableMol(mol)
    atoms_to_del = [a.GetIdx() for a in mol.GetAtoms()
                    if a.GetPDBResidueInfo() and a.GetPDBResidueInfo().GetChainId().strip() != chain_id]

    # 逆序删除以保持索引一致
    for idx in sorted(atoms_to_del, reverse=True):
        edit_mol.RemoveAtom(idx)

    monomer_mol = edit_mol.GetMol()
    sasa_monomer = calculate_sasa(monomer_mol, chain_id, res_num)

    if sasa_complex is None or sasa_monomer is None:
        print(f"Error: 在链 {chain_id} 中未找到残基编号 {res_num}")
        return None

    return {
        "residue": f"{chain_id}:{res_num}",
        "monomer_sasa": sasa_monomer,
        "complex_sasa": sasa_complex,
        "delta_sasa": sasa_monomer - sasa_complex
    }

def main():
    parser = argparse.ArgumentParser(description="计算 PDB 中指定残基的 ΔSASA (绝对值)")
    parser.add_argument("pdb", help="输入的 PDB 文件路径")
    parser.add_argument("chain", help="目标链 ID (例如: A)")
    parser.add_argument("resno", type=int, help="目标残基编号 (例如: 89)")

    args = parser.parse_args()

    result = get_residue_delta_sasa(args.pdb, args.chain, args.resno)

    if result:
        print(f"\n[ SASA Analysis Result ]")
        print(f"Residue:      {result['residue']}")
        print(f"Monomer SASA: {result['monomer_sasa']:.3f} Å²")
        print(f"Complex SASA: {result['complex_sasa']:.3f} Å²")
        print(f"ΔSASA:        {result['delta_sasa']:.3f} Å²")
        print("-" * 25)

if __name__ == "__main__":
    main()
(openfe_env) gkxiao@master:/public/gkxiao/work/beigene/sasa$ python sasa.py -h
python: can't open file '/public/gkxiao/work/beigene/sasa/sasa.py': [Errno 2] No such file or directory
(openfe_env) gkxiao@master:/public/gkxiao/work/beigene/sasa$ python sasa_calc.py -h
usage: sasa_calc.py [-h] pdb chain resno

计算 PDB 中指定残基的 ΔSASA (绝对值)

positional arguments:
  pdb         输入的 PDB 文件路径
  chain       目标链 ID (例如: A)
  resno       目标残基编号 (例如: 89)

options:
  -h, --help  show this help message and exit
