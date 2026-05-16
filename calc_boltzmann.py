#!/usr/bin/env python3
from rdkit import Chem
import numpy as np
import sys
import argparse

def read_sdf_energies(sdf_file):
    """读取SDF文件中所有分子的Psi4_Energy (kcal/mol)"""
    energies = []
    mols = []
    suppl = Chem.SDMolSupplier(sdf_file)
    
    for mol in suppl:
        if mol is None:
            energies.append(np.nan)
            mols.append(None)
            continue
        # 获取Psi4_Energy性质
        if mol.HasProp('Psi4_Energy (kcal/mol)'):
            energy_str = mol.GetProp('Psi4_Energy (kcal/mol)')
            try:
                energy = float(energy_str)
                energies.append(energy)
            except ValueError:
                print(f"警告: 无法转换能量值 '{energy_str}' 为数字，设为NaN")
                energies.append(np.nan)
        else:
            print(f"警告: 分子缺少 'Psi4_Energy (kcal/mol)' 性质，设为NaN")
            energies.append(np.nan)
        mols.append(mol)
    
    return mols, energies

def compute_relative_energies(energies):
    """计算相对能量，忽略NaN"""
    valid_energies = [e for e in energies if not np.isnan(e)]
    if not valid_energies:
        raise ValueError("没有有效的能量数据")
    
    min_energy = min(valid_energies)
    rel_energies = []
    for e in energies:
        if np.isnan(e):
            rel_energies.append(np.nan)
        else:
            rel_energies.append(e - min_energy)
    return rel_energies, min_energy

def compute_boltzmann_weights(rel_energies, temperature=298.15):
    """
    计算Boltzmann权重 (EPOP)
    公式: weight = exp(-E_rel / (R * T))，然后归一化
    R = 0.00198720425864083 kcal/(mol·K) (理想气体常数)
    """
    R = 0.00198720425864083  # kcal/(mol·K)
    kT = R * temperature
    
    # 计算指数权重（只对有效相对能量）
    weights = []
    valid_indices = []
    for i, e in enumerate(rel_energies):
        if np.isnan(e):
            weights.append(np.nan)
        else:
            # E_rel 单位是 kcal/mol
            w = np.exp(-e / kT)
            weights.append(w)
            valid_indices.append(i)
    
    # 归一化（只对有效数据）
    if valid_indices:
        valid_weights = [weights[i] for i in valid_indices]
        sum_weights = sum(valid_weights)
        if sum_weights > 0:
            for i in valid_indices:
                weights[i] = weights[i] / sum_weights
        else:
            # 如果所有权重为0（能量极高），平均分配
            for i in valid_indices:
                weights[i] = 1.0 / len(valid_indices)
    
    return weights

def write_sdf_with_properties(input_file, output_file, rel_energies, boltzmann_weights):
    """将相对能量和EPOP写回SDF文件"""
    # 重新读取分子
    suppl = Chem.SDMolSupplier(input_file)
    writer = Chem.SDWriter(output_file)
    
    for idx, mol in enumerate(suppl):
        if mol is None:
            print(f"警告: 跳过构象 {idx+1}，无法读取")
            continue
        
        # 添加新性质
        rel_e = rel_energies[idx]
        bw = boltzmann_weights[idx]
        
        if not np.isnan(rel_e):
            mol.SetProp('Psi4_Rel_Energy (kcal/mol)', f"{rel_e:.6f}")
        else:
            mol.SetProp('Psi4_Rel_Energy (kcal/mol)', 'NaN')
        
        if not np.isnan(bw):
            mol.SetProp('EPOP', f"{bw:.10f}")
        else:
            mol.SetProp('EPOP', 'NaN')
        
        writer.write(mol)
    
    writer.close()

def main():
    parser = argparse.ArgumentParser(
        description='计算SDF文件中构象的相对能量和Boltzmann权重(EPOP)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  %(prog)s input.sdf                    # 输出为 input_with_rel_energy.sdf
  %(prog)s input.sdf -o output.sdf     # 指定输出文件名
  %(prog)s input.sdf -t 300.0          # 使用300K温度

注意:
  输入文件必须包含 'Psi4_Energy (kcal/mol)' 性质
  能量单位假设为 kcal/mol
        """
    )
    
    parser.add_argument('input', 
                       help='输入SDF文件路径')
    parser.add_argument('-o', '--output', 
                       default=None,
                       help='输出SDF文件路径 (默认: 输入文件名_with_rel_energy.sdf)')
    parser.add_argument('-t', '--temperature', 
                       type=float, 
                       default=298.15,
                       help='温度 (K)，默认: 298.15 K')
    parser.add_argument('-v', '--verbose', 
                       action='store_true',
                       help='显示详细信息')
    
    args = parser.parse_args()
    
    # 设置输出文件名
    if args.output is None:
        if args.input.endswith('.sdf'):
            args.output = args.input.replace('.sdf', '_with_rel_energy.sdf')
        else:
            args.output = args.input + '_with_rel_energy.sdf'
    
    # 读取数据
    if args.verbose:
        print(f"读取文件: {args.input}")
    try:
        mols, energies = read_sdf_energies(args.input)
    except OSError as e:
        print(f"错误: 无法读取文件 '{args.input}'")
        print(f"详细信息: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"错误: 读取文件时发生异常: {e}")
        sys.exit(1)
    
    # 计算相对能量
    try:
        rel_energies, min_energy = compute_relative_energies(energies)
    except ValueError as e:
        print(f"错误: {e}")
        print("请确保输入文件包含有效的 'Psi4_Energy (kcal/mol)' 数据")
        sys.exit(1)
    
    # 计算Boltzmann权重
    boltzmann_weights = compute_boltzmann_weights(rel_energies, args.temperature)
    
    # 写入新文件
    try:
        write_sdf_with_properties(args.input, args.output, rel_energies, boltzmann_weights)
    except Exception as e:
        print(f"错误: 写入文件时发生异常: {e}")
        sys.exit(1)
    
    # 输出统计信息
    valid_indices = [i for i, e in enumerate(rel_energies) if not np.isnan(e)]
    valid_count = len(valid_indices)
    
    print(f"完成！")
    print(f"  最小能量: {min_energy:.6f} kcal/mol")
    print(f"  温度: {args.temperature:.2f} K")
    print(f"  成功处理: {valid_count}/{len(mols)} 个构象")
    print(f"  输出文件: {args.output}")
    
    if args.verbose and valid_count > 0:
        print("\n详细结果:")
        print(f"{'构象':<8} {'相对能量(kcal/mol)':<25} {'EPOP':<15}")
        print("-" * 50)
        for i in valid_indices[:10]:  # 只显示前10个
            print(f"{i+1:<8} {rel_energies[i]:<25.6f} {boltzmann_weights[i]:<15.6f}")
        if valid_count > 10:
            print(f"... 共{valid_count}个构象，仅显示前10个")

if __name__ == "__main__":
    main()
