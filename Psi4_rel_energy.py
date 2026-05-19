#!/usr/bin/env python
"""
sdf_boltzmann.py - 从 SDF 文件计算相对能量和玻尔兹曼权重（基于 Psi4_Energy）
用法: python sdf_boltzmann.py <input.sdf> [温度(K)]
输出: CSV 格式 (Title, Psi4_Rel_Energy (kcal/mol), Boltzmann_Weight)
"""

import sys
import math

# 气体常数 R = 0.00198720425864083 kcal/(mol·K)
R = 0.00198720425864083

def parse_sdf_entries(file_path):
    """解析 SDF，返回 [(title, psi4_energy), ...]"""
    entries = []
    with open(file_path, 'r') as f:
        content = f.read()
    blocks = content.strip().split('$$$$')
    for block in blocks:
        if not block.strip():
            continue
        lines = block.strip().split('\n')
        title = lines[0].strip()
        psi4_energy = None
        for i, line in enumerate(lines):
            # 支持 >  <Psi4_Energy (kcal/mol)> 或 > <Psi4_Energy (kcal/mol)>
            if line.startswith('>  <Psi4_Energy (kcal/mol)>') or line.startswith('> <Psi4_Energy (kcal/mol)>'):
                if i + 1 < len(lines):
                    try:
                        psi4_energy = float(lines[i+1].strip())
                    except ValueError:
                        pass
                break
        if psi4_energy is not None:
            entries.append((title, psi4_energy))
        else:
            print(f"警告: 构象 '{title}' 缺少 Psi4_Energy，已跳过", file=sys.stderr)
    return entries

def main():
    if len(sys.argv) < 2:
        print("用法: python sdf_boltzmann.py <input.sdf> [温度(K)]", file=sys.stderr)
        sys.exit(1)

    sdf_file = sys.argv[1]
    T = 298.15  # 默认室温
    if len(sys.argv) >= 3:
        try:
            T = float(sys.argv[2])
        except ValueError:
            print(f"警告: 温度 '{sys.argv[2]}' 无效，使用默认值 298.15 K", file=sys.stderr)

    entries = parse_sdf_entries(sdf_file)
    if not entries:
        print("错误: 未找到任何有效的 Psi4_Energy 数据", file=sys.stderr)
        sys.exit(1)

    # 提取能量列表，计算相对能量
    energies = [e for _, e in entries]
    min_energy = min(energies)
    rel_energies = [e - min_energy for e in energies]  # kcal/mol

    # 计算玻尔兹曼因子
    kT = R * T
    boltz_factors = [math.exp(-rel / kT) for rel in rel_energies]
    sum_factors = sum(boltz_factors)
    weights = [f / sum_factors for f in boltz_factors]

    # 输出 CSV (带标题行)
    print("Title,Psi4_Rel_Energy (kcal/mol),Boltzmann_Weight")
    for (title, _), rel, w in zip(entries, rel_energies, weights):
        # 格式化：相对能量保留 6 位小数，权重保留科学计数或普通小数（这里保留 8 位小数）
        print(f"{title},{rel:.6f},{w:.8f}")

if __name__ == "__main__":
    main()
