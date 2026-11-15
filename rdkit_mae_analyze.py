#!/usr/bin/env python
"""精准分析 .mae/.maegz 文件中单个 Entry 的成分及其数量。

用法:
    python rd_mae_analyze.py <input.mae|input.maegz>
    python rd_mae_analyze.py -h | --help
"""
import gzip
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, GetMolFrags
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolfiles


def print_help():
    """打印命令行帮助信息"""
    help_text = """
==================================================
 RDKit Schrödinger MAE/MAEGZ 成分分析工具
==================================================
用法:
    python rd_mae_analyze.py <input.mae|input.maegz>
    python rd_mae_analyze.py [-h | --help]

参数说明:
    <input.mae|input.maegz>    需要分析的薛定谔 Maestro 格式文件
    -h, --help                显示此帮助信息并退出
"""
    print(help_text.strip())


def open_mae_stream(path):
    """按扩展名打开输入文件，.maegz/.mae.gz 自动解压，返回二进制流。"""
    if path.endswith((".maegz", ".mae.gz")):
        return gzip.open(path, "rb")
    return open(path, "rb")


def analyze_first_entry(infile):
    try:
        with open_mae_stream(infile) as stream:
            suppl = rdmolfiles.MaeMolSupplier(stream, removeHs=False)
            # 获取第一个结构 (Entry)
            st = next(suppl)
    except StopIteration:
        sys.stderr.write(f"ERROR: 文件中未找到任何结构: {infile}\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"ERROR: 读取文件失败: {e}\n")
        sys.exit(1)

    title = st.GetProp("_Name") if st.HasProp("_Name") else "Unknown"
    print(f"==================================================")
    print(f" 文件分析报告: {infile}")
    print(f" 目标 Entry 标题: {title}")
    print(f"==================================================")

    # 初始化分类计数器
    counts = {
        "Protein (Chains/Fragments)": 0,
        "Ligand": 0,
        "Solvent (Water)": 0,
        "Metal Ion": 0,
        "Cofactor": 0,
        "Other": 0
    }
    
    details = []
    solvent_count = 0

    # 将复合物按非共价键断开，切分为独立的碎片（分子、离子、水）
    frags = GetMolFrags(st, asMols=True, sanitizeFrags=False)

    for idx, f in enumerate(frags):
        num_atoms = f.GetNumAtoms()
        mw = Descriptors.MolWt(f)
        
        # 安全计算分子式
        try:
            formula = rdMolDescriptors.CalcMolFormula(f)
        except Exception:
            formula = f"Atoms_{num_atoms}"

        # 1. 金属离子判断（单原子且为常见金属）
        if num_atoms == 1:
            atom = f.GetAtomWithIdx(0)
            symbol = atom.GetSymbol()
            if symbol in ["ZN", "MG", "CA", "FE", "MN", "CU", "NA", "K", "NI", "CO", "CD"]:
                counts["Metal Ion"] += 1
                details.append(f"  - [金属离子] 元素: {symbol}")
                continue

        # 2. 水分子判断（原子数小且含氧/氢）
        if num_atoms <= 3 and ("O" in formula or "H" in formula):
            solvent_count += 1
            continue

        # 3. 蛋白质残基片段判断（原子数较多，通常 > 50）
        if num_atoms > 50:
            counts["Protein (Chains/Fragments)"] += 1
            details.append(f"  - [蛋白质片段 {counts['Protein (Chains/Fragments)']}] 原子数: {num_atoms}")
        
        # 4. 辅助因子判断（分子量和原子数较大，如 HEME、ATP 等）
        elif mw > 350 and num_atoms > 25:
            counts["Cofactor"] += 1
            details.append(f"  - [辅助因子] 分子式: {formula}, 分子量: {mw:.1f}, 原子数: {num_atoms}")
        
        # 5. 常规小分子配体
        else:
            counts["Ligand"] += 1
            details.append(f"  - [配体] 分子式: {formula}, 分子量: {mw:.1f}, 原子数: {num_atoms}")

    # 汇总水分子
    if solvent_count > 0:
        counts["Solvent (Water)"] = solvent_count
        details.insert(0, f"  - [溶剂水分子] 总计: {solvent_count} 个")

    # 输出统计报告
    print("\n【成分数量统计】")
    for category, count in counts.items():
        print(f"  • {category}: {count} 个")

    print("\n【详细组成清单】")
    for d in details:
        print(d)
    print("==================================================")


def main():
    if len(sys.argv) != 2 or sys.argv[1] in ("-h", "--help"):
        print_help()
        sys.exit(0 if len(sys.argv) == 2 else 1)
    
    analyze_first_entry(sys.argv[1])


if __name__ == "__main__":
    main()
