import argparse
import os
import shutil
import stat

def parse_arguments():
    parser = argparse.ArgumentParser(description="Schrödinger Multi-stage IFD Redocking Automation Script with Title Preservation.")
    parser.add_argument("-complex", required=True, help="输入蛋白-配体复合物结构文件路径 (例如: 4IVC_W1305.maegz)")
    parser.add_argument("-core", default=None, help="MCS 公共子结构核心约束文件路径 (例如: core.maegz，可选)")
    parser.add_argument("-ewin", type=float, default=20.0, help="配体环构象采样能量窗口 / Ringconfcut (默认: 20.0)")
    parser.add_argument("-max-n-pose", type=int, default=5, help="Glide 对接每配体保留的最大构象数 (默认: 5)")
    parser.add_argument("-ncpu", type=int, default=8, help="计算核心数 (默认: 8)")
    parser.add_argument("-host", default="localhost", help="计算节点主机名 (默认: localhost)")
    parser.add_argument("-jobname", default="ifd_redock", help="作业名称 (默认: ifd_redock)")
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    print("=== 1. 检查输入文件与内存结构拆分（含标题继承） ===")
    complex_path = os.path.abspath(args.complex)
    if not os.path.exists(complex_path):
        raise FileNotFoundError(f"Error: 找不到指定的复合物文件: {complex_path}")
    print(f"  -> 输入复合物: {complex_path}")

    receptor_file = f"{args.jobname}_rec.maegz"
    ligand_file = f"{args.jobname}_lig.maegz"

    try:
        from schrodinger import structure
        
        # 读取复合物
        reader = structure.StructureReader(complex_path)
        st_list = [st for st in reader]
        reader.close()
        
        if not st_list:
            raise ValueError(f"{complex_path} 中未找到有效结构")
            
        complex_st = st_list[0]
        original_title = complex_st.title if complex_st.title else os.path.splitext(os.path.basename(complex_path))[0]
        
        # 内部逻辑对齐 pv_convert：假定复合物中配体为最后一个分子
        lig_mol_num = len(complex_st.molecule)
        if lig_mol_num <= 1:
            raise ValueError("复合物中只包含单分子，无法通过受体/配体拆分！")
            
        lig_mol = complex_st.molecule[lig_mol_num]
        lig_st = lig_mol.extractStructure()
        
        # 完美恢复/继承配体标题
        lig_st.title = original_title
        
        # 提取受体（保留水分子、金属离子等）
        recep_atoms = set(complex_st.getAtomIndices()).difference(lig_mol.getAtomIndices())
        recep_st = complex_st.extract(recep_atoms)
        recep_st.title = f"{original_title}_receptor"
        
        # 保存待对接的纯配体文件（带正确 Title）
        lig_st.write(ligand_file)
        
        # 将配体复制一份赋予 Chain Z, Resnum 999 并合并进受体
        site_lig_st = lig_st.copy()
        for atom in site_lig_st.atom:
            atom.chain = 'Z'
            atom.resnum = 999
            atom.resname = 'LIG'
            
        recep_st.extend(site_lig_st)
        recep_st.write(receptor_file)
        
        print(f"  -> 待对接配体已保存 (Title: {lig_st.title}): {ligand_file}")
        print(f"  -> 包含定位标签(Z:999)的受体已保存 (Title: {recep_st.title}): {receptor_file}")

    except Exception as e:
        raise RuntimeError(f"处理复合物拆分时出错: {e}")

    print("\n=== 2. 处理核心约束文件 (Core) ===")
    core_block = ""
    if args.core:
        core_src_path = os.path.abspath(args.core)
        if not os.path.exists(core_src_path):
            raise FileNotFoundError(f"Error: 找不到指定的核心约束文件: {core_src_path}")
        
        local_core_file = f"{args.jobname}_reflig.maegz"
        shutil.copy(core_src_path, local_core_file)
        print(f"  -> 核心约束文件已同步到本地: {local_core_file} (已启用 MCS 约束)")
        
        core_block = f"""
  DOCKING_CORE_RESTRAIN YES
  DOCKING_CORE_FILTER YES
  DOCKING_USE_REF_LIGAND YES
  DOCKING_CORE_POS_MAX_RMSD 1.0
  DOCKING_CORE_DEFINITION mcssmarts
  DOCKING_REF_LIGAND_FILE {local_core_file}"""
    else:
        print("  -> 未指定核心约束文件 (-core)，将不使用核心约束。")
        core_block = "\n  DOCKING_CORE_RESTRAIN NO"

    print(f"  -> 环构象采样能量窗口 (Ringconfcut): {args.ewin}")
    print(f"  -> 最大构象数 (Poses per lig): {args.max_n_pose}")
    print(f"  -> 计算资源: {args.host} ({args.ncpu} 核心)")

    print("\n=== 3. 生成多阶段标准 IFD 配置文件 ===")
    in_filename = f"{args.jobname}.inp"
    
    ifd_in_content = f"""# Automatically generated Schrödinger Multi-stage IFD Redock input file

INPUT_FILE {receptor_file}

STAGE GLIDE_DOCKING2
  BINDING_SITE ligand Z:999
  INNERBOX 10.0
  OUTERBOX auto
  LIGAND_FILE {ligand_file}
  LIGANDS_TO_DOCK all
  GRIDGEN_RECEP_CCUT 0.25
  GRIDGEN_RECEP_VSCALE 0.50
  GRIDGEN_FORCEFIELD OPLS_2005
  DOCKING_PRECISION SP
  DOCKING_LIG_CCUT 0.15
  DOCKING_CV_CUTOFF 100.0
  DOCKING_LIG_VSCALE 0.50
  DOCKING_POSES_PER_LIG {args.max_n_pose}
  DOCKING_FORCEFIELD OPLS_2005
  DOCKING_RINGCONFCUT {args.ewin}
  DOCKING_AMIDE_MODE penal{core_block}

STAGE COMPILE_RESIDUE_LIST
  DISTANCE_CUTOFF 5.0

STAGE PRIME_REFINEMENT
  NUMBER_OF_PASSES 1
  USE_MEMBRANE no
  OPLS_VERSION OPLS_2005

STAGE SORT_AND_FILTER
  POSE_FILTER r_psp_Prime_Energy
  POSE_KEEP 30.0

STAGE SORT_AND_FILTER
  POSE_FILTER r_psp_Prime_Energy
  POSE_KEEP 20

STAGE GLIDE_DOCKING2
  BINDING_SITE ligand Z:999
  INNERBOX 10.0
  OUTERBOX auto
  LIGAND_FILE {ligand_file}
  LIGANDS_TO_DOCK self
  GRIDGEN_RECEP_CCUT 0.25
  GRIDGEN_RECEP_VSCALE 1.00
  GRIDGEN_FORCEFIELD OPLS_2005
  DOCKING_PRECISION SP
  DOCKING_LIG_CCUT 0.15
  DOCKING_CV_CUTOFF 0.0
  DOCKING_LIG_VSCALE 0.80
  DOCKING_POSES_PER_LIG 1
  DOCKING_FORCEFIELD OPLS_2005
  DOCKING_RINGCONFCUT {args.ewin}
  DOCKING_AMIDE_MODE penal{core_block}

STAGE SCORING
  SCORE_NAME r_psp_IFDScore
  TERM 1.0,r_i_glide_gscore,0
  TERM 0.05,r_psp_Prime_Energy,1
  REPORT_FILE report.csv
"""

    with open(in_filename, "w") as f:
        f.write(ifd_in_content.strip() + "\n")
    print(f"  -> 配置文件已生成: {in_filename}")

    print("\n=== 4. 自动生成 Shell 提交脚本 ===")
    sh_filename = f"run_{args.jobname}.sh"
    
    sh_content = f"""#!/bin/bash
#SBATCH -J {args.jobname}

"${{SCHRODINGER}}/ifd" -NGLIDECPU {args.ncpu} -NPRIMECPU {args.ncpu} {in_filename} -NOLOCAL -HOST {args.host} -SUBHOST {args.host} -TMPLAUNCHDIR
"""

    with open(sh_filename, "w") as f:
        f.write(sh_content)
    
    os.chmod(sh_filename, os.stat(sh_filename).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    
    print(f"  -> Shell 提交脚本已生成: {sh_filename}")
    print(f"\n标题修复完成！重新运行脚本即可在 Maestro 中完美显示配体标题：")
    print(f"  $SCHRODINGER/run sch_ifd_redock.py -complex ../4IVC_W1305.maegz -core ../core.maegz -ncpu 8")
    print(f"  ./{sh_filename}")

if __name__ == "__main__":
    main()
