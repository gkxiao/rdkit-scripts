import argparse
import os
import stat

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate Schrödinger Multi-stage IFD input configuration and submission script.")
    parser.add_argument("-receptor", required=True, help="受体结构文件路径 (例如: InducedFit_jak1_rec.mae)")
    parser.add_argument("-site", required=True, help="结合位点参考残基或文件 (例如: A:1201 或 参考配体)")
    parser.add_argument("-dbase", required=True, help="待对接的配体数据库文件路径 (.sdf 或 .maegz)")
    parser.add_argument("-core", default=None, help="MCS 公共子结构核心约束文件路径 (例如: core.maegz，可选)")
    parser.add_argument("-ewin", type=float, default=20.0, help="配体环构象采样能量窗口 / Ringconfcut (默认: 20.0)")
    parser.add_argument("-max-n-pose", type=int, default=5, help="Glide 对接每配体保留的最大构象数 (默认: 5)")
    parser.add_argument("-ncpu", type=int, default=8, help="计算核心数 (默认: 8)")
    parser.add_argument("-host", default="localhost", help="计算节点主机名 (默认: localhost)")
    parser.add_argument("-jobname", default="ifd_run", help="作业名称 (默认: ifd_run)")
    return parser.parse_args()

def main():
    args = parse_arguments()
    
    print("=== 1. 检查输入文件与参数 ===")
    for filepath, desc in [
        (args.receptor, "受体文件"), 
        (args.dbase, "对接数据库")
    ]:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Error: 找不到指定的{desc}: {filepath}")
        print(f"  -> {desc}: {filepath}")

    # 处理核心约束相关参数块
    core_block = ""
    if args.core:
        if not os.path.exists(args.core):
            raise FileNotFoundError(f"Error: 找不到指定的核心约束文件: {args.core}")
        print(f"  -> 核心约束文件: {args.core} (已启用 MCS 约束)")
        core_block = f"""
  DOCKING_CORE_RESTRAIN YES
  DOCKING_CORE_FILTER YES
  DOCKING_USE_REF_LIGAND YES
  DOCKING_CORE_POS_MAX_RMSD 1.0
  DOCKING_CORE_DEFINITION mcssmarts
  DOCKING_REF_LIGAND_FILE {args.core}"""
    else:
        print("  -> 未指定核心约束文件 (-core)，将不使用核心约束。")
        core_block = "\n  DOCKING_CORE_RESTRAIN NO"

    print(f"  -> 环构象采样能量窗口 (Ringconfcut): {args.ewin}")
    print(f"  -> 最大构象数 (Poses per lig): {args.max_n_pose}")
    print(f"  -> 计算资源: {args.host} ({args.ncpu} 核心)")

    print("\n=== 2. 生成多阶段 IFD 任务输入文件 (ifd_job.inp) ===")
    
    in_filename = f"{args.jobname}.inp"
    
    # 构建完整的多阶段 multisim 输入内容
    ifd_in_content = f"""# Automatically generated Schrödinger Multi-stage IFD input file

INPUT_FILE {args.receptor}

STAGE GLIDE_DOCKING2
  BINDING_SITE ligand {args.site}
  INNERBOX 10.0
  OUTERBOX auto
  LIGAND_FILE {args.dbase}
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
  BINDING_SITE ligand {args.site}
  INNERBOX 10.0
  OUTERBOX auto
  LIGAND_FILE {args.dbase}
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

    print("\n=== 3. 自动生成 Shell 提交脚本 (run_ifd.sh) ===")
    sh_filename = f"run_{args.jobname}.sh"
    
    # 按照标准的 Schrödinger Master 格式编写 Shell 脚本
    sh_content = f"""#!/bin/bash
#SBATCH -J {args.jobname}

"${{SCHRODINGER}}/ifd" -NGLIDECPU {args.ncpu} -NPRIMECPU {args.ncpu} {in_filename} -NOLOCAL -HOST {args.host} -SUBHOST {args.host} -TMPLAUNCHDIR
"""

    with open(sh_filename, "w") as f:
        f.write(sh_content)
    
    # 赋予 .sh 文件可执行权限
    os.chmod(sh_filename, os.stat(sh_filename).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    
    print(f"  -> Shell 提交脚本已生成并赋予执行权限: {sh_filename}")
    print("\n您可以直接通过以下命令运行任务：")
    print(f"  ./{sh_filename}")

if __name__ == "__main__":
    main()
