#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import logging
from pathlib import Path
import re

# 配置日志
def setup_logger():
    logger = logging.getLogger('conformer_workflow')
    logger.setLevel(logging.INFO)

    # 创建控制台处理器
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    # 创建格式化器
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)

    # 添加处理器到日志器
    logger.addHandler(ch)
    return logger

# 从energy文件中提取能量
def extract_energy(energy_file):
    try:
        with open(energy_file, 'r') as f:
            content = f.read()

        # 使用正则表达式匹配能量值
        match = re.search(r'\$energy\n\s*\d+\s+(-?\d+\.\d+)', content)
        if match:
            return match.group(1)
        else:
            return None
    except Exception as e:
        return None

# 更新XYZ文件中的能量信息
def update_xyz_energy(original_xyz, new_xyz, energy):
    try:
        with open(original_xyz, 'r') as f:
            lines = f.readlines()

        # 检查XYZ文件格式
        if len(lines) < 2:
            return False

        # 替换第二行
        lines[1] = f" energy: {energy} xtb: 6.7.1 with g-xtb energy\n"

        with open(new_xyz, 'w') as f:
            f.writelines(lines)

        return True
    except Exception as e:
        return False

def main():
    parser = argparse.ArgumentParser(description='Conformation clustering workflow with g-xtb and CREST')
    parser.add_argument('input_sdf', type=str, help='Input SDF file path')
    parser.add_argument('--charge', type=int, required=True, help='Formal charge of the molecule')
    args = parser.parse_args()

    logger = setup_logger()
    logger.info(f"Starting workflow for {args.input_sdf} with charge {args.charge}")

    # 创建目录结构
    dirs = ['1_conf', '2_conf_sp', '3_cluster']
    for d in dirs:
        Path(d).mkdir(exist_ok=True)
        logger.info(f"Created directory: {d}")

    # Step 1: 拆分SDF为单个XYZ构象
    logger.info("Step 1: Splitting SDF into individual conformers")
    try:
        result = subprocess.run(
            ['obabel', '-isdf', args.input_sdf, '-oxyz', '-O', '1_conf/CONF_.xyz', '-m'],
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            logger.error(f"OpenBabel failed: {result.stderr}")
            sys.exit(1)
        logger.info("SDF splitting completed successfully")
    except Exception as e:
        logger.error(f"Error running OpenBabel: {str(e)}")
        sys.exit(1)

    # 获取生成的构象数量
    conf_files = list(Path('1_conf').glob('CONF_*.xyz'))
    num_conformers = len(conf_files)
    if num_conformers == 0:
        logger.error("No conformers generated in 1_conf directory")
        sys.exit(1)
    logger.info(f"Number of conformers generated: {num_conformers}")

    # Step 2: 使用g-xtb计算单点能
    logger.info("Step 2: Calculating single-point energies with g-xtb")
    os.chdir('2_conf_sp')

    for i in range(1, num_conformers + 1):
        conf_dir = f"CONF_{i}"
        Path(conf_dir).mkdir(exist_ok=True)
        logger.info(f"Processing conformer {i}/{num_conformers}")

        # 创建电荷和UHF文件
        try:
            with open(f"{conf_dir}/.CHRG", 'w') as f:
                f.write(str(args.charge))
            with open(f"{conf_dir}/.UHF", 'w') as f:
                f.write("0")
        except Exception as e:
            logger.error(f"Error creating charge files for {conf_dir}: {str(e)}")
            continue

        # 创建XYZ文件的符号链接 - 修正路径错误
        src_xyz = f"../../1_conf/CONF_{i}.xyz"  # 修正为正确的相对路径
        dest_xyz = f"{conf_dir}/CONF_{i}.xyz"
        try:
            if os.path.exists(dest_xyz):
                os.remove(dest_xyz)
            os.symlink(src_xyz, dest_xyz)
            logger.debug(f"Created symlink: {dest_xyz} -> {src_xyz}")
        except Exception as e:
            logger.error(f"Error creating symlink for {conf_dir}: {str(e)}")
            continue

        # 运行g-xtb计算
        try:
            with open(f"{conf_dir}/CONF_{i}.log", 'w') as logfile:
                result = subprocess.run(
                    ['gxtb', '-c', f"CONF_{i}.xyz"],
                    cwd=conf_dir,
                    stdout=logfile,
                    stderr=subprocess.STDOUT
                )
            if result.returncode != 0:
                logger.warning(f"g-xtb calculation for CONF_{i} returned non-zero exit status")
            else:
                logger.info(f"g-xtb calculation completed for CONF_{i}")
        except Exception as e:
            logger.error(f"Error running g-xtb for CONF_{i}: {str(e)}")

    os.chdir('..')
    logger.info("Single-point energy calculations completed")

    # Step 3: 更新构象能量信息
    logger.info("Step 3: Updating conformer energies")
    ensemble_file = Path('2_conf_sp/ensemble_gxtb.xyz')
    if ensemble_file.exists():
        ensemble_file.unlink()

    for i in range(1, num_conformers + 1):
        conf_dir = Path(f"2_conf_sp/CONF_{i}")
        energy_file = conf_dir / 'energy'
        original_xyz = Path(f"1_conf/CONF_{i}.xyz")
        new_xyz = conf_dir / f"CONF_{i}_gxtb.xyz"

        # 提取能量
        energy = extract_energy(energy_file)
        if energy is None:
            logger.warning(f"Energy not found for CONF_{i}, skipping update")
            continue

        # 更新XYZ文件
        if not update_xyz_energy(original_xyz, new_xyz, energy):
            logger.warning(f"Failed to update XYZ for CONF_{i}")
            continue

        # 添加到系综文件
        try:
            with open(new_xyz, 'r') as src, open('2_conf_sp/ensemble_gxtb.xyz', 'a') as dst:
                dst.write(src.read())
            logger.info(f"Updated energy for CONF_{i}: {energy}")
        except Exception as e:
            logger.error(f"Error adding CONF_{i} to ensemble: {str(e)}")

    logger.info("Energy updates completed")

    # Step 4: 使用CREST聚类
    logger.info("Step 4: Clustering conformers with CREST")
    os.chdir('3_cluster')

    # 复制参考构象
    ref_xyz = Path("../2_conf_sp/CONF_1/CONF_1_gxtb.xyz")
    if not ref_xyz.exists():
        logger.warning("Reference conformer (CONF_1_gxtb.xyz) not found, using original CONF_1.xyz")
        ref_xyz = Path("../1_conf/CONF_1.xyz")

    try:
        subprocess.run(['cp', str(ref_xyz), '.'], check=True)
        logger.info(f"Copied reference conformer: {ref_xyz.name}")
    except Exception as e:
        logger.error(f"Error copying reference conformer: {str(e)}")
        sys.exit(1)

    # 运行CREST
    try:
        ensemble_path = Path("../2_conf_sp/ensemble_gxtb.xyz")
        if not ensemble_path.exists():
            logger.error(f"Ensemble file not found: {ensemble_path}")
            sys.exit(1)

        logger.info("Starting CREST clustering...")
        result = subprocess.run(
            ['crest', ref_xyz.name, '--cregen', str(ensemble_path),
             '-chrg', str(args.charge), '-rthr', '0.125', '-ewin', '10'],
            capture_output=True,
            text=True
        )

        # 记录CREST输出
        with open('crest.log', 'w') as f:
            f.write(result.stdout)
            if result.stderr:
                f.write("\n\nERRORS:\n")
                f.write(result.stderr)

        if result.returncode != 0:
            logger.warning(f"CREST returned non-zero exit status: {result.returncode}")
        else:
            logger.info("CREST clustering completed successfully")

        # 检查结果文件
        if Path('crest_ensemble.xyz').exists():
            logger.info("Final ensemble saved as crest_ensemble.xyz")
        else:
            logger.warning("CREST did not generate crest_ensemble.xyz")

    except Exception as e:
        logger.error(f"Error running CREST: {str(e)}")
        sys.exit(1)

    os.chdir('..')
    logger.info("Workflow completed")

if __name__ == "__main__":
    main()
