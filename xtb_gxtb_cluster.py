#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gxtb_cluster.py
===============
SDF  ->  逐构象拆分  ->  g-xTB 单点能  ->  CREST(cregen) 去重聚类

===============================================================================
【重要变更】g-xTB 已不再是独立可执行文件，而是并入 xtb（经 tblite 库调用）
===============================================================================

    旧（已失效）:
        gxtb -c CONF_1.xyz
        产出 Turbomole 风格的 `energy` 文件，需解析其中的 `$energy` 块

    新:
        xtb CONF_1.xyz -c 0 -u 0 --sp --gxtb --cosmo water
        不再产出 `energy` 文件！总能量只打印到标准输出，必须从 stdout 日志读取

能量提取方式（唯一来源：stdout 日志）
    锚定末尾的盒线，每次运行只出现一次，且为全精度：
        | TOTAL ENERGY            -1247.539215742853 Eh   |
    等价命令行验证：  grep 'TOTAL ENERGY' gxtb_cosmo.log

    为什么不匹配泛化的 "total energy"：日志里还有
        行 137: total energy   -1.2475256933249E+03 Eh   <- 截断到 14 位有效数字
        行 126: cycle  total energy  energy error  density error  <- 表头，无单位
    前者精度不足（差 ~1e-10 Eh），后者根本不是能量值。故优先取 TOTAL ENERGY 盒线，
    仅在盒线缺失时才逐级回退到 `:: total energy ... Eh ::` 摘要行。

    ★ 该能量**已包含 ddCOSMO 水相溶剂化贡献**（g-xTB + --cosmo 只打印一个总能量，
      日志中不存在单独的 Gsolv 项）。因此不得再额外叠加任何溶剂化能修正。

===============================================================================
默认计算设置（命令行可覆盖）
===============================================================================
    --gxtb          g-xTB 方法（取代旧 gxtb 独立程序）
    --sp            单点能（不做几何优化，本工作流只做单点能聚类）
    --cosmo water   水相 ddCOSMO 隐式溶剂化
    -u 0            基态，0 个未成对电子

本脚本只做单点能，不含任何几何优化步骤。输入坐标按原样使用。

===============================================================================
实现要点（踩过的坑）
===============================================================================
* 每个构象放在独立子目录 2_conf_sp/CONF_i/ 下运行 xtb。
  xtb 默认会写/读 xtbrestart 并在下次运行打印 "Wavefunction read from restart
  file"（见 gxtb_gas.log 第 114 行）。若所有构象挤在同一目录，第 2 个构象就会
  读入第 1 个构象的波函数，结果不可信。目录隔离是最省事的正确做法。

* --gxtb 走 tblite 库；帮助文档里"solvation is not yet implemented within
  tblite"说的是 --tblite 选项，与 --gxtb 无关。gxtb_cosmo.log 已实证
  `--sp --gxtb --cosmo water` 正常工作，且能量与气相不同（相对更低），
  说明溶剂化确实生效。

* 用 shutil.copy2 复制坐标而非 os.symlink：Windows 下创建符号链接常因权限
  失败，且 xtb 会在同目录写一堆中间文件。

* 传给 CREST 的能量必须是 Hartree，写在 XYZ 第二行注释里，格式 "energy: <x>"。

依赖：xtb (>= 6.7，含 tblite/g-xTB 支持)、obabel、crest
"""

import argparse
import logging
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# 1 Hartree = 627.5094740631 kcal/mol
HARTREE2KCAL = 627.5094740631

# --------------------------------------------------------------------------- #
# 能量行正则
# --------------------------------------------------------------------------- #
# 主匹配：末尾的盒线，全精度，每次运行仅出现一次
#   | TOTAL ENERGY            -1247.539215742853 Eh   |
# 用 "Eh" 收尾做锚点，并允许科学计数法（个别版本/大分子可能溢出到 E 记法）
RE_TOTAL_ENERGY_BOX = re.compile(
    r"\|\s*TOTAL\s+ENERGY\s+(-?\d+\.\d+(?:[eE][-+]?\d+)?)\s*Eh", re.IGNORECASE
)

# 回退 1：SUMMARY 段里的摘要行，同为全精度
#   :: total energy           -1247.539215742853 Eh    ::
RE_TOTAL_ENERGY_SUMMARY = re.compile(
    r"::\s*total\s+energy\s+(-?\d+\.\d+(?:[eE][-+]?\d+)?)\s*Eh", re.IGNORECASE
)

# 回退 2：任何 "total energy ... Eh"（可能命中截断值，仅作最后保险）
RE_TOTAL_ENERGY_LOOSE = re.compile(
    r"total\s+energy\s+(-?\d+\.\d+(?:[eE][-+]?\d+)?)\s*Eh", re.IGNORECASE
)

# g-xTB 正常收尾标记
TERMINATION_MARKERS = ("normal termination", "finished run")

# 元素符号 -> 原子序数（1-118），仅用于电子数奇偶性校验
ELEMENTS = [
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]
SYMBOL2Z = {sym: z for z, sym in enumerate(ELEMENTS) if sym}


# --------------------------------------------------------------------------- #
# 日志
# --------------------------------------------------------------------------- #
def setup_logger():
    logger = logging.getLogger("conformer_workflow")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logger.addHandler(ch)
    return logger


# --------------------------------------------------------------------------- #
# 环境检查
# --------------------------------------------------------------------------- #
def probe_tool(name, logger):
    """确认可执行文件在 PATH 中，返回其路径；缺失则返回 None。"""
    path = shutil.which(name)
    if path is None:
        logger.error(f"未找到可执行文件: {name}（请确认已安装并在 PATH 中）")
    else:
        logger.info(f"{name}: {path}")
    return path


def check_gxtb_support(logger):
    """探测当前 xtb 是否认识 --gxtb（旧版 xtb 没有此选项）。"""
    try:
        proc = subprocess.run(
            ["xtb", "--help"], capture_output=True, text=True, timeout=60
        )
        blob = (proc.stdout or "") + (proc.stderr or "")
        if "gxtb" in blob.lower():
            logger.info("xtb 支持 --gxtb（g-xTB 经 tblite 调用）")
            return True
        logger.error(
            "当前 xtb 的 --help 中未见 gxtb 选项。\n"
            "  g-xTB 需要 xtb >= 6.7 且编译时链接 tblite。\n"
            "  请升级 xtb，或确认这份 xtb 是否为新版。"
        )
        return False
    except Exception as e:
        logger.warning(f"无法探测 xtb 是否支持 --gxtb（{e}），继续尝试")
        return True


def count_electrons(xyz_path):
    """从 XYZ 统计总电子数所需的核电荷和；解析失败返回 None。"""
    try:
        lines = Path(xyz_path).read_text(encoding="utf-8", errors="ignore").splitlines()
        natoms = int(lines[0].split()[0])
        zsum = 0
        for line in lines[2:2 + natoms]:
            sym = line.split()[0]
            sym = sym[:1].upper() + sym[1:].lower() if len(sym) > 1 else sym.upper()
            if sym not in SYMBOL2Z:
                return None
            zsum += SYMBOL2Z[sym]
        return zsum
    except Exception:
        return None


def validate_uhf(xyz_path, charge, uhf, logger):
    """校验 charge 与 uhf 的奇偶性是否自洽（电子数必须与未成对电子同奇偶）。"""
    zsum = count_electrons(xyz_path)
    if zsum is None:
        logger.debug("无法解析元素组成，跳过电子数校验")
        return True
    nelec = zsum - charge
    if nelec < 0:
        logger.error(f"电荷 {charge} 过大：核电荷和 {zsum}，电子数 {nelec} < 0")
        return False
    if nelec % 2 != uhf % 2:
        logger.error(
            f"电荷/自旋设置不自洽：核电荷和 {zsum} - 电荷 {charge} = {nelec} 个电子"
            f"（{'偶' if nelec % 2 == 0 else '奇'}数），"
            f"但 -u {uhf} 要求{'偶' if uhf % 2 == 0 else '奇'}数电子。"
        )
        return False
    logger.info(f"电子数校验通过：Z_sum={zsum}, charge={charge}, N_elec={nelec}, "
                f"unpaired={uhf}")
    return True


# --------------------------------------------------------------------------- #
# 能量解析（新版 xtb --gxtb，无 energy 文件，只能读 stdout）
# --------------------------------------------------------------------------- #
def extract_energy(log_file):
    """从 xtb 的 stdout 日志中提取总能量（Hartree）。失败返回 None。

    g-xTB 经 xtb 调用后不再写 Turbomole 的 `energy` 文件，能量只出现在标准输出。
    提取顺序（严格的精度优先级）：

        1) `| TOTAL ENERGY   <float> Eh |`   末尾盒线，全精度，唯一权威值
        2) `:: total energy   <float> Eh ::` SUMMARY 摘要行，全精度
        3) 任意 `total energy <float> Eh`     兜底（可能是截断值，仅防格式变动）

    返回值为总能量，**已包含 --cosmo 指定的溶剂化贡献**。
    """
    try:
        content = Path(log_file).read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return None

    if not content.strip():
        return None

    for pattern in (RE_TOTAL_ENERGY_BOX,
                    RE_TOTAL_ENERGY_SUMMARY,
                    RE_TOTAL_ENERGY_LOOSE):
        matches = pattern.findall(content)
        if matches:
            try:
                # 取最后一次匹配：单点计算里最后打印的才是最终结果
                return float(matches[-1])
            except ValueError:
                continue

    return None


def extract_energy_detail(log_file):
    """返回 (能量, 命中来源)。来源用于日志告知，便于排查版本差异。

    来源取值: 'TOTAL ENERGY box' / 'SUMMARY line' / 'loose match' / None
    """
    try:
        content = Path(log_file).read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return None, None

    labeled = (
        (RE_TOTAL_ENERGY_BOX, "TOTAL ENERGY box"),
        (RE_TOTAL_ENERGY_SUMMARY, "SUMMARY line"),
        (RE_TOTAL_ENERGY_LOOSE, "loose match"),
    )
    for pattern, label in labeled:
        matches = pattern.findall(content)
        if matches:
            try:
                return float(matches[-1]), label
            except ValueError:
                continue
    return None, None


def check_termination(log_file):
    """确认 xtb 正常收尾（日志尾部有 normal termination 或 finished run）。"""
    try:
        tail = Path(log_file).read_text(encoding="utf-8", errors="ignore")[-4000:]
    except Exception:
        return False
    low = tail.lower()
    return any(marker in low for marker in TERMINATION_MARKERS)


def tail_of(path, n=12):
    try:
        lines = Path(path).read_text(encoding="utf-8", errors="ignore").splitlines()
        return "\n".join("      " + s for s in lines[-n:])
    except Exception:
        return "      <无法读取日志>"


# --------------------------------------------------------------------------- #
# XYZ 注释行写入能量（供 CREST 读取，单位必须为 Eh）
# --------------------------------------------------------------------------- #
def update_xyz_energy(original_xyz, new_xyz, energy_eh, solvent_tag):
    """复制坐标并把能量写进第二行注释。CREST 认 "energy: <float>" 格式。"""
    try:
        with open(original_xyz, "r") as f:
            lines = f.readlines()
        if len(lines) < 2:
            return False
        lines[1] = (
            f" energy: {energy_eh:.12f} | xtb --gxtb --sp | solvent: {solvent_tag}\n"
        )
        Path(new_xyz).parent.mkdir(parents=True, exist_ok=True)
        with open(new_xyz, "w") as f:
            f.writelines(lines)
        return True
    except Exception:
        return False


# --------------------------------------------------------------------------- #
# xtb 单点能调用
# --------------------------------------------------------------------------- #
def build_sp_cmd(geom_name, charge, uhf, solvent, threads=None):
    """组装 g-xTB 单点能命令。

    固定包含: --gxtb（g-xTB 方法）、--sp（单点）、-u <uhf>（基态）
    solvent 非空时追加 --cosmo <solvent>（水相）。

    刻意不加 --json：能量直接从 stdout 的 TOTAL ENERGY 行读取，
    少写一个中间文件，命令形式与手工验证时完全一致。
    """
    cmd = [
        "xtb", str(geom_name),
        "-c", str(charge),
        "-u", str(uhf),
        "--sp",
        "--gxtb",
    ]
    if solvent:
        cmd += ["--cosmo", solvent]
    if threads:
        cmd += ["-P", str(threads)]
    return cmd


def run_xtb(cmd, cwd, log_path):
    """执行 xtb，stdout/stderr 合并写入 log_path。返回 returncode。"""
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "w", encoding="utf-8", errors="replace") as logfile:
        proc = subprocess.run(
            cmd, cwd=str(cwd), stdout=logfile, stderr=subprocess.STDOUT
        )
    return proc.returncode


# --------------------------------------------------------------------------- #
# 单个构象：g-xTB 单点能
# --------------------------------------------------------------------------- #
def process_conformer(idx, args, logger, solvent_tag):
    """对第 idx 个构象做 g-xTB 单点能。

    返回 dict: {index, energy(Eh), xyz(带能量注释的路径), ok, seconds}
    """
    conf_dir = Path("2_conf_sp") / f"CONF_{idx}"
    # 目录隔离：避免 xtb 读到上一个构象的 xtbrestart
    if conf_dir.exists():
        shutil.rmtree(conf_dir, ignore_errors=True)
    conf_dir.mkdir(parents=True, exist_ok=True)

    src_xyz = Path("1_conf") / f"CONF_{idx}.xyz"
    if not src_xyz.is_file():
        logger.error(f"CONF_{idx}: 源坐标不存在 {src_xyz}")
        return {"index": idx, "energy": None, "xyz": None, "ok": False, "seconds": 0}

    work_xyz = conf_dir / f"CONF_{idx}.xyz"
    try:
        shutil.copy2(src_xyz, work_xyz)
    except Exception as e:
        logger.error(f"CONF_{idx}: 复制坐标失败 - {e}")
        return {"index": idx, "energy": None, "xyz": None, "ok": False, "seconds": 0}

    sp_log = conf_dir / f"CONF_{idx}_gxtb.log"
    cmd = build_sp_cmd(work_xyz.name, args.charge, args.uhf, args.solvent,
                       args.xtb_threads)
    logger.info(f"CONF_{idx}/{args.total}: {' '.join(cmd)}")

    t0 = time.time()
    try:
        rc = run_xtb(cmd, conf_dir, sp_log)
    except Exception as e:
        logger.error(f"CONF_{idx}: 调用 xtb 失败 - {e}")
        return {"index": idx, "energy": None, "xyz": None, "ok": False, "seconds": 0}
    elapsed = time.time() - t0

    if rc != 0:
        logger.warning(f"CONF_{idx}: xtb 返回非零状态 {rc}；日志末尾：\n{tail_of(sp_log)}")

    energy, source = extract_energy_detail(sp_log)
    if energy is None:
        logger.error(
            f"CONF_{idx}: 未能从 stdout 解析总能量（应含 '| TOTAL ENERGY ... Eh |'）；"
            f"日志末尾：\n{tail_of(sp_log)}"
        )
        return {"index": idx, "energy": None, "xyz": None, "ok": False,
                "seconds": elapsed}
    if source != "TOTAL ENERGY box":
        logger.warning(
            f"CONF_{idx}: 未找到 TOTAL ENERGY 盒线，改由「{source}」取得能量 "
            f"{energy:.12f} Eh；请核对 {sp_log}"
        )

    if not check_termination(sp_log):
        logger.warning(f"CONF_{idx}: 未见正常收尾标记，请检查 {sp_log}")

    tagged_xyz = conf_dir / f"CONF_{idx}_gxtb.xyz"
    if not update_xyz_energy(work_xyz, tagged_xyz, energy, solvent_tag):
        logger.error(f"CONF_{idx}: 写入能量注释失败")
        return {"index": idx, "energy": None, "xyz": None, "ok": False,
                "seconds": elapsed}

    logger.info(f"CONF_{idx}: E = {energy:.12f} Eh  ({elapsed:.1f}s)")
    return {"index": idx, "energy": energy, "xyz": tagged_xyz, "ok": True,
            "seconds": elapsed}


# --------------------------------------------------------------------------- #
# Step 1：SDF 拆分
# --------------------------------------------------------------------------- #
def split_sdf(input_sdf, logger):
    out_dir = Path("1_conf")
    # 清掉上次拆分残留，否则 CONF_ 数量会累加导致误判构象数
    for old in out_dir.glob("CONF_*.xyz"):
        try:
            old.unlink()
        except Exception:
            pass

    logger.info(f"拆分 SDF: {input_sdf} -> {out_dir}/CONF_*.xyz")
    try:
        proc = subprocess.run(
            ["obabel", "-isdf", str(input_sdf), "-oxyz",
             "-O", str(out_dir / "CONF_.xyz"), "-m"],
            capture_output=True, text=True
        )
    except Exception as e:
        logger.error(f"调用 obabel 失败 - {e}")
        sys.exit(1)

    if proc.returncode != 0:
        logger.error(f"obabel 返回 {proc.returncode}: {proc.stderr.strip()}")
        sys.exit(1)

    conf_files = sorted(out_dir.glob("CONF_*.xyz"))
    if not conf_files:
        logger.error(f"{out_dir} 下没有生成任何构象，检查输入 SDF 内容")
        sys.exit(1)

    logger.info(f"共拆出 {len(conf_files)} 个构象")
    return len(conf_files)


# --------------------------------------------------------------------------- #
# Step 3：汇总系综 + 能量表
# --------------------------------------------------------------------------- #
def build_ensemble(results, logger):
    """把各构象（已带能量注释）拼成 ensemble_gxtb.xyz，并写能量汇总 CSV。"""
    ensemble = Path("2_conf_sp/ensemble_gxtb.xyz")
    ok = [r for r in results if r["ok"]]
    ok.sort(key=lambda r: r["index"])

    if not ok:
        logger.error("没有任何构象成功算出能量，无法构建系综")
        sys.exit(1)

    # 相对能量以本次成功子集的最低者为参考（与 CREST cregen 处理一致）
    emin = min(r["energy"] for r in ok)
    ref_idx = next(r["index"] for r in ok if r["energy"] == emin)

    with open(ensemble, "w") as dst:
        for r in ok:
            dst.write(Path(r["xyz"]).read_text())
    logger.info(f"系综文件已写入 {ensemble}（{len(ok)} 帧）")

    # 能量汇总表，便于核对；ΔE 单位为 kcal/mol
    csv_path = Path("2_conf_sp/energy_summary.csv")
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("conformer,E_Eh,dE_kcal,seconds\n")
        for r in ok:
            f.write(f"CONF_{r['index']},{r['energy']:.12f},"
                    f"{(r['energy'] - emin) * HARTREE2KCAL:.3f},{r['seconds']:.1f}\n")
        for r in results:
            if not r["ok"]:
                f.write(f"CONF_{r['index']},FAILED,,\n")
    logger.info(f"能量汇总表已写入 {csv_path}")
    logger.info(f"参考构象（最低能）: CONF_{ref_idx}，E = {emin:.12f} Eh")
    logger.info(f"能量跨度 ΔE = "
                f"{(max(r['energy'] for r in ok) - emin) * HARTREE2KCAL:.2f} kcal/mol")

    return ok, emin, ref_idx


# --------------------------------------------------------------------------- #
# Step 4：CREST cregen 聚类
# --------------------------------------------------------------------------- #
def run_cregen(ok_results, emin, ref_idx, args, logger):
    work = Path("3_cluster")
    work.mkdir(parents=True, exist_ok=True)
    for old in work.glob("*"):
        if old.is_file():
            try:
                old.unlink()
            except Exception:
                pass

    # 两个路径必须区分开：
    #   ensemble_local —— 相对脚本当前 cwd（项目根目录），仅用于本进程检查文件是否存在
    #   ensemble_arg   —— 相对 work(3_cluster)，因为 crest 以 cwd=3_cluster 启动
    # 混用会导致存在性检查指向项目根目录的上一级，误报"系综文件缺失"。
    ensemble_local = Path("2_conf_sp") / "ensemble_gxtb.xyz"
    ensemble_arg = Path("..") / "2_conf_sp" / "ensemble_gxtb.xyz"
    if not ensemble_local.is_file():
        logger.error(f"系综文件缺失: {ensemble_local.resolve()}")
        sys.exit(1)

    # cregen 需要一份参考坐标（只用于确定原子与电荷），取最低能构象最合理
    ref_src = next(r["xyz"] for r in ok_results if r["index"] == ref_idx)
    ref_name = ref_src.name
    try:
        shutil.copy2(ref_src, work / ref_name)
    except Exception as e:
        logger.error(f"复制参考构象失败 - {e}")
        sys.exit(1)

    cmd = [
        "crest", ref_name,
        "--cregen", str(ensemble_arg),
        "-chrg", str(args.charge),
        "-uhf", str(args.uhf),
        "-rthr", str(args.rthr),
        "-ewin", str(args.ewin),
    ]
    logger.info(f"CREST 聚类: {' '.join(cmd)}")
    logger.info(f"参数: RMSD 阈值 {args.rthr} Å, 能量窗口 {args.ewin} kcal/mol, "
                f"参考构象 CONF_{ref_idx}")

    try:
        proc = subprocess.run(
            cmd, cwd=str(work), capture_output=True, text=True
        )
    except Exception as e:
        logger.error(f"调用 crest 失败 - {e}")
        sys.exit(1)

    log_path = work / "crest.log"
    with open(log_path, "w", encoding="utf-8") as f:
        f.write(proc.stdout or "")
        if proc.stderr:
            f.write("\n\n=== STDERR ===\n")
            f.write(proc.stderr)

    if proc.returncode != 0:
        logger.warning(f"crest 返回 {proc.returncode}，详见 {log_path}")
    else:
        logger.info("CREST 聚类完成")

    # cregen 的产出文件名随版本略有差异，都检查一遍
    candidates = ["crest_ensemble.xyz", "crest_conformers.xyz", "crest_best.xyz"]
    found = [c for c in candidates if (work / c).is_file()]
    if found:
        for c in found:
            logger.info(f"聚类结果: {work / c}")
    else:
        logger.warning(f"未发现聚类结果文件（{', '.join(candidates)}），请检查 {log_path}")
        return None

    return work / found[0]


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    parser = argparse.ArgumentParser(
        description="构象单点能 + CREST 聚类工作流（g-xTB 经 xtb --gxtb 调用）",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python gxtb_cluster.py mol.sdf --charge 0
  python gxtb_cluster.py mol.sdf --charge 0 --jobs 4
  python gxtb_cluster.py mol.sdf --charge -1 --solvent none --ewin 6

默认对每个构象执行:
  xtb CONF_i.xyz -c <charge> -u 0 --sp --gxtb --cosmo water
能量取自 stdout 日志中的 '| TOTAL ENERGY ... Eh |' 行（已含水相溶剂化）
"""
    )
    parser.add_argument("input_sdf", type=str, help="输入 SDF 文件路径")
    parser.add_argument("--charge", type=int, required=True, help="分子的形式电荷")
    parser.add_argument("--uhf", type=int, default=0,
                        help="未成对电子数，基态闭壳层用 0（默认 0）")
    parser.add_argument("--solvent", type=str, default="water",
                        help="ddCOSMO 溶剂名，默认 water；填 none 表示气相")
    parser.add_argument("--rthr", type=float, default=0.125,
                        help="CREST 聚类 RMSD 阈值 (Å)，默认 0.125")
    parser.add_argument("--ewin", type=float, default=10.0,
                        help="CREST 能量窗口 (kcal/mol)，默认 10")
    parser.add_argument("--jobs", type=int, default=1,
                        help="并行处理的构象数，默认 1（xtb 自身已多线程）")
    parser.add_argument("--xtb-threads", type=int, default=None,
                        help="传给 xtb 的 -P 线程数，默认不传（用全部核）")
    args = parser.parse_args()

    logger = setup_logger()

    # 溶剂名归一化：none/gas 等写法统一视为不加溶剂
    if args.solvent.lower() in ("none", "no", "off", "gas", ""):
        args.solvent = None
    solvent_tag = args.solvent if args.solvent else "gas"

    input_sdf = Path(args.input_sdf)
    if not input_sdf.is_file():
        logger.error(f"输入文件不存在: {input_sdf}")
        sys.exit(1)

    logger.info("=" * 70)
    logger.info(f"输入: {input_sdf}")
    logger.info(f"电荷 {args.charge}, 未成对电子 {args.uhf}, 溶剂 {solvent_tag}")
    logger.info("方法: xtb --gxtb  |  任务: --sp 单点能  |  聚类: CREST cregen")
    logger.info("=" * 70)

    # --- 环境检查 ---
    ok_env = True
    for tool in ("xtb", "obabel", "crest"):
        if probe_tool(tool, logger) is None:
            ok_env = False
    if not ok_env:
        logger.error("缺少必要可执行文件，退出")
        sys.exit(1)
    if not check_gxtb_support(logger):
        sys.exit(1)

    # --- 目录 ---
    for d in ("1_conf", "2_conf_sp", "3_cluster"):
        Path(d).mkdir(exist_ok=True)

    # --- Step 1: SDF 拆分 ---
    logger.info("-" * 70)
    logger.info("Step 1: SDF 拆分为单构象 XYZ")
    num_conformers = split_sdf(input_sdf, logger)
    args.total = num_conformers

    # 用第一个构象校验电荷/自旋自洽性（错配时 xtb 会直接失败）
    if not validate_uhf(Path("1_conf") / "CONF_1.xyz", args.charge, args.uhf, logger):
        logger.error("请检查 --charge 与 --uhf 设置")
        sys.exit(1)

    # --- Step 2: g-xTB 单点能 ---
    logger.info("-" * 70)
    logger.info(f"Step 2: g-xTB 单点能（{num_conformers} 个构象，并行 {args.jobs}）")
    logger.info(f"  命令模板: xtb <xyz> -c {args.charge} -u {args.uhf} --sp --gxtb"
                + (f" --cosmo {args.solvent}" if args.solvent else ""))

    t_start = time.time()
    results = []

    if args.jobs <= 1:
        for i in range(1, num_conformers + 1):
            results.append(process_conformer(i, args, logger, solvent_tag))
    else:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {
                pool.submit(process_conformer, i, args, logger, solvent_tag): i
                for i in range(1, num_conformers + 1)
            }
            for fut in as_completed(futures):
                i = futures[fut]
                try:
                    results.append(fut.result())
                except Exception as e:
                    logger.error(f"CONF_{i}: 处理异常 - {e}")
                    results.append({"index": i, "energy": None, "xyz": None,
                                    "ok": False, "seconds": 0})

    results.sort(key=lambda r: r["index"])
    n_ok = sum(1 for r in results if r["ok"])
    logger.info(f"单点能完成: {n_ok}/{num_conformers} 成功，"
                f"总耗时 {time.time() - t_start:.1f}s")
    if n_ok == 0:
        logger.error("全部构象失败，退出。请查看各 2_conf_sp/CONF_i/*.log")
        sys.exit(1)
    if n_ok < num_conformers:
        failed = [f"CONF_{r['index']}" for r in results if not r["ok"]]
        logger.warning(f"失败的构象: {', '.join(failed)}（已从系综中剔除）")

    # --- Step 3: 汇总系综 ---
    logger.info("-" * 70)
    logger.info("Step 3: 汇总系综与能量表")
    ok_results, emin, ref_idx = build_ensemble(results, logger)

    # --- Step 4: CREST 聚类 ---
    logger.info("-" * 70)
    logger.info("Step 4: CREST cregen 去重聚类")
    final = run_cregen(ok_results, emin, ref_idx, args, logger)

    # --- 收尾 ---
    logger.info("-" * 70)
    logger.info("完成。产出：")
    logger.info("  1_conf/CONF_*.xyz                拆分出的原始构象")
    logger.info("  2_conf_sp/CONF_i/CONF_i_gxtb.log 每次 xtb 运行的完整输出")
    logger.info("  2_conf_sp/CONF_i/CONF_i_gxtb.xyz 带能量注释的单构象")
    logger.info("  2_conf_sp/ensemble_gxtb.xyz      全部构象系综")
    logger.info("  2_conf_sp/energy_summary.csv     能量汇总（Eh 与 ΔE kcal/mol）")
    if final:
        logger.info(f"  {final}                         聚类结果系综")


if __name__ == "__main__":
    main()
