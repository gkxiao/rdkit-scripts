#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chembl37_adme_search.py
=======================
RDKit 子结构 / 相似性 / 骨架搜索（ChEMBL CSV），按 standard_type 分组做单位归一化。

单位归一化原则
--------------
standard_type + standard_value + standard_units 是三元组，只有在同一
standard_type 内部换算才有物理意义：
  * 内置默认目标单位（可用 --normalize-to TYPE:UNIT 覆盖）
        CL / CLint / CLrenal / CL/F  -> mL.min-1.kg-1
        Vd / Vdss / Vss              -> L.kg-1
        T1/2 / MRT / Tmax            -> h
        AUC / AUClast / AUCinf       -> h.ng.mL-1
        Cmax                         -> ng.mL-1
        F                            -> %
  * 未在映射表中的 standard_type 不做换算，计入 norm_failed
  * 数值阈值必须绑定 standard_type：`--max-value CL:10`，或先用
    `--value-type CL` 限定。若数据里存在多个 standard_type 而只给了裸数字，
    程序直接报错并列出所有类型，避免把 CL 阈值套到 Vd 上。

示例
----
# 1) 只筛 CL <= 10 mL/min/kg（自动按 CL 的目标单位换算）
python chembl37_adme_search.py -i chembl37_adme.csv -mode similarity \
    -query 'COc1ccc(NS(=O)(=O)C(F)(F)F)cc1' -o hits.csv -cutoff 0.7 \
    --max-hits 100 --normalize-units --max-value CL:10

# 2) 自定义某类型的归一化目标
python chembl37_adme_search.py -i chembl37_adme.csv -mode substruct \
    -query 'c1ccc2[nH]ccc2c1' -o hits.csv --normalize-units \
    --normalize-to 'CL:L.h-1.kg-1' --normalize-to 'Vd:L.kg-1'

# 3) 一个全局目标单位（仅建议在只含单一 standard_type 时使用）
python chembl37_adme_search.py -i chembl37_adme.csv -mode similarity \
    -query 'CCO' -o hits.csv --normalize-units --normalize-to 'mL.min-1.kg-1'
"""

from __future__ import annotations

import argparse
import csv
import heapq
import math
import os
import sys
import time
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

try:
    from rdkit.Chem import rdFingerprintGenerator as _fpg
    _HAS_FPG = True
except Exception:
    _HAS_FPG = False

RDLogger.DisableLog("rdApp.*")

SEARCH_EXTRA_FIELDS = ("tanimoto", "match_atoms", "murcko_scaffold")
NORM_EXTRA_FIELDS = ("standard_value_norm", "standard_units_norm")

# standard_type -> 默认目标单位（键统一小写）
DEFAULT_NORM_TARGETS = {
    "cl": "mL.min-1.kg-1",
    "clint": "mL.min-1.kg-1",
    "clint,u": "mL.min-1.kg-1",
    "clrenal": "mL.min-1.kg-1",
    "cl/f": "mL.min-1.kg-1",
    "vd": "L.kg-1",
    "vdss": "L.kg-1",
    "vss": "L.kg-1",
    "t1/2": "h",
    "mrt": "h",
    "tmax": "h",
    "auc": "h.ng.mL-1",
    "auclast": "h.ng.mL-1",
    "aucinf": "h.ng.mL-1",
    "cmax": "ng.mL-1",
    "f": "%",
}

_G: Dict[str, object] = {}


# --------------------------------------------------------------------------- #
# 单位解析 / 换算
# --------------------------------------------------------------------------- #
# 单位 -> (量纲名, 在基准单位下的大小)。基准：L，min，kg，"(1 = 100%)"
_UNIT_SIZE = {
    "nl": ("vol", 1e-9), "ul": ("vol", 1e-6), "ml": ("vol", 1e-3),
    "l": ("vol", 1.0),
    "s": ("time", 1.0 / 60.0), "min": ("time", 1.0), "h": ("time", 60.0),
    "hr": ("time", 60.0), "hrs": ("time", 60.0), "d": ("time", 1440.0),
    "pg": ("mass", 1e-15), "ng": ("mass", 1e-12), "ug": ("mass", 1e-9),
    "mg": ("mass", 1e-6), "g": ("mass", 1e-3), "kg": ("mass", 1.0),
    "%": ("pct", 1.0), "percent": ("pct", 1.0),
}
_DIM_ORDER = ("vol", "time", "mass", "pct")


def parse_unit(unit) -> Optional[Tuple[float, Tuple[int, int, int, int]]]:
    """'mL.min-1.g-1' -> (1.0, (1, -1, -1, 0))；返回 (倍率, 量纲元组)。"""
    if unit is None:
        return None
    s = str(unit).strip().lower()
    for a, b in (("µ", "u"), ("μ", "u"), ("·", "."), ("⋅", "."), ("*", ".")):
        s = s.replace(a, b)
    s = s.replace("⁻¹", "-1").replace("^", "").replace(" ", "")
    if not s:
        return None

    parts = s.split("/")
    tokens: List[Tuple[str, int]] = []
    for i, part in enumerate(parts):
        sign = 1 if i == 0 else -1
        for tok in part.split("."):
            if tok:
                tokens.append((tok, sign))

    factor = 1.0
    dims = {k: 0 for k in _DIM_ORDER}
    for tok, sign in tokens:
        exp = sign
        if tok.endswith("-1"):
            tok, exp = tok[:-2], -sign
        if tok not in _UNIT_SIZE:
            return None
        kind, size = _UNIT_SIZE[tok]
        dims[kind] += exp
        factor *= size ** exp
    return factor, tuple(dims[k] for k in _DIM_ORDER)


def to_float(text) -> Optional[float]:
    try:
        v = float(str(text).strip())
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def normalize_row_value(row, cfg):
    """返回 (归一化值, 目标单位字符串, 错误标记)。成功时错误标记为 None。"""
    v = to_float(row.get(cfg["value_col"]))
    if v is None:
        return None, None, "no_value"

    typ = (row.get(cfg["type_col"]) or "").strip().lower()
    entry = cfg["norm_map"].get(typ)
    if entry is None:
        return None, None, "unmapped_type"

    target_unit, target_factor, target_dims = entry
    parsed = parse_unit(row.get(cfg["units_col"]))
    if parsed is None:
        return None, None, "unknown_unit"
    factor, dims = parsed
    if dims != tuple(target_dims):
        return None, None, "dim_mismatch"
    return v * factor / target_factor, target_unit, None


# --------------------------------------------------------------------------- #
# 预扫描 standard_type / standard_units
# --------------------------------------------------------------------------- #
def prescan(args):
    types_count: Dict[str, int] = {}
    units_count: Dict[str, Dict[str, int]] = {}
    n = 0
    if not os.path.exists(args.input):
        raise SystemExit("输入文件不存在: %s" % args.input)
    with open(args.input, "r", newline="", encoding=args.encoding) as fh:
        reader = csv.DictReader(fh, delimiter=args.delimiter)
        if not reader.fieldnames:
            raise SystemExit("CSV 没有表头: %s" % args.input)
        for i, row in enumerate(reader, start=1):
            if args.limit and i > args.limit:
                break
            n += 1
            t = (row.get(args.type_col) or "").strip()
            u = (row.get(args.units_col) or "").strip()
            types_count[t] = types_count.get(t, 0) + 1
            bucket = units_count.setdefault(t, {})
            bucket[u] = bucket.get(u, 0) + 1
    return types_count, units_count, n


# --------------------------------------------------------------------------- #
# 归一化映射 / 阈值
# --------------------------------------------------------------------------- #
def parse_normalize_to(items):
    overrides, global_spec = {}, None
    for s in items or []:
        s = s.strip()
        if ":" in s:
            typ, unit = s.split(":", 1)
            overrides[typ.strip().lower()] = unit.strip()
        else:
            if global_spec is not None:
                raise SystemExit("只能给一个不带 TYPE: 前缀的全局 --normalize-to")
            global_spec = s
    return overrides, global_spec


def build_norm_map(types_found, overrides, global_spec):
    norm_map: Dict[str, tuple] = {}
    plan: Dict[str, Optional[str]] = {}
    for typ in types_found:
        key = typ.strip().lower()
        spec = overrides.get(key)
        if spec is None:
            spec = DEFAULT_NORM_TARGETS.get(key)
        if spec is None and global_spec is not None:
            spec = global_spec
        if spec is None:
            plan[typ] = None
            continue
        parsed = parse_unit(spec)
        if parsed is None:
            raise SystemExit("无法解析目标单位 %r（standard_type=%r）" % (spec, typ))
        if key in norm_map and norm_map[key][0] != spec:
            raise SystemExit("standard_type=%r 的目标单位冲突" % typ)
        norm_map[key] = (spec, parsed[0], parsed[1])
        plan[typ] = spec
    return norm_map, plan


def build_thresholds(args, types_found):
    raw: Dict[str, List[Optional[float]]] = {}
    for token, is_min in ((args.min_value, True), (args.max_value, False)):
        if token is None:
            continue
        token = str(token).strip()
        if ":" in token:
            typ, num = token.split(":", 1)
            key = typ.strip().lower()
        else:
            key, num = "*", token
        try:
            v = float(num)
        except ValueError:
            raise SystemExit("数值阈值格式错误: %r（应为 10 或 CL:10）" % token)
        slot = raw.setdefault(key, [None, None])
        slot[0 if is_min else 1] = v

    if not raw:
        return {}

    if "*" in raw:
        glob = raw.pop("*")
        if len(types_found) == 1:
            raw[next(iter(types_found)).strip().lower()] = list(glob)
        elif args.value_type:
            for t in args.value_type:
                raw.setdefault(t.strip().lower(), list(glob))
        else:
            listing = ", ".join(repr(t) for t in types_found)
            raise SystemExit(
                "检测到多个 standard_type (%s)，单一数值阈值无法跨类型比较。\n"
                "  请改用：--max-value CL:10  或  --value-type CL --max-value 10\n"
                "  或先用 -filter standard_type=CL 限定数据。" % listing)

    known = {t.strip().lower() for t in types_found}
    for k in raw:
        if k not in known:
            print("[warn] 阈值指向的 standard_type=%r 在数据中不存在" % k,
                  file=sys.stderr)
    return raw


# --------------------------------------------------------------------------- #
# 指纹 / 骨架 / 查询
# --------------------------------------------------------------------------- #
class _LegacyFPGen:
    def __init__(self, kind="morgan", radius=2, nbits=2048, chirality=False):
        self.kind, self.radius = kind, radius
        self.nbits, self.chirality = nbits, chirality

    def GetFingerprint(self, mol):
        if self.kind == "morgan":
            return AllChem.GetMorganFingerprintAsBitVect(
                mol, self.radius, nBits=self.nbits, useChirality=self.chirality)
        if self.kind == "rdkit":
            return AllChem.RDKFingerprint(mol, fpSize=self.nbits)
        return AllChem.GetMACCSKeysFingerprint(mol)


def make_fpgen(cfg):
    kind = cfg["fp_type"]
    if _HAS_FPG:
        if kind == "morgan":
            return _fpg.GetMorganGenerator(
                radius=cfg["radius"], fpSize=cfg["fp_bits"],
                includeChirality=cfg["chirality"])
        if kind == "rdkit":
            return _fpg.GetRDKitFPGenerator(fpSize=cfg["fp_bits"])
        return _fpg.GetMACCSGenerator()
    return _LegacyFPGen(kind, radius=cfg["radius"],
                        nbits=cfg["fp_bits"], chirality=cfg["chirality"])


def scaffold_of(mol, cfg):
    try:
        scaf = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception:
        return None
    if scaf is None:
        return None
    if cfg["generic_scaffold"] and scaf.GetNumAtoms() > 0:
        try:
            scaf = MurckoScaffold.MakeScaffoldGeneric(scaf)
        except Exception:
            pass
    return scaf


def prepare_query(cfg) -> Dict[str, object]:
    mode, q = cfg["mode"], cfg["query"]
    g: Dict[str, object] = {"fpgen": make_fpgen(cfg)}

    if mode == "substruct":
        qmol = None
        if cfg["query_format"] in ("auto", "smarts"):
            qmol = Chem.MolFromSmarts(q)
        if qmol is None and cfg["query_format"] in ("auto", "smiles"):
            qmol = Chem.MolFromSmiles(q)
        if qmol is None:
            raise ValueError("无法解析 query（SMARTS/SMILES 均失败）: %r" % q)
        g["qmol"] = qmol
        return g

    qmol = Chem.MolFromSmiles(q)
    if qmol is None:
        raise ValueError("无法解析 query SMILES: %r" % q)

    if mode == "similarity":
        g["qfp"] = g["fpgen"].GetFingerprint(qmol)
        return g

    if cfg["query_is_scaffold"]:
        qscaf = qmol
        if cfg["generic_scaffold"]:
            try:
                qscaf = MurckoScaffold.MakeScaffoldGeneric(qscaf)
            except Exception:
                pass
    else:
        qscaf = scaffold_of(qmol, cfg)
    if qscaf is None or qscaf.GetNumAtoms() == 0:
        raise ValueError("query 没有环系 / 无法生成 Murcko 骨架: %r" % q)

    g["qscaf"] = qscaf
    g["qscaf_smiles"] = Chem.MolToSmiles(qscaf)
    if cfg["scaffold_match"] == "similar":
        g["qscaf_fp"] = g["fpgen"].GetFingerprint(qscaf)
    return g


def describe_query(g, cfg) -> str:
    if cfg["mode"] == "substruct":
        return "SMARTS = " + Chem.MolToSmarts(g["qmol"])
    if cfg["mode"] == "similarity":
        return "SMILES = %s | %s radius=%d bits=%d cutoff=%.3f" % (
            cfg["query"], cfg["fp_type"], cfg["radius"],
            cfg["fp_bits"], cfg["cutoff"])
    return "scaffold SMILES = %s | match=%s cutoff=%.3f" % (
        g["qscaf_smiles"], cfg["scaffold_match"], cfg["cutoff"])


# --------------------------------------------------------------------------- #
# 过滤
# --------------------------------------------------------------------------- #
def col_filters_pass(row, cfg) -> bool:
    for col, want in cfg["filters"]:
        if (row.get(col) or "").strip() != want:
            return False
    return True


# --------------------------------------------------------------------------- #
# 单条记录处理
# --------------------------------------------------------------------------- #
def _init_worker(cfg):
    global _G
    _G = {"cfg": cfg}
    _G.update(prepare_query(cfg))


def process_chunk(rows):
    cfg = _G["cfg"]
    mode, fpgen = cfg["mode"], _G["fpgen"]
    norm_on = cfg["normalize_units"]
    has_thr = cfg["has_thresholds"]

    hits: List[Tuple[float, int, dict, dict]] = []
    n_total = n_filtered = n_invalid = n_match = n_fail = n_other = 0

    for idx, row in rows:
        n_total += 1
        if not col_filters_pass(row, cfg):
            n_filtered += 1
            continue

        # ---- 数值：按 standard_type 归一化 ----
        base_extra: Dict[str, object] = {}
        val: Optional[float] = None
        if norm_on:
            val, tgt_unit, err = normalize_row_value(row, cfg)
            if err is not None:
                n_fail += 1
            else:
                base_extra["standard_value_norm"] = round(val, 6)
                base_extra["standard_units_norm"] = tgt_unit
        else:
            val = to_float(row.get(cfg["value_col"]))

        # ---- 数值阈值（必须绑定 standard_type）----
        if has_thr:
            typ = (row.get(cfg["type_col"]) or "").strip().lower()
            thr = cfg["thresholds"].get(typ)
            if thr is None:
                n_other += 1          # 该类型没有阈值 -> 不参与筛选
                continue
            if val is None:
                n_filtered += 1
                continue
            if (thr[0] is not None and val < thr[0]) or \
               (thr[1] is not None and val > thr[1]):
                n_filtered += 1
                continue

        # ---- 结构匹配 ----
        smi = (row.get(cfg["smiles_col"]) or "").strip()
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            n_invalid += 1
            continue

        if mode == "substruct":
            match = mol.GetSubstructMatch(_G["qmol"])
            if not match:
                continue
            n_match += 1
            extra = dict(base_extra)
            extra["match_atoms"] = len(match)
            hits.append((0.0, idx, row, extra))

        elif mode == "similarity":
            sim = DataStructs.TanimotoSimilarity(
                _G["qfp"], fpgen.GetFingerprint(mol))
            if sim + 1e-12 < cfg["cutoff"]:
                continue
            n_match += 1
            sim = round(sim, 4)
            extra = dict(base_extra)
            extra["tanimoto"] = sim
            hits.append((sim, idx, row, extra))

        else:  # scaffold
            scaf = scaffold_of(mol, cfg)
            if scaf is None or scaf.GetNumAtoms() == 0:
                n_invalid += 1
                continue
            scaf_smi = Chem.MolToSmiles(scaf)
            how = cfg["scaffold_match"]
            extra = dict(base_extra)

            if how == "exact":
                if scaf_smi != _G["qscaf_smiles"]:
                    continue
                score = 1.0
                extra["murcko_scaffold"] = scaf_smi
            elif how == "similar":
                score = DataStructs.TanimotoSimilarity(
                    _G["qscaf_fp"], fpgen.GetFingerprint(scaf))
                if score + 1e-12 < cfg["cutoff"]:
                    continue
                score = round(score, 4)
                extra["murcko_scaffold"] = scaf_smi
                extra["tanimoto"] = score
            else:
                if not scaf.HasSubstructMatch(_G["qscaf"]):
                    continue
                score = 0.0
                extra["murcko_scaffold"] = scaf_smi

            n_match += 1
            hits.append((score, idx, row, extra))

    if cfg["max_hits"] > 0 and len(hits) > cfg["max_hits"]:
        hits = heapq.nlargest(cfg["max_hits"], hits, key=lambda h: (h[0], -h[1]))

    stats = {"total": n_total, "filtered": n_filtered, "invalid": n_invalid,
             "matched": n_match, "norm_failed": n_fail, "other_type": n_other}
    return stats, hits


# --------------------------------------------------------------------------- #
# 迭代工具 / 主搜索
# --------------------------------------------------------------------------- #
def iter_rows(reader, limit=0) -> Iterator[Tuple[int, dict]]:
    for i, row in enumerate(reader, start=1):
        if limit and i > limit:
            return
        yield i, row


def chunked(it, size: int):
    buf = []
    for x in it:
        buf.append(x)
        if len(buf) >= size:
            yield buf
            buf = []
    if buf:
        yield buf


def parse_filters(items: Sequence[str]) -> List[Tuple[str, str]]:
    out = []
    for s in items or []:
        if "=" not in s:
            raise SystemExit("--filter 需要 COL=VAL 格式，收到: %r" % s)
        col, val = s.split("=", 1)
        out.append((col.strip(), val.strip()))
    return out


def run_search(cfg, args):
    stats = {"total": 0, "filtered": 0, "invalid": 0,
             "matched": 0, "norm_failed": 0, "other_type": 0}
    heap_mode = cfg["max_hits"] > 0
    heap: List[tuple] = []
    all_hits: List[tuple] = []
    t0 = time.time()

    with open(args.input, "r", newline="", encoding=args.encoding) as fh:
        reader = csv.DictReader(fh, delimiter=args.delimiter)
        header = list(reader.fieldnames)
        chunks = chunked(iter_rows(reader, args.limit), max(1, args.chunk_size))
        next_mark = args.progress_every

        def consume(result):
            nonlocal next_mark
            chunk_stats, chunk_hits = result
            for k in stats:
                stats[k] += chunk_stats[k]
            if heap_mode:
                for score, idx, row, extra in chunk_hits:
                    item = (score, -idx, idx, row, extra)
                    if len(heap) < cfg["max_hits"]:
                        heapq.heappush(heap, item)
                    elif item[:2] > heap[0][:2]:
                        heapq.heapreplace(heap, item)
            else:
                all_hits.extend(chunk_hits)
            if not args.quiet and stats["total"] >= next_mark:
                print("\r[scan] %d rows | %d matched | %.0fs"
                      % (stats["total"], stats["matched"], time.time() - t0),
                      end="", file=sys.stderr, flush=True)
                while stats["total"] >= next_mark:
                    next_mark += args.progress_every

        if args.jobs and args.jobs > 1:
            from concurrent.futures import ProcessPoolExecutor
            with ProcessPoolExecutor(max_workers=args.jobs,
                                     initializer=_init_worker,
                                     initargs=(cfg,)) as ex:
                for res in ex.map(process_chunk, chunks, chunksize=1):
                    consume(res)
        else:
            _init_worker(cfg)
            for res in map(process_chunk, chunks):
                consume(res)

    if not args.quiet:
        print("", file=sys.stderr)
    hits = [(s, i, r, e) for s, _n, i, r, e in heap] if heap_mode else all_hits
    return stats, hits, header


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser():
    p = argparse.ArgumentParser(
        prog="chembl37_adme_search.py",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="RDKit 子结构 / 相似性 / 骨架搜索 + 按 standard_type 归一化",
        epilog=__doc__)

    g = p.add_argument_group("输入 / 输出")
    g.add_argument("-i", "--input", default="/public/gkxiao/software/chembl/chembl_37_adme_stability.csv")
    g.add_argument("-o", "--output", default="hits.csv")
    g.add_argument("--smiles-col", default="canonical_smiles")
    g.add_argument("--id-col", default="compound_chembl_id")
    g.add_argument("--type-col", default="standard_type")
    g.add_argument("--delimiter", default=",")
    g.add_argument("--encoding", default="utf-8-sig")
    g.add_argument("--limit", type=int, default=0)

    g = p.add_argument_group("搜索")
    g.add_argument("-mode", "--mode", required=True,
                   choices=["substruct", "similarity", "scaffold"])
    g.add_argument("-query", "--query", required=True)
    g.add_argument("--query-format", default="auto",
                   choices=["auto", "smiles", "smarts"])
    g.add_argument("-cutoff", "--cutoff", type=float, default=0.7)
    g.add_argument("-max-hits", "--max-hits", type=int, default=None)
    g.add_argument("--sort", default="auto", choices=["auto", "input", "score"])
    g.add_argument("--dedupe", action="store_true")

    g = p.add_argument_group("指纹")
    g.add_argument("--fp-type", default="morgan", choices=["morgan", "rdkit", "maccs"])
    g.add_argument("--radius", type=int, default=2)
    g.add_argument("--fp-bits", type=int, default=2048)
    g.add_argument("--chirality", action="store_true")

    g = p.add_argument_group("骨架选项")
    g.add_argument("--scaffold-match", default="substruct",
                   choices=["substruct", "exact", "similar"])
    g.add_argument("--generic-scaffold", action="store_true")
    g.add_argument("--query-is-scaffold", action="store_true")

    g = p.add_argument_group("数值过滤 / 单位归一化（按 standard_type）")
    g.add_argument("-filter", "--filter", action="append", metavar="COL=VAL",
                   default=[])
    g.add_argument("--value-col", default="standard_value")
    g.add_argument("--units-col", default="standard_units")
    g.add_argument("--min-value", default=None,
                   help="10 或 CL:10")
    g.add_argument("--max-value", default=None,
                   help="10 或 CL:10")
    g.add_argument("--value-type", action="append", default=[],
                   help="裸数字阈值作用到的 standard_type，可重复")
    g.add_argument("--normalize-units", action="store_true",
                   help="按 standard_type 把 standard_value 换算到各自目标单位")
    g.add_argument("--normalize-to", action="append", default=None,
                   metavar="[TYPE:]UNIT",
                   help="覆盖目标单位，如 CL:L.h-1.kg-1；不带 TYPE: 则为全局默认")

    g = p.add_argument_group("性能")
    g.add_argument("-j", "-jobs", "--jobs", dest="jobs", type=int, default=1)
    g.add_argument("--chunk-size", type=int, default=2000)
    g.add_argument("--progress-every", type=int, default=20000)
    g.add_argument("-q", "--quiet", action="store_true")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    if not 0.0 <= args.cutoff <= 1.0:
        raise SystemExit("--cutoff 必须在 [0,1] 区间")

    max_hits = args.max_hits
    if max_hits is None:
        max_hits = 0 if args.mode == "substruct" else 100
    if max_hits < 0:
        max_hits = 0

    # ---- 预扫描：standard_type x standard_units ----
    types_count, units_count, n_prescanned = prescan(args)
    if not args.quiet:
        print("[types] %d standard_type(s), %d 行" % (len(types_count), n_prescanned),
              file=sys.stderr)

    # ---- 归一化映射 ----
    norm_on = bool(args.normalize_units)
    norm_map: Dict[str, tuple] = {}
    plan: Dict[str, Optional[str]] = {}
    if norm_on:
        overrides, global_spec = parse_normalize_to(args.normalize_to)
        norm_map, plan = build_norm_map(types_count.keys(), overrides, global_spec)
        if not norm_map:
            raise SystemExit(
                "没有任何 standard_type 能确定目标单位；请用 "
                "--normalize-to TYPE:UNIT 指定，或检查 --type-col。")
        if not args.quiet:
            for typ, cnt in sorted(types_count.items(), key=lambda kv: -kv[1]):
                top = ", ".join(
                    "%s:%d" % (u or "<empty>", c)
                    for u, c in sorted(units_count[typ].items(),
                                       key=lambda kv: -kv[1])[:4])
                print("  %-12s %7d rows | %s -> %s"
                      % (typ or "<empty>", cnt, top,
                         plan.get(typ) or "不归一化"), file=sys.stderr)

    # ---- 数值阈值（绑定 standard_type）----
    thresholds = build_thresholds(args, types_count)
    if thresholds and not norm_on:
        multi_unit = any(len(v) > 1 for v in units_count.values())
        if multi_unit:
            print("[warn] 未启用 --normalize-units，但同一 standard_type 下存在多种 "
                  "standard_units，阈值比较可能失真。", file=sys.stderr)

    cfg = {
        "mode": args.mode,
        "query": args.query,
        "query_format": args.query_format,
        "cutoff": args.cutoff,
        "max_hits": max_hits,
        "smiles_col": args.smiles_col,
        "id_col": args.id_col,
        "type_col": args.type_col,
        "fp_type": args.fp_type,
        "radius": args.radius,
        "fp_bits": args.fp_bits,
        "chirality": bool(args.chirality),
        "scaffold_match": args.scaffold_match,
        "generic_scaffold": bool(args.generic_scaffold),
        "query_is_scaffold": bool(args.query_is_scaffold),
        "filters": parse_filters(args.filter),
        "value_col": args.value_col,
        "units_col": args.units_col,
        "normalize_units": norm_on,
        "norm_map": norm_map,
        "thresholds": thresholds,
        "has_thresholds": bool(thresholds),
    }

    try:
        qinfo = prepare_query(cfg)
    except ValueError as e:
        raise SystemExit("[query error] %s" % e)

    if not args.quiet:
        print("[query] " + describe_query(qinfo, cfg), file=sys.stderr)
        print("[mode ] %s | max_hits=%s | jobs=%d"
              % (args.mode, "unlimited" if max_hits == 0 else max_hits, args.jobs),
              file=sys.stderr)
        if thresholds:
            print("[thr  ] " + ", ".join(
                "%s: [%s, %s]" % (k, v[0], v[1]) for k, v in thresholds.items()),
                file=sys.stderr)

    stats, hits, header = run_search(cfg, args)

    if args.dedupe:
        best = {}
        for s, i, row, extra in hits:
            key = (row.get(cfg["id_col"]) or "").strip() or ("ROW_%d" % i)
            if key not in best or s > best[key][0]:
                best[key] = (s, i, row, extra)
        hits = list(best.values())

    score_mode = (cfg["mode"] == "similarity" or
                  (cfg["mode"] == "scaffold" and cfg["scaffold_match"] == "similar"))
    by_score = score_mode if args.sort == "auto" else (args.sort == "score")
    hits.sort(key=(lambda t: (-t[0], t[1])) if by_score else (lambda t: t[1]))

    out_fields = list(header)
    extra_fields = (list(NORM_EXTRA_FIELDS) if norm_on else []) + \
                   list(SEARCH_EXTRA_FIELDS)
    for f in extra_fields:
        if f not in out_fields:
            out_fields.append(f)

    with open(args.output, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=out_fields, extrasaction="ignore")
        w.writeheader()
        for _s, _i, row, extra in hits:
            rec = dict(row)
            rec.update(extra)
            w.writerow(rec)

    if not args.quiet:
        msg = ("[done ] scanned=%d | filtered=%d | invalid_smiles=%d | "
               "matched=%d | written=%d"
               % (stats["total"], stats["filtered"], stats["invalid"],
                  stats["matched"], len(hits)))
        if norm_on:
            msg += " | norm_failed=%d" % stats["norm_failed"]
        if thresholds:
            msg += " | other_type=%d" % stats["other_type"]
        print(msg + " -> " + args.output, file=sys.stderr)


if __name__ == "__main__":
    main()
