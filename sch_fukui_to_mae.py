#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fukui_to_mae.py

把 Jaguar/fukui 主输出 (*.out) 里的 "Atomic Fukui indices" 原子表
写成结构文件的 Maestro 原子属性 (r_ 前缀, 可在 Project Table 直接排序/上色)。

    f- = HOMO 块 f_NN   -> r_fukui_fminus
    f+ = LUMO 块 f_NN   -> r_fukui_fplus
    f0 = (f- + f+)/2    -> r_fukui_f0

用法:
    $SCHRODINGER/run fukui_to_mae.py \
        -fukui-out CONF_17_1.out \
        -imae      CONF_17.mae \
        -omae      CONF_17-out.mae
"""

import argparse
import re
import sys

LABEL_RE = re.compile(r"^([A-Z][a-z]?)(\d+)$")


def parse_fukui_table(path):
    """解析 'Atomic Fukui indices' 表, 返回 {atom_index: {...}}"""
    try:
        with open(path, "r", errors="replace") as fh:
            lines = fh.read().splitlines()
    except OSError as exc:
        sys.exit("ERROR: 打不开 %s: %s" % (path, exc))

    start = next((i for i, l in enumerate(lines)
                  if "Atomic Fukui indices" in l), None)
    if start is None:
        sys.exit("ERROR: %s 中找不到 'Atomic Fukui indices' 表" % path)

    table = {}
    for line in lines[start + 1:]:
        s = line.strip()
        if not s:                     # 表结束
            if table:
                break
            continue
        parts = s.split()
        # 数据行: Atom  f_NN f_NS f_SN f_SS | f_NN f_NS f_SN f_SS  (9 列)
        if len(parts) != 9:
            continue
        m = LABEL_RE.match(parts[0])
        if not m:
            continue
        try:
            v = [float(x) for x in parts[1:]]
        except ValueError:
            continue
        idx = int(m.group(2))
        f_minus, f_plus = v[0], v[4]      # 两个块的 f_NN
        table[idx] = {
            "elem":    m.group(1),
            "fminus":  f_minus,
            "fplus":   f_plus,
            "f0":      0.5 * (f_minus + f_plus),
            "f0_tab":  v[3],             # 表内 f_SS, 用于对账
        }

    if not table:
        sys.exit("ERROR: 找到表但没解析出原子行, 请检查 %s 的格式" % path)
    return table


def main():
    ap = argparse.ArgumentParser(description="把 fukui 原子指数写进 mae")
    ap.add_argument("-fukui-out", dest="fukui_out", required=True,
                    help="含 Atomic Fukui indices 表的 *.out")
    ap.add_argument("-imae", dest="imae", required=True, help="输入结构")
    ap.add_argument("-omae", dest="omae", required=True, help="输出结构")
    args = ap.parse_args()

    try:
        from schrodinger import structure
    except ImportError:
        sys.exit("ERROR: 找不到 schrodinger 模块, 请用 "
                 "$SCHRODINGER/run fukui_to_mae.py ... 运行")

    table = parse_fukui_table(args.fukui_out)

    st = next(structure.StructureReader(args.imae), None)
    if st is None:
        sys.exit("ERROR: %s 里没有结构" % args.imae)

    # --- 逐原子写入 ---
    missing, mismatch = [], []
    for at in st.atom:
        rec = table.get(at.index)
        if rec is None:
            missing.append(at.index)
            continue
        if at.element.lower() != rec["elem"].lower():
            mismatch.append("%d(%s!=%s)" % (at.index, at.element, rec["elem"]))
        at.property["r_fukui_fminus"] = rec["fminus"]
        at.property["r_fukui_fplus"]  = rec["fplus"]
        at.property["r_fukui_f0"]     = rec["f0"]

    extra = sorted(set(table) - {a.index for a in st.atom})

    # --- 对账: 表内 f_SS 应 ≈ (f-+f+)/2 ---
    diff = [abs(r["f0"] - r["f0_tab"]) for r in table.values()]
    if diff and max(diff) > 1e-3:
        print("WARN: f0 与表内 f_SS 最大差 %.4g, 请确认列序" % max(diff))

    if missing:
        print("WARN: %d 个原子表中无对应, 以 -1 填充: %s"
              % (len(missing), missing[:20]))
        for i in missing:
            st.atom[i].property["r_fukui_fminus"] = -1.0
            st.atom[i].property["r_fukui_fplus"]  = -1.0
            st.atom[i].property["r_fukui_f0"]     = -1.0
    if mismatch:
        print("WARN: 元素符号不一致(原子顺序可能没对齐): " + ", ".join(mismatch))
    if extra:
        print("WARN: 表中有结构中不存在的原子编号: %s" % extra)

    # --- 写文件 ---
    with structure.StructureWriter(args.omae) as writer:
        writer.append(st)

    # --- 摘要 ---
    def top(key):
        return sorted(table.items(), key=lambda kv: -kv[1][key])[:5]

    def fmt(items, key):
        return ", ".join("%s%d=%.3f" % (r["elem"], i, r[key]) for i, r in items)

    print("已写入 %d 个原子的 f-/f+/f0 -> %s" % (st.atom_total, args.omae))
    print("  f- 最大(亲电进攻): " + fmt(top("fminus"), "fminus"))
    print("  f+ 最大(亲核进攻): " + fmt(top("fplus"),  "fplus"))
    print("  f0 最大(自由基)  : " + fmt(top("f0"),     "f0"))
    print("完成。")


if __name__ == "__main__":
    main()
