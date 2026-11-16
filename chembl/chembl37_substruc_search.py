#!/usr/bin/env python
import argparse
import csv
import pickle
from rdkit import Chem

def main():
    parser = argparse.ArgumentParser(description="用 SMARTS 在 ChEMBL parent 子结构库中做子结构搜索")
    parser.add_argument("-query", required=True, help="查询 SMARTS")
    parser.add_argument("-o", required=True, help="输出 CSV 文件路径")
    parser.add_argument("--lib", default="/public/gkxiao/software/chembl/chembl37_parent_ssslib.pkl", help="子结构库 pickle 文件路径")
    parser.add_argument("--records", default="chembl37_parent_fps.pkl", help="元数据 records pickle 文件路径")
    parser.add_argument("--max", type=int, default=None, help="最多输出多少条命中，默认输出全部")
    parser.add_argument("--chirality", action="store_true", help="启用立体化学匹配（默认忽略手性）")
    args = parser.parse_args()

    query = Chem.MolFromSmarts(args.query)
    if query is None:
        raise SystemExit(f"无法解析 SMARTS: {args.query}")

    with open(args.lib, "rb") as f:
        library = pickle.load(f)

    match_indices = library.GetMatches(query, useChirality=args.chirality)

    with open(args.records, "rb") as f:
        records = pickle.load(f)

    if args.max is not None:
        match_indices = match_indices[:args.max]

    hits = []
    for idx in match_indices:
        r = records[idx]
        hits.append({
            "molregno": r["molregno"],
            "chembl_id": r["chembl_id"],
            "smiles": r["smiles"],
        })

    with open(args.o, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["molregno", "chembl_id", "smiles"])
        writer.writeheader()
        writer.writerows(hits)

    print(f"命中 {len(match_indices)} 条，已写入 {args.o}")

if __name__ == "__main__":
    main()
