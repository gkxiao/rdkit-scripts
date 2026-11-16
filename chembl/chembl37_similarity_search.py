#!/usr/bin/env python
import argparse
import csv
import pickle
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

def main():
    parser = argparse.ArgumentParser(description="用 Morgan 指纹在 ChEMBL parent 库中做相似性搜索")
    parser.add_argument("-query", required=True, help="查询分子的 SMILES")
    parser.add_argument("-o", required=True, help="输出 CSV 文件路径")
    parser.add_argument("--pkl", default="/public/gkxiao/software/chembl/chembl37_parent_fps.pkl", help="预计算指纹 pickle 文件路径")
    parser.add_argument("--threshold", type=float, default=0.7, help="Tanimoto 相似度阈值，默认 0.7")
    parser.add_argument("--topn", type=int, default=None, help="只输出前 N 条（按相似度降序），默认输出全部命中")
    args = parser.parse_args()

    # 指纹参数必须与保存时一致
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    query_mol = Chem.MolFromSmiles(args.query)
    if query_mol is None:
        raise SystemExit(f"无法解析查询 SMILES: {args.query}")
    query_fp = morgan_gen.GetFingerprint(query_mol)

    with open(args.pkl, "rb") as f:
        records = pickle.load(f)
    fps = [r["fp"] for r in records]

    sims = DataStructs.BulkTanimotoSimilarity(query_fp, fps)

    hits = []
    for i, sim in enumerate(sims):
        if sim >= args.threshold:
            hits.append({
                "score": sim,
                "molregno": records[i]["molregno"],
                "chembl_id": records[i]["chembl_id"],
                "smiles": records[i]["smiles"],
            })

    hits.sort(key=lambda x: x["score"], reverse=True)
    if args.topn is not None:
        hits = hits[:args.topn]

    with open(args.o, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["score", "molregno", "chembl_id", "smiles"])
        writer.writeheader()
        writer.writerows(hits)

    print(f"阈值 {args.threshold} 下命中 {len(hits)} 条，已写入 {args.o}")

if __name__ == "__main__":
    main()
