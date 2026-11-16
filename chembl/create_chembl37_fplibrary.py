import csv
import pickle
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

CSV_PATH = "chembl_37_parent_molecules.csv"
PKL_PATH = "chembl_parent_fps.pkl"

morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

records = []
skipped = 0

with open(CSV_PATH, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if mol is None:
            skipped += 1
            continue
        records.append({
            "molregno": int(row["molregno"]),
            "chembl_id": row["chembl_id"],
            "smiles": row["canonical_smiles"],
            "fp": morgan_gen.GetFingerprint(mol),
        })

with open(PKL_PATH, "wb") as f:
    pickle.dump(records, f)

print(f"已保存 {len(records)} 条指纹，跳过无效 SMILES {skipped} 条")
