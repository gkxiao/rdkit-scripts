import sqlite3
import csv

DB_PATH = "/public/gkxiao/chembl/chembl_37_sqlite3.db"
OUTPUT_CSV = "chembl_37_parent_molecules.csv"

conn = sqlite3.connect(DB_PATH)
cursor = conn.cursor()

query = """
SELECT
    md.molregno          AS parent_molregno,
    cs.canonical_smiles  AS parent_smiles,
    md.chembl_id         AS parent_chembl_id
FROM molecule_hierarchy mh
JOIN molecule_dictionary md
    ON mh.parent_molregno = md.molregno
JOIN compound_structures cs
    ON md.molregno = cs.molregno
WHERE mh.molregno = mh.parent_molregno
  AND cs.canonical_smiles IS NOT NULL
ORDER BY md.molregno
"""

cursor.execute(query)
rows = cursor.fetchall()

with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["molregno", "smiles", "chembl_id"])
    writer.writerows(rows)

print(f"已导出 {len(rows)} 条 parent 分子记录至 {OUTPUT_CSV}")

cursor.close()
conn.close()
