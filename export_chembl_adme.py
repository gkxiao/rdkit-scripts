import sqlite3
import pandas as pd

# 1. 连接到你的 ChEMBL 35 SQLite 数据库
db_path = "/public/gkxiao/software/chembl/chembl_35_sqlite3.db"
conn = sqlite3.connect(db_path)

# 2. 编写 SQL 查询语句
sql_query = """
SELECT
    md.chembl_id AS compound_chembl_id,
    cs.canonical_smiles,
    act.standard_type,
    act.standard_value,
    act.standard_units,
    act.standard_relation,
    ast.chembl_id AS assay_chembl_id,
    ast.description AS assay_description,
    ast.assay_type
FROM
    activities act
JOIN
    assays ast ON act.assay_id = ast.assay_id
JOIN
    molecule_dictionary md ON act.molregno = md.molregno
JOIN
    compound_structures cs ON md.molregno = cs.molregno
WHERE
    ast.assay_type = 'A'
    AND LOWER(act.standard_type) IN ('t1/2', 'clint', 'cl', 'remaining', '% remaining', '%remaining');
"""

print("正在从 ChEMBL 35 SQLite 数据库中执行查询...")
# 使用 pandas 读取并直接执行
df = pd.read_sql_query(sql_query, conn)

# 3. 关闭连接
conn.close()

# 4. 保存为 CSV
output_file = "chembl_35_adme_stability.csv"
df.to_csv(output_file, index=False, encoding='utf-8')
print(f"查询完成！共获取 {len(df)} 条记录，已成功保存至: {output_file}")
