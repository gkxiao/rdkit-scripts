## 导出母体化合物

从SQLite3数据库里将母体化合物的登记号，SMILES与CHEMBL_ID导出为CSV格式，以便后续查询。

```bash
python export_chembl_smiles.py
```

## 创建搜索用数据库

```bash
# 创建指纹相似性搜索数据库
python create_chembl37_fplibrary.py

# 创建子结构数搜索据库
python create_chembl37_substructlibrary.py
```

## 相似性用法示例

```bash
# 基本用法
python chembl37_similarity_search.py -query 'Cc1cc(-n2ncc(=O)[nH]c2=O)ccc1C(=O)c1ccccc1Cl' -o hits.csv

# 调整阈值并只输出前 20 条
python chembl37_similarity_search.py -query 'CCO' -o hits.csv --threshold 0.6 --topn 20
```

### 几个实现上的注意点

- **SMILES 里的特殊字符要加引号**：`-query` 的值如果包含 `#`、`=` 等字符，在 shell 里必须用单引号包住，否则会被 shell 解释。上面的示例已经加了单引号。
- **输出顺序**：`hits` 先按阈值过滤再排序，`--topn` 在排序之后截断，所以拿到的是全局最相似的前 N 条，而不是先截断再排序。
- **指纹参数一致性**：脚本里 `radius=2, fpSize=2048` 必须和生成 `chembl37_parent_fps.pkl` 时完全一致，否则相似度无意义。如果以后改了参数，建议把参数也写进 pickle 文件名，避免混用。
- **无效查询**：`MolFromSmiles` 返回 `None` 时直接退出并提示，避免后续 `GetFingerprint` 报错。


## 子结构搜索用法示例

```bash
# 基本用法：查苯环
python chembl37_substruct_search.py -query 'c1ccccc1' -o hits.csv

# 查酰胺，限制最多 100 条，并启用手性匹配
python chembl37_substruct_search.py -query '[CX3](=O)[NX3]' -o hits.csv --max 100 --chirality
```

### 几个实现上的注意点

- **SMARTS 要加单引号**：`-query` 的值里常含 `[`、`]`、`=`、`#` 等字符，shell 会解释它们，必须用单引号包住，否则会报错或参数被截断。
- **索引对齐依赖顺序**：`match_indices` 是 `SubstructLibrary` 内部的索引，回查 `records[idx]` 能正确对应，前提是构建库时按 `records` 的同一顺序遍历添加。如果之后对 `records` 重新排序或过滤，索引就会错位。
- **`--chirality` 默认关闭**：RDKit 子结构匹配默认忽略立体化学，需要区分手性时加这个开关。注意它只影响匹配阶段，如果构建库时分子已丢失手性信息，开启也无济于事。
- **`--max` 在匹配后截断**：先拿到全部命中索引再截断，顺序是 `SubstructLibrary` 返回的自然顺序，不是按相似度排序（子结构搜索没有相似度分数）。如果需要特定排序，可以后续再处理。
- **输出无 score 列**：子结构搜索是布尔匹配，没有 Tanimoto 分数，所以 CSV 里只有 `molregno`、`chembl_id`、`smiles` 三列。
