import pickle
from rdkit import Chem
from rdkit.Chem import rdSubstructLibrary

# 复用你已有的 records（含 smiles, molregno, chembl_id）
with open("chembl37_parent_fps.pkl", "rb") as f:
    records = pickle.load(f)

# 创建分子容器和 Pattern 指纹容器
mol_holder = rdSubstructLibrary.CachedTrustedSmilesMolHolder()
fp_holder = rdSubstructLibrary.PatternHolder(2048)  # Pattern 指纹，非 Morgan

for r in records:
    smi = r["smiles"]
    # 添加 SMILES（Trusted 模式要求是 RDKit 规范化的 SMILES）
    idx_mol = mol_holder.AddSmiles(smi)
    # 生成并添加 PatternFingerprint
    mol = Chem.MolFromSmiles(smi)
    fp = Chem.PatternFingerprint(mol, fpSize=2048)
    idx_fp = fp_holder.AddFingerprint(fp)
    assert idx_mol == idx_fp  # 索引必须对齐

# 组装 SubstructLibrary
library = rdSubstructLibrary.SubstructLibrary(mol_holder, fp_holder)

# 保存：直接用 pickle
with open("chembl37_parent_ssslib.pkl", "wb") as f:
    pickle.dump(library, f)

print(f"已保存子结构库，包含 {len(library)} 个分子")
