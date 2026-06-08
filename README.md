# rdkit-scripts

## 计算PPI界面残基的SASA变化

通过以下方式调用：
`python sasa.py input.pdb A 89`

其中，input.pdb含有A、B两个Chain，计算复合物与单体状态时指定残基的&Delta;SASA。

---

### 使用说明

1.  **基本运行**：
    ```bash
    python sasa_calc.py 1IAR_prepared.pdb A 89
    ```
2.  **批量处理示例**：
    如果你想遍历一个包含多个残基的列表（例如残基 89, 90, 91），可以使用简单的 Shell 循环：
    ```bash
    for res in 89 90 91; do python sasa_calc.py protein.pdb A $res; done
    ```

### 提示：
如果你需要处理大规模的虚筛结果或长程 MD 轨迹的每一帧，可以将 `Chem.MolFromPDBFile` 移到循环外，仅在内存中操作 `EditableMol`，这样可以极大减少磁盘 I/O 开销。

## 用xTB进行几何优化

命令行：

```
# 使用 -o 参数直接指定输出文件（如果 xtb 支持）
# 或者直接使用 namespace 生成的文件名

# XTB 优化（使用 namespace）
xtb test.sdf \
  --namespace CONF_1 \
  --opt tight \
  --alpb water \
  --gfn 2 \
  -c 0 \
  -u 0 \
  --parallel 8
```
优化后新的文件目录如下
```
.
├── test.sdf                    # 原始输入
├── CONF_1.xtbopt.sdf           # 优化后的几何（用于后续计算）
├── CONF_1.xtbopt.log           # 优化日志
├── CONF_1.charges              # 电荷
├── CONF_1.wbo                  # Wiberg键级
├── CONF_1.xtbrestart           # 重启文件
└── CONF_1.xtbtopo.sdf          # 拓扑文件
```

其中，`CONF_1.xtbopt.sdf`是优化过的结果文件。该文件不适合直接使用, 需要转化为合理的SDF格式：

```
python fix_xtb_sdf.py CONF_1.xtbopt.sdf -t "CONF_1" -o CONF_1_opt_fixed.sdf
```

## 计算单点能

用`MayaChemTools`的Psi4工具进行单点能计算，命令行如下：
```shell
Psi4CalculateEnergy.py -i CONF_1_opt_fixed.sdf --ov -o CONF_1_spe.sdf \
  --methodName r2scan-3c \
  --basisSet DEF2-mTZVPP \
  --psi4DDXSolvation yes \
  --psi4DDXSolvationParams "solvent water" \
  --mp NO \
  --psi4RunParams "NumThreads, 16"
```
如果对多个分子进行并行计算，则使用`--mp YES`

检查计算结果是否收敛，了解关键信息：
```shell
psi4_analysis.py CONF_1_opt_fixed_Psi4.out
```

结果如下：

```
======================================================================
Psi4 计算结果分析: CONF_1_opt_fixed_Psi4.out
======================================================================

✅ SCF 收敛成功

📚 计算方法: r2SCAN-3c
📚 基组: DEF2-MTZVPP
💻 计算资源: 16 线程, 1.0 GB 内存

⚡ 单点能 (总能量):
   -768.94201175 Hartree
   -482703.744 kcal/mol

💧 溶剂化能:
   -0.00853204 Hartree
   -5.354 kcal/mol

🔌 偶极矩:
   X 分量: -0.0442 a.u. = -0.1123 Debye
   Y 分量: -0.6831 a.u. = -1.7361 Debye
   Z 分量: 0.0970 a.u. = 0.2465 Debye
   总大小: 0.691 Debye

🔄 SCF 迭代次数: 16

⏱️  计算时间:
   总时间 (wall time): 46.0 秒 (0.77 分钟)
   CPU 时间 (user): 662.7 秒 (11.04 分钟)
   系统时间 (system): 19.8 秒

======================================================================
```

## 构象系综的相对能量计算

根据`Psi4CalculateEnergy.py`计算的单点能SDF结果文件中的`Psi4_Energy (kcal/mol)`进一步计算相对能量：
```
usage: calc_rel_energy.py [-h] -i INPUT -o OUTPUT
                          [--prop PROP]
                          [--outprop OUTPROP]

Calculate relative energies for conformers in an SDF file using a specified energy property.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input multi-conformer SDF file
  -o OUTPUT, --output OUTPUT
                        Output SDF file
  --prop PROP
                        Energy property name in SDF
                        Default: "Psi4_Energy (kcal/mol)"
  --outprop OUTPROP
                        Output relative energy property name
                        Default: "Psi4_Rel_Energy"

Example:
  python calc_rel_energy.py -i input.sdf -o output.sdf
```

也可以根据指定的性质列计算相对能量，并指定能量单位：

Psi4(kcal/mol)

```bash
python calc_rel_energy.py \
    -i input.sdf \
    -o output.sdf \
    --prop "Psi4_Energy (kcal/mol)" \
    --unit kcal \
    --outprop "Psi4_Rel_Energy"
```

xTB（Hartree转kcal/mol）
```bash
python calc_rel_energy.py \
    -i input.sdf \
    -o output.sdf \
    --prop "Energy_xTB" \
    --unit hartree \
    --outprop "Rel_Energy_xTB"
```

## 用Crest补充构象系综

常规的方法很多时候会遗漏重要构象，`crest`可以作为补充方法，合并多种来源的构象后进行分析：
```
crest input.xyz --v3 \
--gfn2 \
--chrg 0 \
--uhf 0 \
--ewin 10 \
--mrest 10 \
--T 16 \
```

其中:
- `--mrest 10`增加 metadynamics restart 次数更容易找到隐藏构象
- `-T 16` 使用 16 线程


## 从XYZ转SDF

量化计算最常见的一个问题是，如何将XYZ转化为SDF，并归属正确的键类型、原子类型与Formal Charge。
假设你的计算是从一个SDF文件（start.sdf）开始，在进行QM计算时使用了从这个SDF而转化得到的start.xyz, 计算之后得到xtbopt.xyz：

```
RDKit_xyz2sdf.py -i start.sdf -x xtbopt.xyz -o xtbopt.sdf
```

这个脚本可以保留原有的拓扑，而将坐标更新为优化后的坐标（从xyz文件读入）来实现格式转化，这可以确保结构正确。

## Schrodinger
1. 构象搜索
```
sch_mmod_csearch.py -i 4zlz_ligand.maegz -o 1_mmod_csearch -m LMCS --ff opls2005 --run
```
2. 构象合并
```
cd 1_mmod_csearch
sch_confs_merge.py -i `ls *-out.maegz` -o 4zlz_ligand_confs.maegz
$SCHRODINGER/utilities/structconvert 4zlz_ligand_confs.maegz 4zlz_ligand_confs.sdf
```
3. xtb几何优化
```
mkdir 2_xtbopt
$SCHRODINGER/utilities/obabel -isdf 4zlz_ligand_confs.sdf -oxyz -O 2_xtbopt/CONF_.xyz -m
cd 2_xtbopt
n=`ls *.xyz|wc -l`
for i in `seq 1 ${n}`
do
xtb CONF_${i}.xyz --opt tight -c 0 -u 0 --alpb water --namespace CONF_${i}
cat CONF_${i}.xtbopt.xyz >> ensemble.xtbopt.xyz
done
```
得到构象系综`ensemble.xtbopt.xyz`，接下来要进行构象聚类


