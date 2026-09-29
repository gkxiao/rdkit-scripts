# rdkit-scripts

## 计算PPI界面残基的SASA变化

通过以下方式调用：
`python sasa_calc.py input.pdb A 89`

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
其中，xtb可以使用$SCHRODINGER/run xtb来代替。
这个优化过程，可以写成脚本`xtbopt_batch.sh`来实现。最后得到构象系综`ensemble.xtbopt.xyz`，接下来要进行构象聚类。

## 混合精度的构象系综分析

混合精度构象系综采样与分析采用多级别能量评估方案：初始阶段依靠低精度力场快速筛除不可行构象；随后通过半经验量子化学（SQM）或其它低计算的方法实施几何弛豫与构象去重；最后执行 DFT 单点能计算，依据高精度能量对构象系综重新排序。在这个过程中，需要考虑两个问题：

### 正在进行DFT计算：早期停止评估

需要评估继续计算的潜在收益，停止计算的潜在损失，剩余全量计算的成本，当前DFT计算的系综已经覆盖了多大的空间？

输入类似：
```
CREST/xTB：1000 个 conformers
DFT SPE：100 / 1000 已完成
```
此时关心的是：

- 继续算下去还有多少潜在价值？
- 剩余 conformers 中还可能有多少个 DFT 目标？
- 如果现在停止，漏掉目标的风险多大？
- 继续计算还需要多少资源？
- 当前已经完成的 DFT 结果，对 xTB 预筛选的支持程度如何？

### 已经全量完全DFT计算：构象覆盖度评估

输入类似：
```
CREST/xTB：1000
DFT SPE：1000 / 1000
```
此时关心的是：

- 初始的构象系综是否覆盖了DFT低能构象空间？
- 原始 CREST/xTB ensemble 是否足够？
- 如果不够，是 DFT 计算不够，还是更根本的起始构象（力场阶段、SQM阶段）空间不够？

这两个问题构成一个闭环:

```
                 CREST / xTB
                     │
                     ▼
             Initial Ensemble
                 Coverage
                     │
                     ▼
              DFT calculations
                     │
          ┌──────────┴──────────┐
          │                     │
       中途状态                完成状态
          │                     │
          ▼                     ▼
   Early Stopping          Final Coverage
       Analysis              Analysis
          │                     │
          └──────────┬──────────┘
                     ▼
             下一轮 CREST 参数
             / DFT 计算策略优化

```

### 示例

Flare进行初始的基于力场的构象搜索、Flare Ligand QM进行几何优化(xTB GFN2)、构象去重，得到起始的构象系综（xtb_ensemble.sdf）, 接着在R2SCAN-3c理论水平进行单点能计算（ensemble_spe.sdf）, 现在要评估DFT ewin=3 kcal/mol 构象系综。

```bash
./dft_ensemble_ana.py --crest xtb_ensemble.sdf --spe ensemble_spe.sdf --ewin 3

======================================================================
Input Ensemble Summary
======================================================================
Total conformers                  : 25
xTB energy window                 : 7.78 kcal/mol

Conformers <= 1 kcal/mol        : 6
Conformers <= 2 kcal/mol        : 6
Conformers <= 3 kcal/mol        : 13
Conformers <= 4 kcal/mol        : 13
Conformers <= 5 kcal/mol        : 13
Conformers <= 6 kcal/mol        : 20

======================================================================
DFT Progress
======================================================================
DFT completed                     : 25 / 25
Completion                        : 100.0%
Current Frontier                  : 7.78 kcal/mol
Remaining Search Space            : 0.00 kcal/mol

======================================================================
Recovery Analysis
======================================================================
Target Window                     : 3.00 kcal/mol
DFT Hits                          : 11
Recovery Frontier                 : 2.55 kcal/mol
Coverage Beyond Recovery Frontier : 5.23 kcal/mol
Consecutive Misses                : 14

======================================================================
Recovery Efficiency Analysis
======================================================================
DFT Window Compression Ratio      : 0.846

Recovery Efficiency
--------------------------------------------------
xTB <=  1 kcal/mol :    6/11   ( 54.5%)
xTB <=  2 kcal/mol :    6/11   ( 54.5%)
xTB <=  3 kcal/mol :   11/11   (100.0%)
xTB <=  4 kcal/mol :   11/11   (100.0%)
xTB <=  5 kcal/mol :   11/11   (100.0%)
xTB <=  6 kcal/mol :   11/11   (100.0%)
xTB <=  7 kcal/mol :   11/11   (100.0%)
xTB <=  8 kcal/mol :   11/11   (100.0%)

Among DFT Hits
--------------------------------------------------
90th percentile E_xTB_Rel        : 2.49 kcal/mol
95th percentile E_xTB_Rel        : 2.52 kcal/mol
Maximum E_xTB_Rel                : 2.55 kcal/mol

======================================================================
Coverage Risk Assessment
======================================================================
Recovery Frontier Occupancy      : 32.7%
Coverage Risk Level              : LOW

======================================================================
Ranking Consistency
======================================================================
Pearson R²                        : 0.970
Spearman rho                      : 0.988

======================================================================
Bayesian Recovery Predictor
======================================================================
Suggested Bayesian Cutoff         : 6.27 kcal/mol
Probability Threshold             : 1.0%
Expected Remaining Hits           : 0.00

======================================================================
Recommendation
======================================================================
DFT calculations are essentially complete.

Current xTB search window appears sufficient.
```

我个人最喜欢`Recovery Analysis`这个部分。在 3 kcal/mol 的 DFT target window 内，目前发现了 11 个 target；最后一个 target 位于 xTB 相对能量 2.55 kcal/mol。此后，沿 xTB 能量升高方向又观察了 14 个连续的 DFT non-target conformers，覆盖了额外 5.23 kcal/mol 的 xTB 能量范围，但没有产生新的 DFT target。而对于正在进行DFT计算的项目，这不是证明“后面没有 target”，而是根据目前观察到的结果，判断继续计算的边际收益是否已经越来越小。

在`Coverage Risk Assessment`部分, 给出Coverage Risk Level的评估是“Low”，Low这种等级依据何来？是否合理需要进一步权衡。

注意：`Bayesian Recovery Predictor`这部分，在代码里是`LogisticRegression`，严格来说不是`Bayesian model`，虽然在过程上采用了“Bayesian/probabilistic decision thinking”。也许改为`Bayesian-Inspired Recovery Predictor`更合适：不是 Bayesian inference，但采用了 Bayesian-style probabilistic decision framework。属于模型外推，而不是 Recovery Analysis 的直接观测结果。
