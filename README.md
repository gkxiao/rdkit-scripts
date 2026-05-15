# rdkit-scripts

## 计算PPI界面残基的SASA变化

通过以下方式调用：
`python sasa.py input.pdb A 89`

其中，input.pdb含有A、B两个残基，计算复合物与单体状态时指定残基的&Delta;SASA。

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
├── CONF_1.xtbtopo.sdf          # 拓扑文件
└── CONF_1_energy.sdf           # Psi4单点能结果
···

其中，`CONF_1.xtbopt.sdf`是优化过的结果文件。该文件不适合直接使用, 需要转化为合理的SDF格式：

```
python fix_xtb_sdf.py CONF_1.xtbopt.sdf -t "CONF_1" -o CONF_1_opt_fixed.sdf
```

## 用mayachemtools计算单点能

命令行：
```shell
Psi4CalculateEnergy.py -i CONF_1_opt_fixed.sdf --ov -o test_spe.sdf \
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
psi4_analysis.py test_Psi4.out
```

结果如下：

```
======================================================================
Psi4 计算结果分析: test_Psi4.txt
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
