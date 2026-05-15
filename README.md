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

## 用mayachemtools计算单点能

命令行：
`Psi4CalculateEnergy.py -i test.sdf --ov -o test_energy.sdf --methodName r2scan-3c --basisSet DEF2-mTZVPP --psi4DDXSolvation yes --psi4DDXSolvationParams "solvent water" --mp NO --psi4RunParams "NumThreads, 16"`

如果对多个分子进行并行计算，则使用`--mp YES`

检查计算结果是否收敛，了解关键信息：
`psi4_analysis.py test_Psi4.out`

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
