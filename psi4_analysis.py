#!/usr/bin/env python3
"""
Psi4 结果文件分析脚本 - 修正版
修复偶极矩和计算时间提取问题
"""

import re
import sys
import argparse
from pathlib import Path
from typing import Dict, Optional, Tuple, List

class Psi4OutputAnalyzer:
    """Psi4 输出文件分析器 - 修正版"""

    def __init__(self, filename: str):
        self.filename = filename
        self.content = ""
        self.results = {}

    def read_file(self) -> bool:
        """读取输出文件"""
        try:
            with open(self.filename, 'r', encoding='utf-8', errors='ignore') as f:
                self.content = f.read()
            return True
        except FileNotFoundError:
            print(f"错误：文件 '{self.filename}' 不存在")
            return False
        except Exception as e:
            print(f"错误：读取文件失败 - {e}")
            return False

    def extract_method_and_basis(self) -> Tuple[Optional[str], Optional[str]]:
        """提取计算方法和基组"""
        method = None
        basis = None

        # 方法识别
        if re.search(r'r2SCAN', self.content, re.IGNORECASE):
            method = 'r2SCAN-3c'
        elif re.search(r'wb97x[-_]d3bj', self.content, re.IGNORECASE):
            method = 'WB97X-D3BJ'
        elif re.search(r'b3lyp', self.content, re.IGNORECASE):
            method = 'B3LYP'
        elif re.search(r'pbe0', self.content, re.IGNORECASE):
            method = 'PBE0'

        # 基组提取
        basis_match = re.search(r'Basis Set:\s+(\S+)', self.content)
        if basis_match:
            basis = basis_match.group(1)

        return method, basis

    def extract_total_energy(self) -> Optional[float]:
        """提取总能量 (Hartree)"""
        # 从 "Total Energy =" 提取
        pattern = r'Total Energy\s*=\s*([-\d.]+)'
        match = re.search(pattern, self.content)
        if match:
            return float(match.group(1))

        # 从 "DF-RKS Final Energy" 提取
        pattern2 = r'@DF-[A-Z]+\s+Final Energy:\s+([-\d.]+)'
        match = re.search(pattern2, self.content)
        if match:
            return float(match.group(1))

        return None

    def extract_solvation_energy(self) -> Optional[float]:
        """提取溶剂化能 (Hartree)"""
        # DD Solvation Energy
        pattern = r'DD Solvation Energy\s*=\s*([-\d.]+)'
        match = re.search(pattern, self.content)
        if match:
            return float(match.group(1))
        return None

    def check_convergence(self) -> Tuple[bool, str]:
        """检查计算是否收敛"""
        if re.search(r'Energy and wave function converged', self.content):
            return True, "SCF 收敛成功"
        elif re.search(r'SCF\s+failed|Maximum number of iterations exceeded', self.content, re.IGNORECASE):
            return False, "SCF 未收敛"
        elif self.extract_total_energy() is not None:
            return True, "计算完成"
        return False, "计算失败"

    def extract_iterations(self) -> Optional[int]:
        """提取SCF迭代次数"""
        iterations = re.findall(r'@DF-[A-Z]+\s+iter\s+(\d+):', self.content)
        if iterations:
            return int(iterations[-1])
        return None

    def extract_dipole_moment(self) -> Optional[Dict[str, float]]:
        """提取偶极矩 - 修正版"""
        dipole = {}

        # 查找 Multipole Moments 表格
        # 匹配 Dipole X, Y, Z 行，提取 Total (a.u.) 值
        # 格式: Dipole X            :        -80.6613884           80.6172283           -0.0441601
        # 第三个值是总偶极矩 (a.u.)

        # 方法1：从 Multipole Moments 表格提取
        dipole_pattern = r'Dipole\s+([XYZ])\s*:\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)'
        matches = re.findall(dipole_pattern, self.content)

        if matches:
            dipole_components = {}
            for axis, elec, nucl, total in matches:
                dipole_components[axis] = float(total)

            if dipole_components:
                dipole['x'] = dipole_components.get('X', 0)
                dipole['y'] = dipole_components.get('Y', 0)
                dipole['z'] = dipole_components.get('Z', 0)

                # 计算总偶极矩大小
                dipole['magnitude_au'] = (dipole['x']**2 + dipole['y']**2 + dipole['z']**2)**0.5
                # 转换为 Debye (1 a.u. = 2.5417464519 Debye)
                dipole['magnitude_debye'] = dipole['magnitude_au'] * 2.5417464519

                return dipole

        # 方法2：直接提取 Magnitude 行
        magnitude_match = re.search(r'Magnitude\s*:\s*([\d.]+)', self.content)
        if magnitude_match:
            dipole['magnitude_debye'] = float(magnitude_match.group(1))
            return dipole

        return None

    def extract_computation_time(self) -> Optional[Dict[str, float]]:
        """提取计算时间 - 修正版"""
        time_info = {}

        # 从 tstop() 信息中提取
        # 格式: total time  =         46 seconds =       0.77 minutes
        total_match = re.search(r'total time\s*=\s*([\d.]+)\s+seconds', self.content)
        if total_match:
            time_info['total_seconds'] = float(total_match.group(1))
            time_info['total_minutes'] = time_info['total_seconds'] / 60

        # 用户时间 (CPU time)
        user_match = re.search(r'user time\s*=\s*([\d.]+)\s+seconds', self.content)
        if user_match:
            time_info['user_seconds'] = float(user_match.group(1))
            time_info['user_minutes'] = time_info['user_seconds'] / 60

        # 系统时间
        sys_match = re.search(r'system time\s*=\s*([\d.]+)\s+seconds', self.content)
        if sys_match:
            time_info['system_seconds'] = float(sys_match.group(1))

        return time_info if time_info else None

    def extract_threads_memory(self) -> Tuple[Optional[int], Optional[float]]:
        """提取线程数和内存"""
        threads = None
        memory = None

        # 线程数
        thread_match = re.search(r'Threads:\s+(\d+)', self.content)
        if thread_match:
            threads = int(thread_match.group(1))

        # 内存 (GB)
        mem_match = re.search(r'Memory:\s+([\d.]+)\s+GiB', self.content)
        if mem_match:
            memory = float(mem_match.group(1))

        return threads, memory

    def extract_warnings_errors(self) -> List[str]:
        """提取警告和错误信息"""
        issues = []

        warnings = re.finditer(r'Warning:\s*(.+?)(?:\n|$)', self.content, re.IGNORECASE)
        for warning in warnings:
            issues.append(f"⚠️  {warning.group(1).strip()}")

        errors = re.finditer(r'Error:\s*(.+?)(?:\n|$)', self.content, re.IGNORECASE)
        for error in errors:
            issues.append(f"❌  {error.group(1).strip()}")

        return issues

    def analyze(self) -> Dict:
        """执行完整分析"""
        if not self.read_file():
            return {}

        method, basis = self.extract_method_and_basis()
        total_energy = self.extract_total_energy()
        solvation_energy = self.extract_solvation_energy()
        converged, conv_msg = self.check_convergence()
        iterations = self.extract_iterations()
        dipole = self.extract_dipole_moment()
        comp_time = self.extract_computation_time()
        threads, memory = self.extract_threads_memory()
        issues = self.extract_warnings_errors()

        # 转换能量单位
        energy_kcal = None
        if total_energy is not None:
            energy_kcal = total_energy * 627.509474

        solvation_kcal = None
        if solvation_energy is not None:
            solvation_kcal = solvation_energy * 627.509474

        self.results = {
            'filename': self.filename,
            'method': method,
            'basis_set': basis,
            'total_energy_hartree': total_energy,
            'total_energy_kcal': energy_kcal,
            'solvation_energy_hartree': solvation_energy,
            'solvation_energy_kcal': solvation_kcal,
            'converged': converged,
            'convergence_message': conv_msg,
            'scf_iterations': iterations,
            'dipole_moment': dipole,
            'computation_time': comp_time,
            'threads': threads,
            'memory_gb': memory,
            'warnings_errors': issues
        }

        return self.results

    def print_summary(self):
        """打印分析结果摘要"""
        if not self.results:
            print("没有可用的分析结果")
            return

        print("\n" + "="*70)
        print(f"Psi4 计算结果分析: {Path(self.filename).name}")
        print("="*70)

        # 收敛状态
        if self.results['converged']:
            print(f"\n✅ {self.results['convergence_message']}")
        else:
            print(f"\n❌ {self.results['convergence_message']}")

        # 方法和基组
        if self.results['method']:
            print(f"\n📚 计算方法: {self.results['method']}")
        if self.results['basis_set']:
            print(f"📚 基组: {self.results['basis_set']}")

        # 资源信息
        if self.results['threads']:
            print(f"💻 计算资源: {self.results['threads']} 线程, {self.results['memory_gb']:.1f} GB 内存")

        # 能量信息
        if self.results['total_energy_hartree'] is not None:
            print(f"\n⚡ 单点能 (总能量):")
            print(f"   {self.results['total_energy_hartree']:.8f} Hartree")
            print(f"   {self.results['total_energy_kcal']:.3f} kcal/mol")

        if self.results['solvation_energy_hartree'] is not None:
            print(f"\n💧 溶剂化能:")
            print(f"   {self.results['solvation_energy_hartree']:.8f} Hartree")
            print(f"   {self.results['solvation_energy_kcal']:.3f} kcal/mol")

        # 偶极矩 - 修正显示
        if self.results['dipole_moment']:
            dipole = self.results['dipole_moment']
            print(f"\n🔌 偶极矩:")
            if 'x' in dipole:
                print(f"   X 分量: {dipole['x']:.4f} a.u. = {dipole['x'] * 2.5417464519:.4f} Debye")
                print(f"   Y 分量: {dipole['y']:.4f} a.u. = {dipole['y'] * 2.5417464519:.4f} Debye")
                print(f"   Z 分量: {dipole['z']:.4f} a.u. = {dipole['z'] * 2.5417464519:.4f} Debye")
            if 'magnitude_debye' in dipole:
                print(f"   总大小: {dipole['magnitude_debye']:.3f} Debye")
            elif 'magnitude_au' in dipole:
                print(f"   总大小: {dipole['magnitude_au']:.4f} a.u. = {dipole['magnitude_au'] * 2.5417464519:.3f} Debye")

        # SCF 迭代
        if self.results['scf_iterations']:
            print(f"\n🔄 SCF 迭代次数: {self.results['scf_iterations']}")

        # 计算时间 - 修正显示
        if self.results['computation_time']:
            ct = self.results['computation_time']
            print(f"\n⏱️  计算时间:")
            if 'total_seconds' in ct:
                total_min = ct['total_seconds'] / 60
                print(f"   总时间 (wall time): {ct['total_seconds']:.1f} 秒 ({total_min:.2f} 分钟)")
            if 'user_seconds' in ct:
                user_min = ct['user_seconds'] / 60
                print(f"   CPU 时间 (user): {ct['user_seconds']:.1f} 秒 ({user_min:.2f} 分钟)")
            if 'system_seconds' in ct:
                print(f"   系统时间 (system): {ct['system_seconds']:.1f} 秒")

        # 警告和错误
        if self.results['warnings_errors']:
            print(f"\n⚠️  发现 {len(self.results['warnings_errors'])} 个问题:")
            for issue in self.results['warnings_errors'][:5]:
                print(f"   {issue}")

        print("\n" + "="*70 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='分析 Psi4 输出文件，提取能量、收敛状态等信息',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('files', nargs='+', help='Psi4 输出文件')
    parser.add_argument('-o', '--output', help='输出CSV文件')
    parser.add_argument('-v', '--verbose', action='store_true', help='详细输出')

    args = parser.parse_args()

    # 处理文件列表
    files = []
    for pattern in args.files:
        path = Path(pattern)
        if path.is_file():
            files.append(str(path))
        else:
            matched = list(Path('.').glob(pattern))
            files.extend([str(f) for f in matched])

    if not files:
        print("错误：未找到匹配的文件")
        sys.exit(1)

    for file in files:
        analyzer = Psi4OutputAnalyzer(file)
        results = analyzer.analyze()
        if results:
            analyzer.print_summary()
        else:
            print(f"无法分析文件: {file}")

if __name__ == "__main__":
    main()
