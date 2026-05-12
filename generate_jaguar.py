#!/usr/bin/env python
import os
import sys
import subprocess

def print_help():
    help_text = """
================================================================================
===  Jaguar 输入文件自动生成脚本 (r2SCAN-3c + COSMO 水溶剂)  ===
================================================================================
功能：
  1. 从 SDF 文件生成 Jaguar 输入文件 (.in)
  2. 自动写入计算参数：r2SCAN-3c / isolv=7(COSMO) / nogas=2(仅溶液相)
  3. 自动生成可执行的提交脚本 (.sh)

用法：
  python generate_jaguar.py 文件名.sdf
  python generate_jaguar.py -h | --help

示例：
  python generate_jaguar.py test.sdf
================================================================================
"""
    print(help_text)
    sys.exit(0)

def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print_help()

    sdf_full = sys.argv[1]
    name = os.path.splitext(sdf_full)[0]

    in_file = f"jag_{name}_spe_r2SCAN-3c.in"
    sh_file = f"jag_{name}_spe_r2SCAN-3c.sh"

    SCHRODINGER = os.getenv("SCHRODINGER")
    if not SCHRODINGER:
        print("错误：未找到环境变量 $SCHRODINGER！")
        sys.exit(1)

    obabel = os.path.join(SCHRODINGER, "utilities", "obabel")

    print(f"==> 生成输入文件: {in_file}")
    subprocess.run([obabel, "-isdf", sdf_full, "-ojin", "-O", in_file], check=True)

    with open(in_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    new_lines = []
    skip_mode = True
    for line in lines:
        if skip_mode:
            if "&zmat" in line:
                new_lines.append("&gen\n")
                new_lines.append("isolv=7\n")
                new_lines.append("pcm_model=cosmo\n")
                new_lines.append("dftname=r2SCAN-3c\n")
                new_lines.append("nogas=2\n")
                new_lines.append("&\n")
                new_lines.append(f"entry_name: {name}\n")
                new_lines.append("&zmat\n")
                skip_mode = False
            continue
        new_lines.append(line)

    with open(in_file, "w", encoding="utf-8") as f:
        f.writelines(new_lines)

    # ====================== ✅ 输出执行脚本 ======================
    # 硬编码使用32核心计算:
    # -max_threads 32
    # 按需要更改
    sh_content = f"${{SCHRODINGER}}/jaguar run -jobname=jag_{name}_spe_r2SCAN-3c {in_file} -HOST localhost -PARALLEL 1 -max_threads 32 -TMPLAUNCHDIR\n"

    with open(sh_file, "w", encoding="utf-8") as f:
        f.write(sh_content)

    os.chmod(sh_file, 0o755)

    print("==> 完成！")
    print(f"==> 运行：bash {sh_file}")

if __name__ == "__main__":
    main()
