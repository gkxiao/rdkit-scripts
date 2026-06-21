#!/bin/bash
# 脚本名称：xtbopt_batch.sh
# 用法：./xtbopt_batch.sh [ -h | --help ] <formal_charge>

# 显示帮助信息的函数
show_help() {
    cat << EOF
用法：$0 [选项] <formal_charge>

对当前目录下所有 CONF_*.xyz 文件进行 xTB 几何优化，
并将优化后的结构合并到 ensemble.xtbopt.xyz。

选项：
  -h, --help     显示此帮助信息并退出

参数：
  formal_charge  形式电荷，必须为整数（例如 0, -1, 2）

示例：
  $0 0           # 电荷为 0，优化所有 CONF_*.xyz
  $0 -1          # 电荷为 -1

注意：
  - 如果 ensemble.xtbopt.xyz 已存在，脚本会报错退出，请手动处理。
  - 需要 xTB 可执行文件在 PATH 中。
EOF
}

# 解析命令行参数
if [ $# -eq 0 ]; then
    echo "错误：缺少形式电荷参数"
    echo "使用 -h 或 --help 查看帮助"
    exit 1
fi

if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
    show_help
    exit 0
fi

# 此时第一个参数应为电荷
charge="$1"

# 检查形式电荷是否为整数（简单判断，允许正负号）
if ! [[ "$charge" =~ ^-?[0-9]+$ ]]; then
    echo "错误：形式电荷必须是整数"
    echo "使用 -h 或 --help 查看帮助"
    exit 1
fi

# 检查是否存在 CONF_*.xyz 文件
if ! ls CONF_*.xyz 1>/dev/null 2>&1; then
    echo "错误：当前目录下没有找到 CONF_*.xyz 文件"
    exit 1
fi

# 统计文件数量
n=$(ls CONF_*.xyz | wc -l)

# 检查 ensemble.xtbopt.xyz 是否已存在
if [ -f ensemble.xtbopt.xyz ]; then
    echo "错误：ensemble.xtbopt.xyz 已存在，请手动删除或重命名后重新运行"
    exit 1
fi

# 循环优化每个文件
for i in $(seq 1 "$n"); do
    echo "正在优化 CONF_${i}.xyz，电荷 = ${charge}"
    xtb "CONF_${i}.xyz" -c "$charge" -u 0 --opt tight --alpb water --namespace "CONF_${i}"

    # 检查优化后的文件是否生成，并追加到 ensemble.xtbopt.xyz
    if [ -f "CONF_${i}.xtbopt.xyz" ]; then
        cat "CONF_${i}.xtbopt.xyz" >> ensemble.xtbopt.xyz
    else
        echo "警告：CONF_${i}.xtbopt.xyz 未生成，跳过该文件"
    fi
done

echo "全部优化完成！结果已合并至 ensemble.xtbopt.xyz"
