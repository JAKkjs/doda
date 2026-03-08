#!/bin/bash

# ==============================================================================
# 差异化表达分析bash
# 调用方式：bash run_deg.sh <用户ID> <矩阵文件路径> [分组文件路径]
# ==============================================================================

# 1. 接收 Java 传来的参数
USER_ID=${1}
MATRIX_FILE=${2}
GROUP_FILE=${3} # 可选参数

# 2. 基础路径
BASE_DIR=$(cd "$(dirname "$0")";pwd) #这个就是当前位置
RESULT_DIR="${BASE_DIR}/results/${USER_ID}/deg"
R_SCRIPT="${BASE_DIR}/Rscripts/diff_analysis.R"
LOG_FILE="${RESULT_DIR}/deg_analysis.log"

RSCRIPT_BIN=$(conda run -n r_doda which Rscript)

# 3. 执行前检查
echo "[$(date)] 开始任务: 用户 ${USER_ID}"

# 检查输入文件是否存在
if [ ! -f "$MATRIX_FILE" ]; then
    echo "错误: 找不到表达矩阵文件 $MATRIX_FILE" >&2
    exit 1
fi

# 创建结果目录 (-p 确保父目录不存在时也能创建)
mkdir -p "$RESULT_DIR"

# 4. 构建并执行 R 命令
echo "开始差异化表达分析..."

if [ -z "$GROUP_FILE" ]; then
    # 模式 B: 自动识别
    # 将标准输出和错误都记录到 log 文件中，方便排查
    $RSCRIPT_BIN "$R_SCRIPT" "$MATRIX_FILE" "$RESULT_DIR" > "$LOG_FILE" 2>&1
else
    # 模式 A: 指定分组
    if [ ! -f "$GROUP_FILE" ]; then
        echo "错误: 找不到分组文件 $GROUP_FILE" >&2
        exit 1
    fi
    $RSCRIPT_BIN "$R_SCRIPT" "$MATRIX_FILE" "$RESULT_DIR" "$GROUP_FILE" > "$LOG_FILE" 2>&1
fi

# 5. 检查 R 运行结果
if [ $? -eq 0 ]; then
    echo "STATUS:SUCCESS"

    echo "RESULT_FILE:$(realpath ${RESULT_DIR}/diff_result.csv)"

    echo "PLOT_FILE:$(realpath ${RESULT_DIR}/volcano_plot.png)"
    exit 0
else
    echo "STATUS:FAILED"
    echo "LOG:$(realpath $LOG_FILE)" >&2
    exit 1
fi
