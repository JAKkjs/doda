#!/bin/bash

# ==============================================================================
# ESTIMATE 免疫浸润评分分析 Bash
# 调用方式：bash run_estimate.sh <用户ID> <矩阵文件路径> [分组文件路径]
# ==============================================================================

# 1. 接收参数
USER_ID=${1}
MATRIX_FILE=${2}
GROUP_FILE=${3} # 可选参数

# 2. 配置基础路径
BASE_DIR=$(cd "$(dirname "$0")";pwd)
# 结果存放目录，与差异分析对齐
RESULT_DIR="${BASE_DIR}/results/${USER_ID}/estimate"
# 确保 R 脚本路径正确
R_SCRIPT="${BASE_DIR}/Rscripts/estimate_analysis.R"
LOG_FILE="${RESULT_DIR}/estimate_analysis.log"

RSCRIPT_BIN=$(conda run -n r_doda which Rscript)

# 3. 执行前检查
echo "[$(date)] 开始任务: 用户 ${USER_ID} (ESTIMATE)"

# 检查输入矩阵文件是否存在
if [ ! -f "$MATRIX_FILE" ]; then
    echo "错误: 找不到表达矩阵文件 $MATRIX_FILE" >&2
    exit 1
fi

# 创建结果目录
mkdir -p "$RESULT_DIR"

# 4. 构建并执行 R 命令
echo "开始免疫评分分析..."

if [ -z "$GROUP_FILE" ]; then
    # 模式 B: 自动识别
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
    echo "SUCCESS:"
    echo "RESULT_CSV:${RESULT_DIR}/estimate_result.csv"
    exit 0
else
    echo "FAILED, Check: $LOG_FILE" >&2
    exit 1
fi