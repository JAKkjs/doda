#!/bin/bash

# ==============================================================================
# CIBERSORT 免疫细胞浸润分析 Bash
# 调用方式：bash run_cibersort.sh <用户ID> <矩阵文件路径> [分组文件路径]
# ==============================================================================

# 1. 接收参数
USER_ID=${1}
MATRIX_FILE=${2}
GROUP_FILE=${3}   # 可选

# 2. 配置基础路径（与 run_deg / run_estimate 一致）
BASE_DIR=$(cd "$(dirname "$0")"; pwd)
RESULT_DIR="${BASE_DIR}/results/${USER_ID}/cibersort"
R_SCRIPT="${BASE_DIR}/Rscripts/cibersort_analysis.R"
LOG_FILE="${RESULT_DIR}/cibersort_analysis.log"

RSCRIPT_BIN=$(conda run -n r_doda which Rscript)

# 3. 执行前检查
echo "[$(date)] 开始任务: 用户 ${USER_ID} (CIBERSORT)"

if [ ! -f "$MATRIX_FILE" ]; then
    echo "错误: 找不到表达矩阵文件 $MATRIX_FILE" >&2
    exit 1
fi

mkdir -p "$RESULT_DIR"

# 4. 执行 R：传入 矩阵 输出目录 [分组] 项目根目录（用于找 LM22.txt）
echo "开始 CIBERSORT 免疫浸润分析..."

# 参数顺序：矩阵 输出目录 [分组文件] 项目根目录(供 R 找 LM22.txt)
if [ -z "$GROUP_FILE" ] || [ ! -f "$GROUP_FILE" ]; then
    $RSCRIPT_BIN "$R_SCRIPT" "$MATRIX_FILE" "$RESULT_DIR" "" "$BASE_DIR" > "$LOG_FILE" 2>&1
else
    $RSCRIPT_BIN "$R_SCRIPT" "$MATRIX_FILE" "$RESULT_DIR" "$GROUP_FILE" "$BASE_DIR" > "$LOG_FILE" 2>&1
fi

# 5. 检查退出码
if [ $? -eq 0 ]; then
    echo "SUCCESS:"
    echo "RESULT_CSV:${RESULT_DIR}/cibersort_result.csv"
    exit 0
else
    echo "FAILED, Check: $LOG_FILE" >&2
    exit 1
fi
