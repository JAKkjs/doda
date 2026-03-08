#!/bin/bash

# ==============================================================================
# 免疫检查点差异表达分析 Bash
# 调用：bash run_immune_checkpoint.sh <用户ID> <矩阵文件路径> [分组文件路径]
# ==============================================================================

USER_ID=${1}
MATRIX_FILE=${2}
GROUP_FILE=${3}

BASE_DIR=$(cd "$(dirname "$0")"; pwd)
RESULT_DIR="${BASE_DIR}/results/${USER_ID}/immune_checkpoint"
R_SCRIPT="${BASE_DIR}/Rscripts/immune_checkpoint_de.R"
LOG_FILE="${RESULT_DIR}/immune_checkpoint_analysis.log"

RSCRIPT_BIN=$(conda run -n r_doda which Rscript)

echo "[$(date)] 开始任务: 用户 ${USER_ID} (免疫检查点)"

if [ ! -f "$MATRIX_FILE" ]; then
    echo "错误: 找不到表达矩阵文件 $MATRIX_FILE" >&2
    exit 1
fi

mkdir -p "$RESULT_DIR"

if [ -z "$GROUP_FILE" ] || [ ! -f "$GROUP_FILE" ]; then
    $RSCRIPT_BIN "$R_SCRIPT" "$MATRIX_FILE" "$RESULT_DIR" "" "$BASE_DIR" > "$LOG_FILE" 2>&1
else
    $RSCRIPT_BIN "$R_SCRIPT" "$MATRIX_FILE" "$RESULT_DIR" "$GROUP_FILE" "$BASE_DIR" > "$LOG_FILE" 2>&1
fi

if [ $? -eq 0 ]; then
    echo "SUCCESS:"
    echo "RESULT_CSV:${RESULT_DIR}/icg_expression.csv"
    exit 0
else
    echo "FAILED, Check: $LOG_FILE" >&2
    exit 1
fi
