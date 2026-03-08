# ====================================================================
# ESTIMATE 分析
# ====================================================================

args <- commandArgs(trailingOnly = TRUE)
input_matrix_file <- args[1]
output_dir        <- args[2]
group_file        <- if(length(args) >= 3) args[3] else NA

suppressMessages({
  library(estimate)
  library(limma)
  library(ggplot2)
  library(ggpubr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

# --- 1. ID 转换与矩阵读取 (修复重复行名报错版) ---
message(">>> [1/5] 读取矩阵并处理 ID 转换...")

# 【修复点 1】：先不设 row.names，避免因重复 ID 导致读取失败
raw_data <- read.csv(input_matrix_file, check.names = FALSE)

# 确保至少有两列（ID列 + 数据列）
if (ncol(raw_data) < 2) stop("错误：表达矩阵格式不正确，列数不足。")

# 提取原始 ID（第一列）
gene_names_raw <- as.character(raw_data[, 1])
data_mat_raw <- as.matrix(raw_data[, -1])

# 检查基因名类型
if (all(grepl("^[0-9]+$", gene_names_raw[1:min(10, length(gene_names_raw))]))) {
    stop("【数据格式错误】：检测到纯数字基因名。请确保行名为 Ensembl ID 或 Gene Symbol。")
}

# --- 核心逻辑：区分 Ensembl 还是 Symbol 并处理重复 ---

if (all(grepl("^ENSG[0-9]+", gene_names_raw[1:min(5, length(gene_names_raw))]))) {
    message(">>> 检测到 Ensembl ID，开始转换...")
    ensembl_ids <- gsub("\\..*", "", gene_names_raw)
    
    # 转换 ID
    gene_map <- bitr(ensembl_ids, fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)
    
    # 构建带 Symbol 的临时框进行合并
    tmp_df <- data.frame(ENSEMBL = ensembl_ids, data_mat_raw, check.names = FALSE)
    merged_data <- merge(gene_map, tmp_df, by = "ENSEMBL")
    
    # 提取数值矩阵（SYMBOL 在第 2 列，数值从第 3 列开始）
    final_mat <- as.matrix(merged_data[, 3:ncol(merged_data)])
    
    # 【修复点 2】：转换后可能多个 Ensembl 对应一个 Symbol，必须再次合并
    message(">>> 正在合并重复的 Gene Symbol...")
    data_final <- avereps(final_mat, ID = merged_data$SYMBOL)
    
} else {
    message(">>> 检测到 Gene Symbol，处理重复项...")
    # 如果已经是 Symbol，直接处理原始 ID 中的重复
    data_final <- avereps(data_mat_raw, ID = gene_names_raw)
}

# 确保数值锁定
storage.mode(data_final) <- "numeric"

# 写入临时文件供 ESTIMATE 包读取
temp_txt <- file.path(output_dir, "temp_estimate_input.txt")
write.table(data.frame(GeneSymbol=rownames(data_final), data_final), 
            file=temp_txt, quote=F, sep="\t", row.names=F)

message(">>> 矩阵处理完成，剩余有效基因数: ", nrow(data_final))

# --- 2. 运行 ESTIMATE ---
message(">>> [2/5] 运行 ESTIMATE 算法...")
common_gct <- file.path(output_dir, "common_genes.gct")
score_gct  <- file.path(output_dir, "estimate_scores.gct")

filterCommonGenes(input.f = temp_txt, output.f = common_gct, id = "GeneSymbol")
estimateScore(input.ds = common_gct, output.ds = score_gct, platform = "affymetrix")

scores_raw <- read.table(score_gct, skip = 2, header = 1)
scores <- t(scores_raw[, 3:ncol(scores_raw)])
colnames(scores) <- c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity")
rownames(scores) <- colnames(data_final) # 确保样本名完全一致
write.csv(scores, file.path(output_dir, "estimate_result.csv"))

# --- 3. 识别样本分组 (强化匹配逻辑) ---
message(">>> [3/5] 识别样本分组...")
sample_names <- rownames(scores)

if (!is.na(group_file) && file.exists(group_file)) {
    group_df <- read.csv(group_file, header = TRUE, check.names = FALSE)
    
    # 模糊匹配优化：统一转为大写并把点/横杠全部去掉再比对
    clean_matrix_names <- toupper(gsub("[._-]", "", sample_names))
    clean_group_names  <- toupper(gsub("[._-]", "", as.character(group_df[,1])))
    
    match_idx <- match(clean_matrix_names, clean_group_names)
    
    if (any(is.na(match_idx))) {
        message("精准匹配失败，尝试前缀匹配...")
        # 针对 TCGA 这种 ID 比较长的情况，取前 12 位匹配
        match_idx <- match(substr(clean_matrix_names, 1, 12), substr(clean_group_names, 1, 12))
    }

    if (any(is.na(match_idx))) stop("错误：分组文件无法匹配矩阵样本名！请检查样本 ID 命名。")
    
    final_group <- factor(group_df[match_idx, 2])
    ref_level <- grep("normal|control|healthy|0", levels(final_group), ignore.case = TRUE, value = TRUE)
    if (length(ref_level) > 0) final_group <- relevel(final_group, ref = ref_level[1])
} else {
    # 自动识别逻辑 (保持之前的增强版)
    clean_names <- gsub("\\.", "-", sample_names)
    is_tcga <- any(grepl("-(0[1-9]|1[0-9])[A-Z]?$", clean_names))
    if (is_tcga) {
        codes <- as.numeric(substr(clean_names, nchar(clean_names)-2, nchar(clean_names)-1))
        group_vec <- ifelse(!is.na(codes) & codes < 10, "Tumor", "Normal")
    } else {
        group_vec <- rep("Unknown", length(sample_names))
        lower_n <- swallow_names <- tolower(clean_names)
        group_vec[grepl("normal|control|healthy|adj", lower_n)] <- "Normal"
        group_vec[grepl("tumor|cancer|treat|case", lower_n)] <- "Tumor"
    }
    final_group <- factor(group_vec)
    if ("Normal" %in% levels(final_group)) final_group <- relevel(final_group, ref = "Normal")
}

# --- 4. 绘图可视化 ---
message(">>> [4/5] 生成图表...")
plot_data <- data.frame(
    Value = as.vector(scores),
    ScoreType = rep(colnames(scores), each = length(sample_names)),
    Group = rep(final_group, 4)
)

p <- ggboxplot(plot_data, x = "Group", y = "Value", fill = "Group", 
               palette = "npg", facet.by = "ScoreType", scales = "free") +
     geom_jitter(shape = 21, size = 0.8, width = 0.2, alpha = 0.3) +
     stat_compare_means(method = "t.test", label = "p.format") +
     theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")

ggsave(file.path(output_dir, "estimate_plot.png"), p, width = 10, height = 8, bg = "white")

# --- 5. 清理 ---
message(">>> [5/5] 清理临时文件...")
file.remove(temp_txt, common_gct, score_gct)
message(">>> 分析任务已完成！")