# ====================================================================
#  CIBERSORT 免疫细胞浸润分析
#  调用方式：
#    方式1 (指定分组): Rscript cibersort_analysis.R "matrix.csv" "./result" "group_info.csv"
#    方式2 (仅浸润):   Rscript cibersort_analysis.R "matrix.csv" "./result"
#  输入矩阵：行=基因(Symbol 或 Ensembl)，列=样本；建议 FPKM/TPM 等表达量。
# ====================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("错误：参数不足！至少需要：输入表达矩阵路径、输出目录。")
}

input_matrix_file <- args[1]
output_dir        <- args[2]
group_file        <- if (length(args) >= 3 && nzchar(trimws(args[3]))) args[3] else NA
base_dir          <- if (length(args) >= 4 && nzchar(trimws(args[4]))) args[4] else getwd()

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
orig_wd <- getwd()
on.exit(setwd(orig_wd), add = TRUE)

suppressMessages({
  library(limma)
  library(reshape2)
  library(ggpubr)
  library(ggplot2)
})

if (!requireNamespace("CIBERSORT", quietly = TRUE)) {
  stop("请安装 CIBERSORT 包。若使用 IOBR，请确保已安装并加载 CIBERSORT 兼容接口。")
}
suppressMessages(library(CIBERSORT))

# --------------------------------------------------------------------
# 1. 加载 LM22 签名矩阵（base_dir 由 Bash 传入，即项目根目录）
# --------------------------------------------------------------------
lm22_path <- file.path(base_dir, "Rscripts", "LM22.txt")

if (file.exists(lm22_path)) {
  message(">>> 使用脚本同目录 LM22.txt")
  lm22_raw <- read.delim(lm22_path, check.names = FALSE)
  gene_col <- which(tolower(names(lm22_raw)) %in% c("gene symbol", "gene_symbol", "genesymbol", "gene"))[1]
  if (is.na(gene_col)) gene_col <- 1
  LM22 <- as.matrix(lm22_raw[, -gene_col])
  rownames(LM22) <- as.character(lm22_raw[[gene_col]])
} else {
  message(">>> 尝试从 CIBERSORT 包加载 data(LM22)")
  if (exists("LM22", envir = .GlobalEnv)) {
    LM22 <- get("LM22", envir = .GlobalEnv)
  } else {
    tryCatch({
      data(LM22, package = "CIBERSORT")
      LM22 <- as.matrix(LM22)
    }, error = function(e) {
      stop("未找到 LM22.txt 且 CIBERSORT 包无 LM22 数据，请将 LM22.txt 放在 Rscripts/ 目录下。")
    })
  }
}

# --------------------------------------------------------------------
# 2. 读取表达矩阵并统一为 Gene Symbol（允许重复行名，读入后合并）
# --------------------------------------------------------------------
message(">>> [1/4] 读取表达矩阵...")
raw <- read.csv(input_matrix_file, check.names = FALSE)
if (ncol(raw) < 2) stop("错误：表达矩阵至少需要两列（基因名 + 至少一个样本）。")
gene_names <- as.character(raw[, 1])
data_raw <- as.matrix(raw[, -1, drop = FALSE])
rownames(data_raw) <- gene_names
# 若有重复基因名，按行取平均
if (any(duplicated(gene_names))) {
  message(">>> 检测到重复基因名，已按平均合并。")
  data_raw <- avereps(data_raw)
}
gene_names <- rownames(data_raw)

if (all(grepl("^ENSG[0-9]+", gene_names))) {
  message(">>> 检测到 Ensembl ID，转换为 Gene Symbol...")
  suppressMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
  })
  ensembl_ids <- gsub("\\..*", "", gene_names)
  gene_map <- bitr(ensembl_ids, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  data_with_symbol <- merge(gene_map, cbind(ENSEMBL = ensembl_ids, data_raw), by = "ENSEMBL")
  data_expr <- as.matrix(data_with_symbol[, 3:ncol(data_with_symbol)])
  rownames(data_expr) <- data_with_symbol$SYMBOL
  data_expr <- avereps(data_expr)
} else {
  data_expr <- as.matrix(data_raw)
  if (any(duplicated(rownames(data_expr)))) {
    data_expr <- avereps(data_expr)
  }
}

# 与 LM22 取基因交集
common_genes <- intersect(rownames(LM22), rownames(data_expr))
if (length(common_genes) < 100) {
  stop("错误：与 LM22 共有基因数过少 (", length(common_genes), ")，请检查表达矩阵行名是否为 Gene Symbol 或 Ensembl。")
}
message(">>> 与 LM22 共有基因数: ", length(common_genes))

LM22_sub <- LM22[common_genes, , drop = FALSE]
mix_sub <- data_expr[common_genes, , drop = FALSE]

# --------------------------------------------------------------------
# 3. 运行 CIBERSORT
# --------------------------------------------------------------------
message(">>> [2/4] 运行 CIBERSORT...")
setwd(output_dir)
results <- tryCatch({
  cibersort(LM22_sub, mix_sub, perm = 100, QN = FALSE)
}, error = function(e) {
  stop("CIBERSORT 运行失败: ", e$message)
})

# 结果可能带 P-value、Correlation 等列，前 n 列为 22 种细胞比例
cell_cols <- intersect(colnames(results), colnames(LM22))
if (length(cell_cols) == 0) cell_cols <- colnames(results)[1:min(22, ncol(results))]
result_core <- results[, cell_cols, drop = FALSE]

out_csv <- file.path(output_dir, "cibersort_result.csv")
write.csv(result_core, out_csv)
message(">>> 结果已保存: ", out_csv)

# 若有 P-value 等列也一并保存完整表
write.csv(results, file.path(output_dir, "cibersort_result_full.csv"))

# --------------------------------------------------------------------
# 4. 分组与可视化（若有分组则按组比较细胞比例）
# --------------------------------------------------------------------
message(">>> [3/4] 识别分组...")
sample_names <- rownames(result_core)
final_group <- NULL

if (!is.na(group_file) && file.exists(group_file)) {
  group_df <- read.csv(group_file, header = TRUE, check.names = FALSE)
  clean_matrix <- toupper(gsub("[._-]", "", sample_names))
  clean_group  <- toupper(gsub("[._-]", "", as.character(group_df[, 1])))
  match_idx <- match(clean_matrix, clean_group)
  if (any(is.na(match_idx))) {
    match_idx <- match(substr(clean_matrix, 1, 12), substr(clean_group, 1, 12))
  }
  if (any(is.na(match_idx))) {
    message(">>> 分组文件与样本未完全匹配，跳过分组图。")
  } else {
    final_group <- factor(group_df[match_idx, 2])
    ref_level <- grep("normal|control|healthy|0", levels(final_group), ignore.case = TRUE, value = TRUE)
    if (length(ref_level) > 0) final_group <- relevel(final_group, ref = ref_level[1])
  }
} else {
  clean_names <- gsub("\\.", "-", sample_names)
  is_tcga <- any(grepl("-(0[1-9]|1[0-9])[A-Z]?$", clean_names))
  if (is_tcga) {
    codes <- as.numeric(substr(clean_names, nchar(clean_names) - 2, nchar(clean_names) - 1))
    group_vec <- ifelse(!is.na(codes) & codes < 10, "Tumor", "Normal")
  } else {
    group_vec <- rep("Unknown", length(sample_names))
    lower_n <- tolower(clean_names)
    group_vec[grepl("normal|control|healthy|adj", lower_n)] <- "Normal"
    group_vec[grepl("tumor|cancer|treat|case", lower_n)] <- "Tumor"
  }
  if (!all(group_vec == "Unknown")) {
    final_group <- factor(group_vec)
    if ("Normal" %in% levels(final_group)) final_group <- relevel(final_group, ref = "Normal")
  }
}

message(">>> [4/4] 生成图表...")
if (!is.null(final_group) && length(levels(final_group)) >= 2) {
  sample_results <- data.frame(result_core, Group = final_group, check.names = FALSE)
  data_long <- melt(sample_results, id.vars = "Group", variable.name = "CellType", value.name = "Proportion")
  # 按细胞类型算 Wilcoxon p 值，避免 ggpubr::stat_compare_means 在 facet 下报错
  pval_df <- do.call(rbind, lapply(unique(data_long$CellType), function(ct) {
    d <- data_long[data_long$CellType == ct, ]
    pv <- tryCatch(wilcox.test(Proportion ~ Group, data = d)$p.value, error = function(e) NA)
    data.frame(CellType = ct, p = pv, y = max(d$Proportion, na.rm = TRUE) * 1.05)
  }))
  pval_df$label <- sprintf("p = %.2g", pval_df$p)
  p <- ggplot(data_long, aes(x = Group, y = Proportion, fill = Group)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_jitter(shape = 21, size = 0.6, width = 0.2, alpha = 0.4) +
    geom_text(data = pval_df, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 2.5) +
    facet_wrap(~ CellType, scales = "free_y", ncol = 4) +
    scale_fill_brewer(palette = "Set2") +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none",
                      axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(output_dir, "cibersort_plot.png"), p, width = 14, height = 10, dpi = 150, bg = "white")
  message(">>> 分组箱线图已保存: ", file.path(output_dir, "cibersort_plot.png"))
} else {
  # 无分组：画各样本 22 细胞类型堆叠条形图或热图
  mat_plot <- as.matrix(result_core)
  pdf(file.path(output_dir, "cibersort_heatmap.pdf"), width = max(8, ncol(mat_plot) * 0.15), height = 6)
  heatmap(t(mat_plot), scale = "none", col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
          margins = c(8, 10), main = "CIBERSORT immune cell proportions")
  dev.off()
  message(">>> 热图已保存: ", file.path(output_dir, "cibersort_heatmap.pdf"))
}

message(">>> CIBERSORT 分析完成。")
