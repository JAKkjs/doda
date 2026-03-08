# ====================================================================
#  GSVA 基因集变异分析
#  调用：Rscript GSVA.R "matrix.csv" "./result" [group_file] [base_dir] [gmt_file]
#  输入：表达矩阵(行=基因 Symbol, 列=样本)；GMT 默认 Rscripts/c2.cp.kegg.*.gmt
# ====================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("错误：至少需要 表达矩阵路径、输出目录。")

input_matrix_file <- args[1]
output_dir        <- args[2]
group_file        <- if (length(args) >= 3 && nzchar(trimws(args[3]))) args[3] else NA
base_dir          <- if (length(args) >= 4 && nzchar(trimws(args[4]))) args[4] else getwd()
gmt_file          <- if (length(args) >= 5 && nzchar(trimws(args[5]))) args[5] else NA

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
orig_wd <- getwd()
on.exit(setwd(orig_wd), add = TRUE)

suppressMessages({
  library(limma)
  library(GSEABase)
  library(GSVA)
  library(ggplot2)
})
if (requireNamespace("clusterProfiler", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
}

# --------------------------------------------------------------------
# 1. GMT 文件路径
# --------------------------------------------------------------------
if (is.na(gmt_file) || !file.exists(gmt_file)) {
  cand <- list.files(file.path(base_dir, "Rscripts"), pattern = "\\.gmt$", full.names = TRUE)
  if (length(cand) > 0) gmt_file <- cand[1] else gmt_file <- file.path(base_dir, "Rscripts", "c2.cp.kegg.v2023.1.Hs.symbols.gmt")
}
if (!file.exists(gmt_file)) {
  stop("错误：未找到 GMT 文件。请将基因集 GMT 文件置于 Rscripts/ 目录下（如 MSigDB c2.cp.kegg 或 hallmark），或通过第 5 个参数指定路径。下载示例: https://www.gsea-msigdb.org/gsea/msigdb/")
}

message(">>> 使用 GMT 文件: ", gmt_file)
geneSets <- getGmt(gmt_file, geneIdType = SymbolIdentifier())

# --------------------------------------------------------------------
# 2. 读取表达矩阵（修正版：处理重复行名与类型转换）
# --------------------------------------------------------------------
message(">>> [1/4] 读取表达矩阵...")
# 先不设置 row.names，避免读取时因重复 ID 报错
raw <- read.csv(input_matrix_file, check.names = FALSE) 

if (ncol(raw) < 2) stop("错误：表达矩阵至少需要两列（基因名 + 至少一个样本）。")

# 提取基因名（第一列）和数据（其他列）
gene_names_raw <- as.character(raw[, 1])
data_raw <- as.matrix(raw[, -1, drop = FALSE])

# 确保数据列是数值型（非常重要！防止因读取错误变成 character）
if (!is.numeric(data_raw)) {
  suppressWarnings(storage.mode(data_raw) <- "numeric")
}

# 处理原始 ID 的重复值（如重复的 ENSG ID）
if (any(duplicated(gene_names_raw))) {
  message(">>> 检测到原始 ID 重复，已按平均合并。")
  rownames(data_raw) <- gene_names_raw
  data_raw <- limma::avereps(data_raw)
  gene_names_raw <- rownames(data_raw)
} else {
  rownames(data_raw) <- gene_names_raw
}

# Ensembl 转 Symbol 逻辑
if (all(grepl("^ENSG[0-9]+", gene_names_raw[1:min(5, length(gene_names_raw))]))) {
  message(">>> 检测到 Ensembl ID，转换为 Gene Symbol...")
  
  ensembl_ids <- gsub("\\..*", "", gene_names_raw)
  gene_map <- bitr(ensembl_ids, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  
  # 构建带 Symbol 的临时数据框
  tmp_df <- data.frame(ENSEMBL = ensembl_ids, data_raw, check.names = FALSE)
  merged_data <- merge(gene_map, tmp_df, by = "ENSEMBL")
  
  # 提取数值矩阵：第 1 列是 ENSEMBL, 第 2 列是 SYMBOL，数值从第 3 列开始
  final_mat <- as.matrix(merged_data[, 3:ncol(merged_data)])
  rownames(final_mat) <- merged_data$SYMBOL
  
  # 对转换后的 Symbol 再次合并重复值（多个 ENSG 对应一个 Symbol）
  data <- limma::avereps(final_mat)
  message(">>> 转换完成，剩余基因数: ", nrow(data))
} else {
  data <- data_raw
}

# 最终检查：确保没有非数值内容
if (any(is.na(data))) {
  message(">>> 警告：矩阵中包含 NA，已替换为 0。")
  data[is.na(data)] <- 0
}
storage.mode(data) <- "numeric" # 强制最后锁定类型

# --------------------------------------------------------------------
# 3. 运行 GSVA
# --------------------------------------------------------------------
message(">>> [2/4] 运行 GSVA...")
param <- GSVA::gsvaParam(exprData = data, geneSets = geneSets)
gsvaResult <- gsva(param, verbose = FALSE)
# gsvaResult: 行=通路, 列=样本
write.csv(t(gsvaResult), file.path(output_dir, "gsva_scores.csv"))
message(">>> GSVA 得分已保存: gsva_scores.csv")

# --------------------------------------------------------------------
# 4. 分组
# --------------------------------------------------------------------
message(">>> [3/4] 识别分组...")
sample_names <- colnames(gsvaResult)
final_group <- NULL

if (!is.na(group_file) && file.exists(group_file)) {
  group_df <- read.csv(group_file, header = TRUE, check.names = FALSE)
  clean_matrix <- toupper(gsub("[._-]", "", sample_names))
  clean_group  <- toupper(gsub("[._-]", "", as.character(group_df[, 1])))
  match_idx <- match(clean_matrix, clean_group)
  if (any(is.na(match_idx))) match_idx <- match(substr(clean_matrix, 1, 12), substr(clean_group, 1, 12))
  if (any(is.na(match_idx))) message(">>> 分组文件与样本未完全匹配，跳过分组差异与图。") else {
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

# --------------------------------------------------------------------
# 5. 有分组时：limma 通路差异 + 发散条形图
# --------------------------------------------------------------------
message(">>> [4/4] 通路差异与绘图...")
if (!is.null(final_group) && length(levels(final_group)) >= 2) {
  design <- model.matrix(~ 0 + final_group)
  colnames(design) <- levels(final_group)
  fit <- lmFit(gsvaResult, design)
  contrast_formula <- paste(colnames(design)[2], "-", colnames(design)[1])
  contrast.matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  allDiff <- topTable(fit2, coef = 1, n = Inf, sort.by = "P")
  write.csv(allDiff, file.path(output_dir, "gsva_diff_all.csv"))
  diff_sig <- allDiff[abs(allDiff$logFC) > 0.1 & allDiff$adj.P.Val < 0.05, ]
  write.csv(diff_sig, file.path(output_dir, "gsva_diff_sig.csv"))
  message(">>> 通路差异结果已保存: gsva_diff_all.csv, gsva_diff_sig.csv")

  # 发散条形图（上下调前 N 条通路）
  p_cutoff <- 0.05
  degs <- allDiff[order(allDiff$logFC, decreasing = TRUE), ]
  n_show <- 20
  up_n <- min(n_show, sum(degs$logFC > 0))
  down_n <- min(n_show, sum(degs$logFC < 0))
  idx_up <- head(which(degs$logFC > 0), up_n)
  idx_down <- tail(which(degs$logFC < 0), down_n)
  Diff <- degs[c(idx_up, idx_down), , drop = FALSE]
  dat_plot <- data.frame(
    id = rownames(Diff),
    p = Diff$P.Value,
    lgfc = Diff$logFC,
    stringsAsFactors = FALSE
  )
  dat_plot$group <- ifelse(dat_plot$lgfc > 0, 1, -1)
  dat_plot$lg_p <- -log10(dat_plot$p + 1e-10) * dat_plot$group
  dat_plot$threshold <- factor(
    ifelse(dat_plot$p <= p_cutoff, ifelse(dat_plot$lgfc > 0, "Up", "Down"), "Not"),
    levels = c("Up", "Down", "Not")
  )
  dat_plot <- dat_plot[order(dat_plot$lg_p), ]
  dat_plot$id <- factor(dat_plot$id, levels = dat_plot$id)

  p <- ggplot(dat_plot, aes(x = id, y = lg_p, fill = threshold)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("Up" = "#36638a", "Not" = "#cccccc", "Down" = "#7bcd7b")) +
    geom_hline(yintercept = c(-log10(p_cutoff), log10(p_cutoff)), color = "gray50", size = 0.5, linetype = "dashed") +
    xlab("") + ylab("-log10(P) * sign(logFC)") +
    theme_bw() + theme(panel.grid = element_blank(), legend.position = "none",
                      axis.text.y = element_text(size = 7))
  ggsave(file.path(output_dir, "gsva_barplot.png"), p, width = 8, height = max(5, nrow(dat_plot) * 0.2), dpi = 150, bg = "white")
  message(">>> 发散条形图已保存: gsva_barplot.png")
} else {
  # 无分组：热图（通路较多时只画前 50 条或聚类后取代表）
  mat <- as.matrix(gsvaResult)
  n_path <- min(50, nrow(mat))
  mat_show <- mat[order(apply(mat, 1, sd), decreasing = TRUE)[1:n_path], , drop = FALSE]
  pdf(file.path(output_dir, "gsva_heatmap.pdf"), width = max(8, ncol(mat) * 0.12), height = max(5, n_path * 0.08))
  heatmap(mat_show, scale = "row", col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
          margins = c(6, 8), main = "GSVA pathway scores")
  dev.off()
  message(">>> 热图已保存: gsva_heatmap.pdf")
}

message(">>> GSVA 分析完成。")
