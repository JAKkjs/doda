# ====================================================================
#  智能免疫检查点差异表达分析 (最终修复版)
# ====================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("错误：至少需要 表达矩阵路径、输出目录。")

input_matrix_file <- args[1]
output_dir        <- args[2]
group_file        <- if (length(args) >= 3 && nzchar(trimws(args[3]))) args[3] else NA
base_dir          <- if (length(args) >= 4 && nzchar(trimws(args[4]))) args[4] else getwd()

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

suppressMessages({
  library(limma)
  library(reshape2)
  library(ggplot2)
})

# --------------------------------------------------------------------
# 1. 资源加载与矩阵读取 (修复重复行名报错)
# --------------------------------------------------------------------
icg_path <- file.path(base_dir, "Rscripts", "ICG_genes.txt")
if (!file.exists(icg_path)) stop("错误：未找到 ICG_genes.txt。")
icg_names <- unique(as.character(read.table(icg_path, header = TRUE)[,1]))

message(">>> [1/4] 读取并预处理矩阵...")
# 先不设 row.names，避免重复基因名导致坠机
raw_data <- read.csv(input_matrix_file, check.names = FALSE, stringsAsFactors = FALSE)
gene_names_vec <- as.character(raw_data[, 1])
numeric_mat <- as.matrix(raw_data[, -1, drop = FALSE])

# 在这里处理重复基因名
if (any(duplicated(gene_names_vec))) {
  message(">>> 检测到重复基因名，正在合并...")
  rownames(numeric_mat) <- gene_names_vec
  numeric_mat <- avereps(numeric_mat)
} else {
  rownames(numeric_mat) <- gene_names_vec
}

# 自动处理 Ensembl 转 Symbol
if (all(grepl("^ENSG", rownames(numeric_mat)[1:min(5, nrow(numeric_mat))]))) {
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    suppressMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
    gene_map <- bitr(gsub("\\..*", "", rownames(numeric_mat)), fromType="ENSEMBL", toType="SYMBOL", OrgDb=org.Hs.eg.db)
    # 匹配并转换
    numeric_mat <- numeric_mat[match(gene_map$ENSEMBL, gsub("\\..*", "", rownames(numeric_mat))), ]
    rownames(numeric_mat) <- gene_map$SYMBOL
    numeric_mat <- avereps(numeric_mat)
  }
}

# 提取交集基因
icg_in_matrix <- intersect(icg_names, rownames(numeric_mat))
if (length(icg_in_matrix) == 0) stop("错误：表达矩阵中未检测到任何免疫检查点基因。")
icg_exp <- numeric_mat[icg_in_matrix, , drop = FALSE]
message(">>> 检测到 ICG 基因数: ", length(icg_in_matrix))

# --------------------------------------------------------------------
# 2. 智能分组识别
# --------------------------------------------------------------------
message(">>> [2/4] 智能识别分组...")
sample_names <- colnames(icg_exp)
final_group <- NULL

# 逻辑 A：如果有分组文件
if (!is.na(group_file) && file.exists(group_file)) {
  group_df <- read.csv(group_file, header = TRUE, check.names = FALSE)
  
  # 寻找含有 risk/group/type 的列
  target_col <- which(tolower(colnames(group_df)) %in% c("risk", "group", "type", "condition", "status"))
  group_col_idx <- if(length(target_col) > 0) target_col[1] else ncol(group_df)
  
  # 模糊匹配 ID：去掉 . - _ 并取前 12 位
  clean_mtx_ids <- toupper(gsub("[._-]", "", substr(sample_names, 1, 12)))
  clean_grp_ids <- toupper(gsub("[._-]", "", substr(as.character(group_df[,1]), 1, 12)))
  
  match_idx <- match(clean_mtx_ids, clean_grp_ids)
  valid_mask <- !is.na(match_idx)
  
  if (sum(valid_mask) >= 2) {
    icg_exp <- icg_exp[, valid_mask, drop = FALSE]
    final_group <- factor(group_df[match_idx[valid_mask], group_col_idx])
    sample_names <- colnames(icg_exp)
    message(">>> 分组文件匹配成功，样本数: ", ncol(icg_exp))
  }
}

# 逻辑 B：若无分组文件，开启关键词+TCGA语义识别
if (is.null(final_group)) {
  group_vec <- rep("Unknown", length(sample_names))
  # 语义匹配
  group_vec[grepl("tumor|cancer|case|treat|high", sample_names, ignore.case=T)] <- "Tumor"
  group_vec[grepl("normal|control|healthy|adj|low", sample_names, ignore.case=T)] <- "Normal"
  # TCGA 编码匹配
  tcga_mask <- grepl("^TCGA-", sample_names, ignore.case=T)
  if (any(tcga_mask)) {
    codes <- as.numeric(substr(sample_names[tcga_mask], 14, 15))
    group_vec[tcga_mask] <- ifelse(!is.na(codes) & codes < 10, "Tumor", "Normal")
  }
  
  if (length(unique(group_vec[group_vec != "Unknown"])) >= 2) {
    keep_idx <- group_vec != "Unknown"
    icg_exp <- icg_exp[, keep_idx, drop = FALSE]
    final_group <- factor(group_vec[keep_idx])
    # 设置对照组
    ref_lv <- grep("normal|control|low", levels(final_group), ignore.case=T, value=T)
    if(length(ref_lv) > 0) final_group <- relevel(final_group, ref = ref_lv[1])
  }
}

# --------------------------------------------------------------------
# 3. 差异分析与绘图 (健壮性增强 + 语义配色版)
# --------------------------------------------------------------------
message(">>> [3/4] 准备统计分析...")

# 转置矩阵并保存
icg_exp_t <- t(icg_exp)
write.csv(icg_exp_t, file.path(output_dir, "icg_expression.csv"))

if (!is.null(final_group) && length(levels(final_group)) >= 2) {
    message(">>> 正在执行组间差异检验...")
    
    # 构造绘图长表
    plot_df <- data.frame(icg_exp_t, Group = final_group, check.names = FALSE)
    data_long <- melt(plot_df, id.vars = "Group", variable.name = "gene", value.name = "val")
    data_long$val <- as.numeric(data_long$val)

    # --- 统计：T-test ---
    stat_res <- do.call(rbind, lapply(unique(as.character(data_long$gene)), function(g){
        d <- data_long[data_long$gene == g,]
        # 增加容错：如果一组全为0或样本太少，跳过
        pv <- tryCatch({
            if(length(unique(d$Group)) < 2) NA 
            else t.test(val ~ Group, data = d)$p.value
        }, error = function(e) NA)
        data.frame(gene = g, p.value = pv)
    }))
    
    # 排序并保存统计结果
    stat_res <- stat_res[order(stat_res$p.value, na.last = TRUE), ]
    write.csv(stat_res, file.path(output_dir, "icg_stats.csv"), row.names = FALSE)

    # --- 筛选绘图基因 ---
    # 取 P 值最显著的前 12 个，如果没有显著的，就取前 12 个
    top_genes <- as.character(stat_res$gene[1:min(12, nrow(stat_res))])
    data_plot <- data_long[data_long$gene %in% top_genes, ]
    data_plot$gene <- factor(data_plot$gene, levels = rev(top_genes))

    # --- 核心：智能语义配色逻辑 ---
    get_smart_colors <- function(lvls) {
        cols <- rep("#3C8DBCFF", length(lvls)) # 默认蓝色系
        names(cols) <- lvls
        
        # 只要名字里包含这些关键词，就强制指定颜色
        # 实验组/高风险/肿瘤 -> 红色
        high_idx <- grepl("high|tumor|cancer|case|treat", lvls, ignore.case = TRUE)
        cols[high_idx] <- "#E64B35FF" 
        
        # 对照组/低风险/正常 -> 蓝色或绿色
        low_idx <- grepl("low|normal|control|healthy|adj", lvls, ignore.case = TRUE)
        cols[low_idx] <- "#4DBBD5FF"
        
        return(cols)
    }
    
    current_colors <- get_smart_colors(levels(final_group))

    # --- 绘图 ---
    message(">>> 生成 ICG 差异表达图...")
    p <- ggplot(data_plot, aes(x = gene, y = val, fill = Group)) +
        geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
        # 散点颜色固定为深灰色，背景填充随组别，这样更美观
        geom_jitter(shape = 21, size = 1.2, width = 0.15, alpha = 0.4, color = "black") + 
        coord_flip() + 
        theme_bw() +
        scale_fill_manual(values = current_colors) +
        labs(title = "Immune Checkpoint Differential Expression", 
             subtitle = paste("Comparison based on:", paste(levels(final_group), collapse = " vs ")),
             x = "", y = "Log2 (Expression + 1)") +
        theme(legend.position = "top", 
              panel.grid.minor = element_blank(),
              axis.text.y = element_text(face = "bold", size = 10),
              plot.title = element_text(hjust = 0.5, face = "bold"))

    ggsave(file.path(output_dir, "icg_plot.png"), p, width = 8, height = 7, dpi = 150, bg = "white")
    message(">>> [4/4] 分析任务已全部圆满完成！")
    
} else {
    message(">>> [!] 警告：有效分组不足（少于 2 组），仅保存表达数据，跳过差异绘图。")
}