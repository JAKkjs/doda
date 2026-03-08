# ====================================================================
#  差异表达分析 
#  调用方式：
#    方式1 (指定分组表): Rscript script.R "counts.csv" "./result" "group_info.csv"
#    方式2 (全自动):     Rscript script.R "counts.csv" "./result"
# ====================================================================

# 1. 接收参数
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("错误：参数不足！至少需要：输入表达矩阵路径、输出目录")
}

input_matrix_file <- args[1]  # 表达矩阵
output_dir        <- args[2]  # 输出目录
group_file        <- NA       # 分组文件路径 (可选)

if (length(args) >= 3) {
  group_file <- args[3]
}

# 检查并创建输出目录
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

suppressMessages(library(edgeR))
suppressMessages(library(limma))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(ggplot2))

print(paste("正在处理表达矩阵:", input_matrix_file))

# 2. 读取表达矩阵
tryCatch({
  rnaCounts <- read.csv(input_matrix_file, row.names = 1, check.names = FALSE)
}, error = function(e) {
  stop("错误：无法读取表达矩阵 CSV。", e$message)
})

sample_names <- colnames(rnaCounts)
print(paste("检测到样本数:", length(sample_names)))

# ====================================================================
#  3. 分组处理逻辑 (双模式切换)
# ====================================================================
final_group <- NULL

# --- 模式 A: 用户上传了分组文件 ---
if (!is.na(group_file) && file.exists(group_file)) {
  print(paste("检测到分组文件，正在读取:", group_file))
  
  # 读取分组表 (假设第一列是样本名，第二列是 Group)
  tryCatch({
    group_df <- read.csv(group_file, header = TRUE)
    # 强制确保只有两列，或者是指定的列名
    # 这里假设第一列就是 ID，第二列就是 Group
    if(ncol(group_df) < 2) stop("分组文件列数不足，需要至少两列(样本名, 分组)")
    
    # 清洗数据
    target_ids   <- as.character(group_df[,1]) # 样本名
    target_group <- as.character(group_df[,2]) # 分组
    
    # 校验：分组文件里的样本是否涵盖了表达矩阵里的样本？
    # 这一步很关键，因为用户上传的文件顺序可能不一样
    match_idx <- match(sample_names, target_ids)
    
    if (any(is.na(match_idx))) {
      stop("错误：表达矩阵中有样本在分组文件中找不到！请检查样本名是否一致。")
    }
    
    # 按表达矩阵的顺序提取分组
    group_vec <- target_group[match_idx]
    
    # 标准化分组名称 (把 0/1/Control/Treat 都映射到 Normal/Tumor)
    # 简单起见，这里直接用用户文件里的内容，但要求用户尽量规范
    # 为了后续计算，我们将这一列转因子，并且尝试将 'Normal/Control' 设为基准
    
    final_group <- factor(group_vec)
    
    # 尝试设置基准水平 (Ref) 为 Normal 或 Control
    ref_level <- grep("normal|control|healthy|0", levels(final_group), ignore.case = TRUE, value = TRUE)
    if (length(ref_level) > 0) {
      final_group <- relevel(final_group, ref = ref_level[1])
      print(paste("已设置参考组为:", ref_level[1]))
    } else {
      print(paste("未检测到标准参考组名，默认参考组为:", levels(final_group)[1]))
    }
    
  }, error = function(e) {
    stop("读取分组文件失败: ", e$message)
  })
  
} else {
  # --- 模式 B: 自动识别 ---
  print("未提供分组文件，启动自动识别 (AUTO)...")
  
  # 策略 1: TCGA 格式
  is_tcga <- all(grepl("TCGA", sample_names, ignore.case = TRUE)) && 
    all(grepl("-\\d{2}[A-Z]?$", sample_names))
  
  if (is_tcga) {
    print("  -> 检测到 TCGA 格式 ID")
    codes <- as.numeric(substring(sample_names, nchar(sample_names)-2, nchar(sample_names)-1))
    group_vec <- ifelse(codes < 10, "Tumor", "Normal")
  } else {
    # 策略 2: 关键词
    print("  -> 尝试关键词识别")
    group_vec <- rep(NA, length(sample_names))
    lower_names <- tolower(sample_names)
    
    idx_normal <- grepl("normal|control|healthy|adj", lower_names)
    idx_tumor  <- grepl("tumor|cancer|treat|case", lower_names)
    
    group_vec[idx_normal] <- "Normal"
    group_vec[idx_tumor]  <- "Tumor"
    
    if (any(is.na(group_vec))) {
      stop("自动分组失败！请编辑正确格式或上传分组文件。")
    }
  }
  
  final_group <- factor(group_vec, levels = c("Normal", "Tumor"))
}

print("分组统计:")
print(table(final_group))

# 校验分组是否有至少两个水平
if (length(levels(final_group)) < 2) {
  stop("错误：分组只有 1 种类别，无法进行差异分析！")
}

# ====================================================================
#  4. 差异分析计算 
# ====================================================================

# (1) 设计矩阵
design <- model.matrix(~0 + final_group)
colnames(design) <- levels(final_group) # 清理列名

# (2) 过滤与归一化
dge <- DGEList(counts = rnaCounts, group = final_group)
keep <- rowSums(cpm(dge) > 1) >= 2
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# (3) 拟合
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)

# (4) 构建对比矩阵 (动态构建)
# 我们需要找到哪两个组进行对比。通常是 Non-Ref vs Ref
# 假设 levels 是 A, B。如果不指定，limma 默认比较复杂的 coef。
# 这里动态生成一个 "B - A" 的对比
grps <- colnames(design)
contrast_formula <- paste(grps[2], "-", grps[1]) # 比如 "Tumor - Normal"
print(paste("正在进行对比:", contrast_formula))

contrast.matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# (5) 结果提取
DEG_Result <- topTable(fit2, coef = 1, n = Inf)

# 1. 把行名（基因 ID）变成独立的一列，起名叫 "gene_symbol"
DEG_Result$gene_symbol <- rownames(DEG_Result)

# 2. 重新排列列顺序，把基因名放第一列，方便 Java 读取
DEG_Result <- DEG_Result[, c("gene_symbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

# 3. 保存 CSV 时，去掉 row.names = FALSE，Java 读到的第一行就是表头
output_csv <- file.path(output_dir, "diff_result.csv")
write.csv(DEG_Result, output_csv, row.names = FALSE, quote = FALSE) 

print(paste("数据已保存:", output_csv))

# ====================================================================
#  5. 绘图 (动态适配颜色)
# ====================================================================
print("正在绘制火山图...")

# 动态生成标题 (防止用户组名不是 Tumor/Normal)
plot_title <- paste(levels(final_group)[2], "vs", levels(final_group)[1])

tryCatch({
  # 设置颜色逻辑
  keyvals <- ifelse(DEG_Result$logFC < -1 & DEG_Result$adj.P.Val < 0.05, 'blue',
                    ifelse(DEG_Result$logFC > 1 & DEG_Result$adj.P.Val < 0.05, 'red', 'grey'))
  
  # 防止全是灰色导致报错
  keyvals[is.na(keyvals)] <- 'grey'
  
  names(keyvals)[keyvals == 'red'] <- 'Up'
  names(keyvals)[keyvals == 'grey'] <- 'NS'
  names(keyvals)[keyvals == 'blue'] <- 'Down'
  
  p <- EnhancedVolcano(DEG_Result,
                       lab = rownames(DEG_Result),
                       x = 'logFC',
                       y = 'adj.P.Val',
                       title = plot_title,           # 动态标题
                       subtitle = 'Differential Expression',
                       pCutoff = 0.05,
                       FCcutoff = 1,
                       pointSize = 2.0,              # 你原本的大小
                       labSize = 3.0,                # 你原本的大小
                       colCustom = keyvals,
                       legendPosition = 'top')       # 你原本的位置
  
  ggsave(file.path(output_dir, "volcano_plot.png"), plot = p, width = 8, height = 8, dpi = 150)
  print(paste("火山图已保存:", file.path(output_dir, "volcano_plot.png")))
  
}, error = function(e) {
  print(paste("绘图警告:", e$message))
})

print("分析完成。")