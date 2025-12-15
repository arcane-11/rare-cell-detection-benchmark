# ===== 0) 加载依赖 =====
library(data.table)
library(sva)

# ===== 1) 文件路径 =====
expr_file  <- "/home/zwye/rare/hippocampus/expression_matrix_full.csv"
meta_file  <- "/home/zwye/rare/hippocampus/cell_metadata.csv"
out_file   <- "/home/zwye/rare/hippocampus/expression_ComBatSeq_corrected.csv"

# ===== 2) 读取表达矩阵（细胞 × 基因）=====
cat("[INFO] 正在读取表达矩阵...\n")
expr_raw <- fread(expr_file)
125
# 假设第一列是 CellID
cell_ids <- expr_raw[[1]]
expr <- as.matrix(expr_raw[,-1,with=FALSE])
mode(expr) <- "numeric"
rownames(expr) <- cell_ids
cat("[INFO] 表达矩阵维度:", dim(expr), "(细胞 × 基因)\n")

# 转置为基因 × 细胞（ComBat-seq 要求）
expr_t <- t(expr)
cat("[INFO] 转置后维度:", dim(expr_t), "(基因 × 细胞)\n")

# ===== 3) 读取元信息 =====
cat("[INFO] 正在读取元信息...\n")
meta <- fread(meta_file)

# 对齐：metadata 的 CellID 和表达矩阵的行名一致
stopifnot("CellID" %in% colnames(meta))
stopifnot(all(rownames(expr) %in% meta$CellID))

# 按照表达矩阵细胞顺序对齐 metadata
meta <- meta[match(rownames(expr), meta$CellID), ]
stopifnot(all(meta$CellID == rownames(expr)))

# ===== 4) 批次向量 =====
# ⚠️ 这里你可以改成 "ChipID" / "Flowcell" / "DonorID"
batch_col <- "SampleID"
batch <- meta[[batch_col]]

# ===== 5) 运行 ComBat-seq =====
cat("[INFO] 运行 ComBat-seq (batch =", batch_col, ")...\n")
expr_adj <- ComBat_seq(counts = expr_t, batch = batch)

# ===== 6) 保存结果（基因 × 细胞）=====
write.csv(
  cbind(gene=rownames(expr_adj), expr_adj),
  file = out_file,
  row.names = FALSE,
  quote = FALSE
)
cat("[INFO] ComBat-seq 完成，已保存到:", out_file, "\n")

)

