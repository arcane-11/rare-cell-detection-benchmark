# 设置环境变量
Sys.setenv(PATH = paste0(
  Sys.getenv("HOME"), "/conda-envs/cs_mcl/bin:",
  Sys.getenv("PATH")
))
mcl_exec <- Sys.which("mcl")
print(mcl_exec)
if (mcl_exec == "") {
  stop("仍然找不到 mcl，请检查 PATH 是否正确设置。")
}

# 1. 先读入数据，不设置 row.names
expr_df <- read.csv(
  "/home/zwye/rare/hippocampus/expression_ComBatSeq_corrected.csv",
  header = TRUE,
  check.names = FALSE
)

# 2. 提取第一列作为基因名
gene_names <- expr_df[[1]]

# 3. 去掉第一列（只保留数值部分）
expr_mat <- as.matrix(expr_df[,-1])

# 4. 设置行名，自动处理重复
rownames(expr_mat) <- make.unique(gsub("_", "-", gene_names))

# 5. 确认列名是否存在，如果缺失则自动生成
if (is.null(colnames(expr_mat)) || any(colnames(expr_mat) == "")) {
  colnames(expr_mat) <- paste0("Cell", seq_len(ncol(expr_mat)))
}

# 6. 确保类型是 matrix
expr_mat <- structure(expr_mat, class = "matrix")

# 7. 检查维度
cat("Matrix dimensions: genes =", nrow(expr_mat), ", cells =", ncol(expr_mat), "\n")

stopifnot(is.matrix(expr_mat))
identical(class(expr_mat), "matrix")

# -----------------------------
# Seurat 初始聚类
# -----------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(CellSIUS)
  library(data.table)
})

set.seed(1234)

seu <- CreateSeuratObject(counts = expr_mat)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 20)
seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.5)
initial_clusters <- as.numeric(Idents(seu))

# -----------------------------
# 修正 melt 参数，避免报错
# -----------------------------
# CellSIUS 内部会调用 melt，我们自己给 expr_dt 构造一个标准格式
expr_dt <- as.data.table(expr_mat)
expr_dt[, gene_id := rownames(expr_mat)]
setcolorder(expr_dt, c("gene_id", setdiff(names(expr_dt), "gene_id")))

expr_long <- melt(
  expr_dt,
  id.vars = "gene_id",
  variable.name = "cell_idx",
  value.name = "expr"
)

# -----------------------------
# 运行 CellSIUS
# -----------------------------
cs_out <- CellSIUS(
  mat.norm          = expr_mat,
  group_id          = initial_clusters,
  min_n_cells       = 3,
  min_fc            = 1,
  corr_cutoff       = 0.7,
  fc_between_cutoff = 1,
  mcl_path          = mcl_exec
)

# -----------------------------
# 整理结果
# -----------------------------
dt <- as.data.table(cs_out)
assign_df <- dt[, .(
  main_cluster = unique(main_cluster),
  sub_cluster  = paste(unique(sub_cluster), collapse = ";")
), by = cell_idx]

setnames(assign_df, "cell_idx", "cell_id")

out_file <- "/home/zwye/rare/hippocampus/cellsius/combat_CellSIUS_rare_cell_assignments.csv"
fwrite(assign_df, out_file)

cat("✅ CellSIUS 结果已保存至:", out_file, "\n")
