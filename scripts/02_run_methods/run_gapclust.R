# 1. 读取表达矩阵（没有行名）
expr_mat <- read.delim(
  "/home/zwye/rare/GSE157783/GSE157783_IPDCO_hg_midbrain_cell/split_by_patient/expr_C1.tsv",
  header = TRUE,        # 第一行是细胞名
  row.names = NULL,     # 不用第一列做行名
  check.names = FALSE
)

# 2. 读取基因文件
gene_info <- read.delim("/home/zwye/rare/GSE157783/IPDCO_hg_midbrain_genes.tsv", sep = "\t", check.names = FALSE)

# 3. 用 gene 列作为表达矩阵行名
rownames(expr_mat) <- gene_info$gene

# 4. 转为矩阵（如果后续需要）
expr_mat <- as.matrix(expr_mat)
expr_mat[is.na(expr_mat)] <- 0
gene_vars <- rowVars(expr_mat)
expr_mat <- expr_mat[gene_vars > 0, , drop = FALSE]
rownames(expr_mat) <- gsub("_", "-", rownames(expr_mat))

expr_sparse <- as(expr_mat, "dgCMatrix")
res <- GapClust(expr_sparse, k = k)
# 选取了总细胞数的1%
# 提取稀有细胞索引
rare_ids   <- unique(unlist(res$rare_cell_indices))
rare_cells <- colnames(expr_sparse)[rare_ids]

# 保存到文件
write.csv(data.frame(cell = rare_cells), 
          file = "rare_cells.csv", 
          row.names = FALSE, quote = FALSE)

