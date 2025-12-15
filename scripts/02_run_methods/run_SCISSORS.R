library(SingleCellExperiment)
library(Seurat)
library(SCISSORS)

sce <- readRDS("/home/zwye/rare/MicroCellClust/version_1/sce_P002.Pre.rds")
# 1. 看看 sce 里有哪些 assay
assayNames(sce)
# 例如： [1] "counts" "logcounts" 或者只有 "counts"
logcounts(sce) <- log1p(assay(sce, "counts"))
# 2a. 如果你已有标准化 / log1p 结果在 assay "logcounts"
seu <- as.Seurat(
  sce,
  counts = "counts",      # 原始计数
  data   = "logcounts",   # 归一化或 log1p 结果
  project = "RareCells"
)
# 3. 初始预处理 + 聚类 + UMAP
seu <- PrepareData(
  seu,
  n.HVG              = 2000,    # 高变基因数量，2000–5000 可调
  n.PC               = 20,      # PCA 维度
  which.dim.reduc    = "umap",  
  initial.resolution = 0.5,     # Louvain 初始分辨率
  use.parallel       = FALSE,
  random.seed        = 42
)
seu <- PrepareData(
  seurat.object   = seu,
  which.dim.reduc = "pca",
  use.parallel    = FALSE,  # 其余仍走默认
  random.seed     = 42
) #多用默认值
#> [1] "Running UMAP on 20 principal components"
#> [1] "Found X unique clusters"

# ─── 4. 计算每个细胞的 Silhouette Score ────────────────
sil_df <- ComputeSilhouetteScores(seu, avg = FALSE)
# sil_df 有两列：Cluster、Score；行名是细胞 ID
sil_df$Cell <- rownames(sil_df)

# ─── 5. 定义稀有阈值（<1% 总细胞数）────────────────────
total_cells <- ncol(seu)
thr         <- ceiling(0.01 * total_cells)

# ─── 6. 收集“初始就很小”的簇里的所有细胞 ───────────────
init_sizes        <- table(Idents(seu))
rare_init_clusters <- names(init_sizes)[ init_sizes < thr ]

rare_cells <- c()
if (length(rare_init_clusters) > 0) {
  rare_cells <- WhichCells(seu, idents = rare_init_clusters)
}

# ─── 7. 对其余簇逐一做子聚类，并收集小于阈值的子簇 ───────
other_clusts <- setdiff(names(init_sizes), rare_init_clusters)

for (cl in other_clusts) {
  reclust_i <- ReclusterCells(
    seu,
    which.clust       = as.integer(cl),
    merge.clusters    = TRUE,
    k.vals            = c(10, 25, 50),     # 可根据该簇规模调
    resolution.vals   = c(0.2, 0.3, 0.4),
    n.HVG             = 2000,
    n.PC              = 20,
    redo.embedding    = TRUE,
    use.parallel      = FALSE,
    random.seed       = 123
  )
  reclust_i <- ReclusterCells(
  seurat.object  = seu,
  which.clust    = as.integer(cl),
  redo.embedding = FALSE,   # 不在子集里再做 UMAP，稳
  use.parallel   = FALSE    # 其余仍走默认
  ) #用默认参数版
  # 如果真的拆出了多个子簇
  sub_ids <- unique(as.character(Idents(reclust_i)))
  if (length(sub_ids) > 1) {
    sub_sizes <- table(Idents(reclust_i))
    rare_sub  <- names(sub_sizes)[ sub_sizes < thr ]
    if (length(rare_sub) > 0) {
      rare_cells <- c(
        rare_cells,
        WhichCells(reclust_i, idents = rare_sub)
      )
    }
  }
}

# ─── 8. 去重并查看结果 ─────────────────────────────────
rare_cells <- unique(rare_cells)
cat("共检测到稀有/异常细胞：", length(rare_cells), "\n")
head(rare_cells, 20)

# ─── 9. （可选）在 seu 元数据里打标记 ───────────────────
seu$IsRare <- colnames(seu) %in% rare_cells
table(seu$IsRare)

# ─── 10. （可选）导出到 CSV ────────────────────────────
write.csv(
  data.frame(Cell = rare_cells),
  file = "SCISSORS_rare_cells.csv",
  row.names = FALSE
)
