import scanpy as sc
import anndata as ad

adata = sc.read_h5ad("/home/zwye/rare/EPI/Cervical_EPI_MAC_sub_Cer20-21-35-36-37.h5ad")

# 1) 用原始或 counts 层做 log1p（不要先 scale）
# 如果你的原始 counts 在 adata.raw 或某个 layer，请先把它放到 X 再操作
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 2) 记录一份未校正的矩阵到 layer（可保留给 scCAD）
adata.layers["log1p"] = adata.X.copy()

# 3) ComBat：选择合适的批次键，比如 'SampleID' 或 'Patients'
sc.pp.combat(adata, key='Patients')  # 会直接作用于 adata.X，可能产生负值

# 5) 保存：combat 后的用于可视化；未校正log1p在 adata.layers["log1p"] 可另存留给 scCAD
adata.write_h5ad("/home/zwye/rare/EPI/Cervical_EPI_MAC_combat_visual.h5ad")

