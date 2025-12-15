import pandas as pd
import anndata
import scanpy as sc
import matplotlib.pyplot as plt
# 1. 读取 tsv，原始矩阵是 行=基因，列=细胞
df = pd.read_csv(
    "/home/zwye/rare/GSE266919/sce_merged_10_expression_matrix.csv",
    index_col=0
)


print("原始矩阵形状 (基因 x 细胞) =", df.shape)

# 2. 转置：变为 行=细胞，列=基因
df = df.T
print("转置后矩阵形状 (细胞 x 基因) =", df.shape)

# 3. 构造 AnnData
adataSC = anndata.AnnData(
    X=df.values,
    obs=pd.DataFrame(index=df.index),   # 细胞 ID
    var=pd.DataFrame(index=df.columns)  # 基因名
)

print(adataSC)

# 3. （可选）做基因过滤，比如保留至少出现在 200 个细胞里的基因：
sc.pp.filter_genes(adataSC, min_cells=200)

# 3. 把 X 转为 float
adataSC.X = adataSC.X.astype("float64")

# 4. 归一化：每个细胞总计数缩放到 10,000
sc.pp.normalize_per_cell(adataSC, counts_per_cell_after=1e4)

# 5. 对数化
sc.pp.log1p(adataSC)
# 注意：后续 GiniClust3 的函数都会直接对 adataSC.X 进行计算，并在 adataSC.obs 中添加对应的标签/结果列
import giniclust3 as gc
# 计算所有基因的 Gini 指数，结果会写入 adataSC.obs['gini']（或类似字段，具体请参照文档）
gc.gini.calGini(adataSC)
# gc.gini.calGini(adataSC, min_gini_value=0.3)
# 传入邻居数 neighbors（邻居数一般不需要太大，默认为 3~10 之间），对高 Gini 基因做聚类
# adataGini = gc.gini.clusterGini(adataSC, neighbors=10)
adataGini = gc.gini.clusterGini(adataSC)
# 运行结束后，adataSC.obs 中会新增一列 'rare'，用于标记基于 Gini 聚类得到的罕见簇标签
# adataGini.obs['rare'] 即为 Gini 聚类结果中每个细胞对应的“罕见簇”编号
#import importlib, giniclust3
#import giniclust3.fano

#importlib.reload(giniclust3)         # 重新加载顶层包
#importlib.reload(giniclust3.fano) 
# 计算所有基因的 Fano 因子，结果会写入 adataSC.obs['fano']（或类似字段）
gc.fano.calFano(adataSC)

# 基于高 Fano 因子的基因做聚类
adataFano = gc.fano.clusterFano(adataSC)
# 运行结束后，adataSC.obs 中会新增一列 'fano'，用于标记基于 Fano 聚类得到的“常见簇”标签
# adataFano.obs['fano'] 即为 Fano 聚类结果中每个细胞对应的簇编号
import numpy as np
# 构建一个字典，存放 Gini 与 Fano 两次聚类得到的簇标签
consensusCluster = {}
# 'rare' 对应 Gini 聚类标签，adataSC.obs['rare'] 是一个 list-like，要转换成 numpy 数组
consensusCluster['giniCluster'] = np.array(adataSC.obs['rare'].values.tolist())
# 'fano' 对应 Fano 聚类标签
consensusCluster['fanoCluster'] = np.array(adataSC.obs['fano'].values.tolist())

# 生成共识矩阵（M-tilde），内部会将 giniCluster 与 fanoCluster 的信息融合，生成一个加权的共识矩阵
gc.consensus.generateMtilde(consensusCluster)

# 基于共识矩阵再做一次聚类（默认使用谱聚类或者 Louvain），得到最终的簇标签
gc.consensus.clusterMtilde(consensusCluster)

# 最终的簇标签保存在 consensusCluster['finalCluster']，形如 array([...], dtype=object)
# 将它写到一个文本文件，或者直接存回 AnnData 中
np.savetxt("final.txt", consensusCluster['finalCluster'], delimiter="\t", fmt='%s')
# 先计算邻接图，再运行 UMAP（使用 Scanpy）
sc.pp.pca(adataSC, n_comps=50)
sc.pp.neighbors(adataSC, n_neighbors=15, n_pcs=50)
sc.tl.umap(adataSC)

# 将最终簇标签加入 obs
# 之前我们已将 adataSC.obs['final_cluster'] 赋好

# 用 Scanpy 画 UMAP，按 final_cluster 着色

sc.pl.umap(adataSC, color='rare', size=20, legend_loc='on data')
plt.show()
labels_final = consensusCluster['finalCluster']  # 这是一个长度为 17149 的字符串数组

# 将它写回到 adataSC.obs，列名可以自定义，这里用 'final_cluster'
adataSC.obs['final_cluster'] = labels_final

# 统计每个 final_cluster 下的细胞数
cluster_counts = adataSC.obs['final_cluster'].value_counts()
total = adataSC.n_obs
cluster_freq = pd.DataFrame({
    'cluster': cluster_counts.index,
    'count': cluster_counts.values,
    'fraction': cluster_counts.values / total
})
# 找出 fraction < 0.01（即 <1%）的那些簇，视为罕见簇
rare_clusters = cluster_freq[cluster_freq['fraction'] < 0.01]['cluster'].tolist()
print("罕见簇标签（<1% 细胞）：", rare_clusters)

