import sys
import anndata
import scanpy as sc
import scvi

data_scaled = anndata.read_h5ad("./_processData/data_all.h5ad")
data = data_scaled

sc.experimental.pp.highly_variable_genes(
    data, 
    flavor="pearson_residuals", # pearson残差的方法来识别高变基因
    layer='counts',
    batch_key = "Sample",
    n_top_genes=5000, #选出变化最大的前5000个基因，其余过滤掉
    subset=True,
    inplace=True,
)

# 这里是去批次效应步骤，方法为SCVI
scvi.model.SCVI.setup_anndata(data, layer = "counts", batch_key="Sample",  # scVI的必要步骤，告诉模型如何解析AnnData中的数据
                             categorical_covariate_keys=["Sample"],  #sample 是一个分类协变量
                             continuous_covariate_keys=['total_counts']) # total_counts 是一个连续协变量，每个细胞的总分子数。

scvi.settings.seed = 7
model = scvi.model.SCVI(data)
model.train()
# 提取训练后模型中的低纬潜在向量
latent = model.get_latent_representation()
data.obsm['X_scVI'] = latent
# 获取模型中的归一化表达
data.layers['scvi_normalized'] = model.get_normalized_expression(library_size = 1e4)
# 保存数据
data.write_h5ad("./_processData/data_integrated.h5ad")
