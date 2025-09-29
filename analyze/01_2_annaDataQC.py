import sys
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

Dpi00_1 =anndata.read_h5ad("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/processData/Dpi00_1_filtered.h5ad")
Dpi00_2 =anndata.read_h5ad("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/processData/Dpi00_2_filtered.h5ad")
Dpi02_1 =anndata.read_h5ad("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/processData/Dpi02_1_filtered.h5ad")
Dpi02_2 =anndata.read_h5ad("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/processData/Dpi02_2_filtered.h5ad")
Dpi05_1 =anndata.read_h5ad("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/processData/Dpi05_1_filtered.h5ad")
Dpi05_2 =anndata.read_h5ad("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/processData/Dpi05_2_filtered.h5ad")
Dpi10_1 =anndata.read_h5ad("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/processData/Dpi10_1_filtered.h5ad")
Dpi10_2 =anndata.read_h5ad("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/processData/Dpi10_2_filtered.h5ad")
Dpi21_1 =anndata.read_h5ad("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/processData/Dpi21_1_filtered.h5ad")
Dpi21_2 =anndata.read_h5ad("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/processData/Dpi21_2_filtered.h5ad")

Dpi00 = Dpi00_1.concatenate(Dpi00_2)
Dpi02 = Dpi02_1.concatenate(Dpi02_2)
Dpi05 = Dpi05_1.concatenate(Dpi05_2)
Dpi10 = Dpi10_1.concatenate(Dpi10_2)
Dpi21 = Dpi21_1.concatenate(Dpi21_2)

samples=["Dpi00", "Dpi02", "Dpi05", "Dpi10", "Dpi21"]
fig, aex = plt.subplots(ncols=2, nrows=5, figsize=(50, 30))
for idx, data in enumerate([Dpi00, Dpi02, Dpi05, Dpi10, Dpi21]):
    data_concatenated = data
    data.write_h5ad(f"processData/{samples[idx]}_concatenated.h5ad")
    data.layers['counts'] = data.X.copy()
    sc.pp.normalize_total(data,target_sum=1e4,inplace=True)
    sc.pp.log1p(data)
    data_scaled = data.copy()
    data.write_h5ad(f"processData/{samples[idx]}_scaled.h5ad")
    sns.histplot(data_concatenated.obs["total_counts"], bins=100,kde=True,ax=aex[idx,0])
    aex[idx,0].set_title(f"{samples[idx]}_Total counts")
    sns.histplot(data_scaled.X.sum(1), bins=100,kde=True,ax=aex[idx,1])
    aex[idx,1].set_title(f"{samples[idx]}_Shifted logarithm")
plt.savefig("figures/hist_counts.png")

adata_list =[]
for sample in samples:
    ad = anndata.read_h5ad(f"processData/{sample}_scaled.h5ad")
    ad.obs["Sample"] = sample
    adata_list.append(ad)

data_all = anndata.concat(adata_list)
data_all.write_h5ad("_processData/data_all.h5ad")
