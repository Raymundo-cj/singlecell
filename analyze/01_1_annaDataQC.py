import sys
import anndata
import pandas as pd
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
import scripts
import scripts.scDblFinder
from scripts.scDblFinder import run_ScDblFinder
import matplotlib.pyplot as plt

Dpi00_1 = sc.read_10x_mtx("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/dpi0_rep1/filtered_feature_bc_matrix", cache=True)
Dpi00_2 = sc.read_10x_mtx("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/dpi0_rep2/filtered_feature_bc_matrix", cache=True)
Dpi02_1 = sc.read_10x_mtx("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/dpi2_rep1/filtered_feature_bc_matrix", cache=True)
Dpi02_2 = sc.read_10x_mtx("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/dpi2_rep2/filtered_feature_bc_matrix", cache=True)
Dpi05_1 = sc.read_10x_mtx("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/dpi5_rep1/filtered_feature_bc_matrix", cache=True)
Dpi05_2 = sc.read_10x_mtx("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/dpi5_rep2/filtered_feature_bc_matrix", cache=True)
Dpi10_1 = sc.read_10x_mtx("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/dpi10_rep1/filtered_feature_bc_matrix", cache=True)
Dpi10_2 = sc.read_10x_mtx("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/dpi10_rep2/filtered_feature_bc_matrix", cache=True)
Dpi21_1 = sc.read_10x_mtx("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/dpi21_rep1/filtered_feature_bc_matrix", cache=True)
Dpi21_2 = sc.read_10x_mtx("/share/org/YZWL/yzwl_caojian/caojian/singlecell/zhongwan/dpi21_rep2/filtered_feature_bc_matrix", cache=True)

def run_plot_scatter(adata, x, y, sample, ax):

    sc.pl.scatter(adata, x=x, y=y, ax=ax)
    ax.set_title(sample)
    ax.set_xlim(0, 50000)
    ax.set_ylim(0, 16000)

samples = ["Dpi00_1", "Dpi00_2","Dpi02_1", "Dpi02_2","Dpi05_1", "Dpi05_2", "Dpi10_1", "Dpi10_2", "Dpi21_1", "Dpi21_2"]
data_list = [Dpi00_1, Dpi00_2, Dpi02_1, Dpi02_2, Dpi05_1, Dpi05_2, Dpi10_1, Dpi10_2, Dpi21_1, Dpi21_2]
raw_data_list = []

fig, aex = plt.subplots(ncols=10, nrows=3, figsize=(50, 30), dpi = 300)
for idx, data in enumerate(data_list):

    data.obs["Sample"] = samples[idx]
    sc.pp.calculate_qc_metrics(data, inplace=True, percent_top=None, log1p=False)

    data_raw = data.copy()
    raw_data_list.append(data_raw)

    # 绘制原始数据散点图
    run_plot_scatter(data, 'total_counts', 'n_genes_by_counts', samples[idx], aex[0, idx])

    # 去除双胞
    run_ScDblFinder(data, copy=False, doubletRatio=0.1)
    # 绘制去除双胞后的散点图
    run_plot_scatter(data, 'total_counts', 'n_genes_by_counts', samples[idx], aex[1, idx])
    # 手动过滤低质量细胞和基因
    sc.pp.filter_cells(data, min_genes=300)
    sc.pp.filter_cells(data, min_counts=600)
    sc.pp.filter_cells(data, max_genes=6000)
    sc.pp.filter_cells(data, max_counts=20000)
    # 绘制过滤后的散点图
    run_plot_scatter(data, 'total_counts', 'n_genes_by_counts', samples[idx], aex[2, idx])
    # 将数据写入文件
    filename = f"processData/{samples[idx]}_filtered.h5ad"
    data.write_h5ad(filename)
# 保存图形
plt.savefig("figures/scatter.pdf")

# 统计过滤前后的细胞数、基因数、分子数
data_pairs = list(zip(raw_data_list, data_list))

####statitcs
colors = ['#1f77b4', '#ff7f0e']

cell_counts = []
sample_names = []

for i, (data_a, data_b) in enumerate(data_pairs):
    cell_counts.append([data_a.n_obs, data_b.n_obs])
    sample_names.append([f'raw', f'filtered'])

fig, axes = plt.subplots(1, 10, figsize=(40, 6), sharey=True, dpi = 300)

for i, ax in enumerate(axes):
    bars = ax.bar(sample_names[i], cell_counts[i], color=colors)
    ax.set_title(f'{samples[i]}')
    if i == 0:
        ax.set_ylabel('Cell counts',fontsize=24)

    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2.0, height, f'{height}', ha='center', va='bottom', fontsize=20)
plt.savefig("figures/cellCount.pdf")

##
gene_counts = []
sample_names = ["raw","filtered"]

for i, (data_a, data_b) in enumerate(data_pairs):
    gene_counts.append(data_a.obs["n_genes_by_counts"].tolist())
    gene_counts.append(data_b.obs["n_genes_by_counts"].tolist())

medianprops = dict(color='black', linewidth=2)

fig, axes = plt.subplots(1, 10, figsize=(40, 6), sharey=True, dpi = 300)

for i, ax in enumerate(axes):
    data_to_plot = [gene_counts[2*i], gene_counts[2*i+1]]
    sns.violinplot(data=data_to_plot, ax=ax, inner=None, palette=colors[:2])

    #for patch, color in zip(box['boxes'], colors):
    #    patch.set_facecolor(color)

    ax.set_ylim(0, 20000)

    medians = [np.median(data) for data in data_to_plot]
    positions = np.arange(len(data_to_plot)+1)
    for pos, median in zip(positions, medians):
        ax.text(pos, median, f'{median}', verticalalignment='center', fontsize=20, color='black', ha='left')

    ax.set_xticks(positions)

    ax.set_title(f'{samples[i]}')
    if i == 0:
        ax.set_ylabel('Gene counts',fontsize=24)

plt.savefig("figures/geneCount.pdf")

##
umi_counts = []
sample_names = ["raw","filtered"]

for i, (data_a, data_b) in enumerate(data_pairs):
    umi_counts.append(data_a.obs["total_counts"].tolist())
    umi_counts.append(data_b.obs["total_counts"].tolist())

medianprops = dict(color='black', linewidth=2)

fig, axes = plt.subplots(1, 10, figsize=(40, 6), sharey=True, dpi = 300)

for i, ax in enumerate(axes):
    data_to_plot = [umi_counts[2*i], umi_counts[2*i+1]]
    sns.violinplot(data=data_to_plot, ax=ax, inner=None, palette=colors)

    #for patch, color in zip(box['boxes'], colors):
    #    patch.set_facecolor(color)

    ax.set_ylim(0, 80000)

    medians = [np.median(data) for data in data_to_plot]
    positions = np.arange(0, len(data_to_plot)+1)
    for pos, median in zip(positions, medians):
        ax.text(pos, median, f'{median}', verticalalignment='center', fontsize=20, color='black', ha='left')

    ax.set_xticks(positions)

    ax.set_title(f'{samples[i]}')
    if i == 0:
        ax.set_ylabel('Umi counts',fontsize=24)

plt.savefig("figures/umiCount.pdf")
