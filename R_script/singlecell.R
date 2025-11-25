library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(cowplot)
library(future)
library(sceasy)
library(reticulate)
library(anndata)
library(glmGamPoi)

options(future.globals.maxSize = 20 * 1024^3)
# 定义路径和标签
sample_dirs <- c(
  "/share/org/YZWL/yzwl_caojian/caojian/singlecell/singlecell_1103/ChiHei/dpi00_rep1/filtered_feature_bc_matrix/",
  "/share/org/YZWL/yzwl_caojian/caojian/singlecell/singlecell_1103/ChiHei/dpi00_rep2/filtered_feature_bc_matrix/",
  "/share/org/YZWL/yzwl_caojian/caojian/singlecell/singlecell_1103/ChiHei/dpi02_rep1/filtered_feature_bc_matrix/",
  "/share/org/YZWL/yzwl_caojian/caojian/singlecell/singlecell_1103/ChiHei/dpi02_rep2/filtered_feature_bc_matrix/",
  "/share/org/YZWL/yzwl_caojian/caojian/singlecell/singlecell_1103/ChiHei/dpi10_rep1/filtered_feature_bc_matrix/",
  "/share/org/YZWL/yzwl_caojian/caojian/singlecell/singlecell_1103/ChiHei/dpi10_rep2/filtered_feature_bc_matrix/",
  "/share/org/YZWL/yzwl_caojian/caojian/singlecell/singlecell_1103/ChiHei/dpi21_rep1/filtered_feature_bc_matrix/",
  "/share/org/YZWL/yzwl_caojian/caojian/singlecell/singlecell_1103/ChiHei/dpi21_rep2/filtered_feature_bc_matrix/"
)

sample_list <- c(
  "dpi00_rep1","dpi00_rep2",
  "dpi02_rep1","dpi02_rep2",
  "dpi10_rep1","dpi10_rep2",
  "dpi21_rep1","dpi21_rep2"
)

stage_labels <- c(
    "dpi00","dpi00",
    "dpi02","dpi02",
    "dpi10","dpi10",
    "dpi21","dpi21")

# 对数据进行质控
objs <- list()
for (i in 1:length(sample_list)){
    data <- Read10X(data.dir = sample_dirs[i])
    TenX_data <- CreateSeuratObject(counts =data,project=sample_list[i],min.cells=3,min.features=200)
    TenX_data$orig.ident <- sample_list[i]
    TenX_data$stage <- stage_labels[i]

    quantile(TenX_data$nFeature_RNA,c(0.05,0.95))
    # normalization
    TenX_data <- subset(TenX_data, subset = nFeature_RNA > 200)
    TenX_data <- NormalizeData(TenX_data, normalization.method = "LogNormalize", scale.factor = 10000)
    
    objs[[i]] <- TenX_data
    assign(sample_list[i],TenX_data)}

merge_data <- merge(objs[[1]],y=objs[-1],add.cell.ids=sample_list)
saveRDS(merge_data,file="./RDS_file/data_all.rds")

tenscRNA <- readRDS("./RDS_file/data_all.rds")
scRNA <- FindVariableFeatures(tenscRNA, selection.method = "vst", nfeatures = 3000)
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA)
scRNA <- RunPCA(scRNA)

# harmony method
seurat_merge <- IntegrateLayers(object=scRNA,method = HarmonyIntegration,
  orig.reduction ="pca",new.reduction ="harmony",verbose = FALSE)
dm =10
seurat_merge <- FindNeighbors(seurat_merge,reduction = "harmony", dims = 1:dm)
seurat_merge <- RunUMAP(seurat_merge, reduction = "harmony",dims = 1:dm, min.dist = 0.1,n.neighbors = 30)


saveRDS(seurat_merge,file="./RDS_file/data_integrated_harmony.rds")  # 保存聚类之前的数据文件作为中间文件
data <- readRDS("./RDS_file/data_integrated_harmony.rds")

# cluster
resolution <- c(0.2, 0.4, 0.6, 0.8, 1.0)
for (res in resolution){
    data <- FindClusters(data,resolution=res,algorithm = 4, random.seed = 0)
    data[[paste0('leiden_',res)]] <- data$seurat_clusters
}

# stage
p1 <- DimPlot(data, reduction = "umap",group.by = "stage",pt.size=0.1,alpha=c(0.6))
ggsave("umap_stage_harmony.pdf",plot = p1, width = 7, height = 6, dpi = 300)

p2 <- DimPlot(data,reduction="umap",group.by="leiden_0.4",label = T,pt.size=0.1,alpha=c(0.6))
p3 <- DimPlot(data,reduction="umap",group.by="leiden_0.6",label = T,pt.size=0.1,alpha=c(0.6))
p4 <- DimPlot(data,reduction="umap",group.by="leiden_0.8",label = T,pt.size=0.1,alpha=c(0.6))
p5 <- DimPlot(data,reduction="umap",group.by="leiden_1",label = T,pt.size=0.1,alpha=c(0.6))

p6 <- plot_grid(p2,p3,p4,p5,ncol=2,nrow=2)
ggsave("umap_clusters_harmony.pdf",plot = p6, width = 14, height = 12, dpi = 300)
saveRDS(data,file="./RDS_file/data_cluster.rds")

tenscRNA <- readRDS("./RDS_file/data_all.rds")
data_cluster <- readRDS("./RDS_file/data_cluster.rds")

scaled_data <- GetAssayData(data_cluster,assay="RNA",layer="scale.data")
tenscRNA <- SetAssayData(object=tenscRNA,assay="RNA",slot="scale.data",new.data=scaled_data)
VariableFeatures(tenscRNA) <- VariableFeatures(data_cluster)
tenscRNA[["pca"]] <- data_cluster[["pca"]]
tenscRNA[["harmony"]] <- data_cluster[["harmony"]]
tenscRNA[["umap"]] <- data_cluster[["umap"]]
tenscRNA@meta.data <- data_cluster@meta.data

saveRDS(tenscRNA,file="./RDS_file/data_cluster_all.rds")

data <- readRDS("./RDS_file/data_cluster_all.rds")

marker_genes <- list(
  Pericycle = c("GmCH-02G027480", "GmCH-10G262450", "GmCH-20G561760"),
  Xylem = c("GmCH-02G053760", "GmCH-03G077330", "GmCH-18G508470"),
  Phloem = c("GmCH-04G092960", "GmCH-06G162110", "GmCH-12G323770"), 
  Epidermis = c("GmCH-08G229910","GmCH-17G456780"),
  Cortex = c("GmCH-06G154030", "GmCH-12G344250", "GmCH-15G418850"),
  Endodermis = c("GmCH-13G354660", "GmCH-18G503850", "GmCH-06G167460"),
  Infected_cell = c("GmCH-10G284570", "GmCH-02G035860", "GmCH-11G315630"),
  Unintected_cell = c("GmCH-13G354640", "GmCH-07G196560", "GmCH-17G455200"))

marker_genes_in_data <- list()

for (ct in names(marker_genes)){
  markers <- marker_genes[[ct]]
  markers_found <- markers[markers %in% rownames(data)]
  marker_genes_in_data[[ct]] <- markers_found
}

print(marker_genes_in_data)
marker_genes_vector <- unlist(marker_genes_in_data)
dotplot_list <- list()
dotplot_1 <- DotPlot(data, features = marker_genes_vector, cols=c("grey","red"),group.by = "leiden_0.4") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5),
        strip.text.x = element_text(angle = 90,hjust = 1),
        strip.text = element_text(margin=margin(b=2, unit="mm")),strip.placement = 'outlet') + labs(x="",y="")

ggsave("markerAnno_filter_dot_v5_own_0.4.pdf",plot = dotplot_1, width = 12, height = 6)
