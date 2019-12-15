# ICPS
Accurate **Identification** of **Cell** type and **Phenotypic** marker genes in **Single** cell transcriptomic data via a data augmentation approach

## Pipeline of ICPS
![image](https://raw.githubusercontent.com/changwn/ICPS/master/fig/fig1.png)

## Installation

```
#install ICPS
install.packages("devtools")
devtools::install_github("changwn/ICPS")
```

## Example
```
#-------------------------------
#
#  seurat clustering, plot tSNE
#
#-------------------------------

#---seurat
MAT <- scData
library(Seurat)
MAT_seurat<-CreateSeuratObject(counts = MAT, project = "MAT", min.cells = 3, min.features = 200)
MAT_seurat[["percent.mt"]] <- PercentageFeatureSet(MAT_seurat, pattern = "^mt-")
#VlnPlot(MAT_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MAT_seurat <- NormalizeData(MAT_seurat)
MAT_seurat<-FindVariableFeatures(MAT_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MAT_seurat)
MAT_seurat <- ScaleData(MAT_seurat, features = all.genes)
MAT_seurat <- RunPCA(MAT_seurat, features = VariableFeatures(object = MAT_seurat))
#MAT_seurat <- RunPCA(MAT_seurat, features = marker_gene) #use default seurat, no marker gene
MAT_seurat <- JackStraw(MAT_seurat, num.replicate = 100)
MAT_seurat <- ScoreJackStraw(MAT_seurat, dims = 1:20)
#JackStrawPlot(MAT_seurat, dims = 1:20)
#ElbowPlot(MAT_seurat)
MAT_seurat <- FindNeighbors(MAT_seurat, dims = 1:20)
MAT_seurat <- FindClusters(MAT_seurat, resolution = 0.5)
MAT_seurat <- RunTSNE(MAT_seurat, dims = 1:20)
DimPlot(MAT_seurat, reduction = "tsne", label=T)
MAT_ident<-as.vector(MAT_seurat@meta.data[,6])
names(MAT_ident) <- rownames(MAT_seurat@meta.data)

library(cluster)
tsne_position <- MAT_seurat@reductions$tsne@cell.embeddings
test_dist<-dist(tsne_position,method = "canberra")
cluster_label <- strtoi(MAT_ident, base = 0L)
sil<-silhouette(cluster_label,test_dist)
sum(sil[,3])

#----------------------
#
# ICPS: simulate using seurat cluster label, run ictd, clustering using ctes3 marker
#
#----------------------

pseudo_BULK <- Bulk_Simu_no_coinfil_v1(scData, MAT_ident, cellNumber=10000, sampleNumber=50)
data_matrix_simu <- pseudo_BULK[[1]]
tProp <- pseudo_BULK[[2]]

library(ICTD)
CTES3 <- ICTD(data_matrix_simu)

MAT <- scData
ictd_genelist <- extract_list(CTES3)
#ictd_genelist <- extract_list(tg_R1_lists)
library(Seurat)
MAT_seurat_mtct<-CreateSeuratObject(counts = MAT, project = "MAT", min.cells = 3, min.features = 200)
MAT_seurat_mtct[["percent.mt"]] <- PercentageFeatureSet(MAT_seurat_mtct, pattern = "^mt-")
#VlnPlot(MAT_seurat_mtct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MAT_seurat_mtct <- NormalizeData(MAT_seurat_mtct)
MAT_seurat_mtct<-FindVariableFeatures(MAT_seurat_mtct, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(MAT_seurat_mtct)
MAT_seurat_mtct <- ScaleData(MAT_seurat_mtct, features = all.genes)
#MAT_seurat <- RunPCA(MAT_seurat, features = VariableFeatures(object = MAT_seurat))
MAT_seurat_mtct <- RunPCA(MAT_seurat_mtct, features = ictd_genelist)
MAT_seurat_mtct <- JackStraw(MAT_seurat_mtct, num.replicate = 100)
MAT_seurat_mtct <- ScoreJackStraw(MAT_seurat_mtct, dims = 1:20)
#JackStrawPlot(MAT_seurat_mtct, dims = 1:20)
#ElbowPlot(MAT_seurat_mtct)
MAT_seurat_mtct <- FindNeighbors(MAT_seurat_mtct, dims = 1:20)
MAT_seurat_mtct <- FindClusters(MAT_seurat_mtct, resolution = 0.5)
MAT_seurat_mtct <- RunTSNE(MAT_seurat_mtct, dims = 1:20)
DimPlot(MAT_seurat_mtct, reduction = "tsne", label=T)

MAT_ident_mtct<-as.vector(MAT_seurat_mtct@meta.data[,6])
names(MAT_ident_mtct) <- rownames(MAT_seurat_mtct@meta.data)
```


## Instruction
Single cell expression profiles exhibit variations of cell type specific genes among different cell types, as well as general variations that contribute to dynamic physiological state changes. We developed ICPS as a biologically informative single cell clustering approach that teases out the local low rank structures in single cell augmented data, aiming to cluster single cells using only genes with strong local co-expression patterns, that are signs of cell type specificity. We demonstrated that the augmented data could more robustly extract informative genes that could be used to characterize variations caused by cell type identify changes, rather than transient state changes. ICPS can be implemented with state-of-arts cell clustering methods via iterative data augmentation, cell type marker identification and cell clustering. ICPS can also be implemented with pre-identified knowledge of cell type specific markers. We have derived comprehensive sets of genes specifically expressed by 11 immune and stromal cell types in the blood, inflammatory and cancer tissues and 9 cell types in the central nervous system of human and mouse. We validated the implementation of ICPS with five state-of-the-art cell clustering methods to increase the accuracy of cell type inference on multiple CITE-seq, scRNA-seq and snRNA-seq data. Application of ICPS also enables the reconstruction of cell lineage and cell alignment between different data sets. 

## Questions & Problems

If you have any questions or problems when using ICPS, please feel free to open a new issue [here](https://github.com/changwn/ICPS/issues). We will fix the new issue ASAP.  You can also email the maintainers and authors below.

## Contact Information

- [Wennan Chang](https://zcslab.github.io/people/wennan/)
(wnchang@iu.edu)

PhD candidate at BDR group, Indiana University School of Medicine

- [Chi Zhang](https://medicine.iu.edu/departments/genetics/faculty/27057/zhang-chi/)
(czhang87@iu.edu)

Assistant Professor

Department of Medical & Molecular Genetics, Indiana University School of Medicine
