# ICPS
Accurate identification of cell type and phenotypic marker genes in single cell transcriptomic data via a data augmentation approach

## Pipeline of ICPS
![image](https://raw.githubusercontent.com/changwn/ICPS/master/fig/fig1.png)

## Installation

```
#install ICPS
install.packages("devtools")
devtools::install_github("changwn/ICPS")
```

## Example


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
