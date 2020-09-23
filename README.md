# biologicSC

```
devtools::install_github("decusinlabore/biologicSC")
library(biologicSC)
library(Seurat)

testObj <- pbmc_small
all.genes <- rownames(testObj)
testObj <- ScaleData(testObj, features = all.genes)
testObj <- RunPCA(testObj, npcs = 3)
testObj <- RunUMAP(testObj, dims = 1:3)

testObj <- FindNeighbors(testObj, dims = 1:3)
testObj <- FindClusters(testObj, resolution = 0.5)

testObj@meta.data[["all"]] <- "all"

params <- scanObjParams(testObj)

seurat2viewer(
    obj = testObj,
    assay = "RNA",
    #slot = "data",
    geneSel = NULL,
    params = params,
    projectName = "pbmc_small_app"
) 

setwd("..")
library(shiny)
runApp("pbmc_small_app")
```
