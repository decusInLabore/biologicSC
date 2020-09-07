



###############################################################################
## Scan Seurat Parameters                                                    ##



#' Scan Parameters
#'
#' TBD.
#'
#' @param obj Seurat object
#' @return paramerer list
#' @export

setGeneric(
  name="scanObjParams",
  def=function(
    obj
  ) {

    addReductions <- function(red = "pca", obj, paramList){
      reds <- names(obj@reductions)
      pos <- grep(red, reds)
      if (length(pos) > 0){
        dfTemp <- data.frame(obj@reductions[[red]]@cell.embeddings)
        if (nrow(dfTemp) > 0){
          paramList[[red]] <- names(dfTemp)

        }
      }
      return(paramList)
    }

    tempList <- list()

    tempList[["meta.data"]] <- c(names(obj@meta.data))
    names(tempList[["meta.data"]]) <- gsub("[.]", "_",gsub("meta.data", "", c(names(obj@meta.data))))
    reds <- names(obj@reductions)

    for (i in 1:length(reds)){
      tempList <- addReductions(
        red = reds[i],
        obj = obj,
        paramList = tempList
      )
    }

    allOptions <- unlist(tempList, use.names = F)
    allOptions <- allOptions[allOptions != "all"]
    uPos <- allOptions[grep("UMAP", allOptions)]
    tPos <- allOptions[grep("tSNE", allOptions)]
    clustPos <- allOptions[grep("cluster", allOptions)]
    rmPos <- c(uPos, tPos, clustPos)
    restPos <- allOptions[!(allOptions %in% rmPos)]
    allOptions <- c(uPos, tPos, clustPos, restPos)
    names(allOptions) <- gsub("[.]", "_", allOptions)



    paramList <- list()
    paramList[["x_axis"]] <- c("log10 Expr" = "lg10Expr", allOptions)
    paramList[["y_axis"]] <- c("log10 Expr" = "lg10Expr", allOptions)

    catOptions <- as.vector(NULL, mode = "character")
    for (i in 1:ncol(obj@meta.data)){
      if (length(unique(obj@meta.data[,i])) <= 20){
        catOptions <- c(
          catOptions,
          names(obj@meta.data)[i]
        )
      }
    }

    names(catOptions) <- gsub("[.]", "_", catOptions)
    #names(catOptions) <- gsub("orig_ident", "sampleID", names(catOptions))
    names(catOptions) <- gsub("seurat_clusters", "Cluster", names(catOptions))
    names(catOptions) <- gsub("all", "None", names(catOptions))


    paramList[["splitPlotsBy"]] <- catOptions

    names(catOptions) <- gsub("None", "Black", names(catOptions))
    paramList[["colorPlotsBy"]] <- catOptions

    ## Create color list ##
    sampleColorList <- list()
    for (i in 1:length(paramList[["colorPlotsBy"]])){
      tag <- paramList[["colorPlotsBy"]][i]

      sampleVec <- as.vector(sort(unique(obj@meta.data[,tag])))
      sampleVec <- na.omit(sampleVec)
      sampleVec <- sampleVec[sampleVec != ""]

      if (length(sampleVec) == 1){
        sampleColVec <- "black"
      }  else if (length(sampleVec) == 2){
        l1 <- length(obj@meta.data[obj@meta.data[,tag] == sampleVec[1],tag])
        l2 <- length(obj@meta.data[obj@meta.data[,tag] == sampleVec[2],tag])
        if (l1 > l2){
          sampleColVec <- c("black", "red")
        } else {
          sampleColVec <- c("red", "black")
        }

      } else {
        library(scales)
        sampleColVec <- hue_pal()(length(sampleVec))

      }

      names(sampleColVec) <- sampleVec
      sampleColorList[[names(paramList[["colorPlotsBy"]])[i]]] <-  sampleColVec
    }
    paramList[["sampleColorList"]] <- sampleColorList
    ## Done creating color list

    return(paramList)

  }

)

###############################################################################

###############################################################################
## Add to Seurat metadata                                                    ##
setGeneric(
  name="addDf2seuratMetaData",
  def=function(obj, dfAdd) {
    #print(paste0("Dims before addition: ", dim(obj@meta.data)))

    for (i in 1:ncol(dfAdd)){
      addVec <- as.vector(dfAdd[,i])
      names(addVec) <- row.names(dfAdd)
      colName <- as.vector(names(dfAdd)[i])
      obj <- Seurat::AddMetaData(
        object = obj,
        metadata = addVec,
        colName
      )
    }

    #print(paste0("Dims after addition: ", dim(obj@meta.data)))
    #print(paste0("Meta data column names: ", paste(names(obj@meta.data), collapse = ", ")))
    return(obj)
  }
)

## Done adding to Seurat metadata                                            ##
###############################################################################

###############################################################################

setGeneric(
  name="createDfCoord",
  def=function(
    obj,
    params = NULL
  ) {
    if (is.null(params)){
      params <- scanObjParams(obj)
    }

    ## Add reductions ##
    reds <- names(obj@reductions)

    for (i in 1:length(reds)){
      dfAdd <- data.frame(obj@reductions[[reds[i]]]@cell.embeddings)
      if (nrow(dfAdd) > 0){
        obj <- addDf2seuratMetaData(obj = obj, dfAdd = dfAdd)
      }
    }

    ##
    dfdbTable <- obj@meta.data
    dfdbTable[["cellID"]] <- row.names(dfdbTable)
    pos <- grep("sampleID", names(dfdbTable))
    pos2 <- grep("orig.ident", names(dfdbTable))
    if (length(pos) == 0 | length(pos2) == 1){
      dfdbTable[["sampleID"]] <- dfdbTable[["orig.ident"]]
    } else {
      dfdbTable[["sampleID"]] <- "sampleID_TBD"
    }

    return(dfdbTable)

  }
)

###############################################################################

###############################################################################

setGeneric(
  name="createDfExpr",
  def=function(
    obj,
    assay = "RNA",
    #slot = "data",
    geneSel = NULL
  ) {
    Seurat::DefaultAssay(obj) <- assay
    dfExpr <- data.frame(obj[[assay]]@data)
    dfExpr[["gene"]] <- row.names(dfExpr)

    if (!is.null(geneSel)){
      dfExpr <- dfExpr[dfExpr$gene %in% geneSel, ]
    }

    dfExpr <- tidyr::gather(
      dfExpr,
      condition,
      expr, 1:(ncol(dfExpr)-1),
      factor_key=TRUE
    )
    dfExpr <- dfExpr[dfExpr$expr != 0,]
    names(dfExpr) <- gsub("condition", "cellID", names(dfExpr))
    names(dfExpr) <- gsub("expr", "lg10Expr", names(dfExpr))
    return(dfExpr)
  }
)



###############################################################################

###############################################################################
#' Make Seurat Viewer
#'
#' TBD2.
#'
#' @param obj Seurat object
#' @return app
#' @export

setGeneric(
  name="seurat2viewer",
  def=function(
    obj,
    assay = "RNA",
    #slot = "data",
    geneSel = NULL,
    params = NULL,
    projectName = "test"
  ) {
    dfCoord <- createDfCoord(obj = testObj, params = params)
    dfExpr <- createDfExpr(obj = testObj, assay = "RNA")

    dfIDTable <- dfExpr
    dfIDTable[["gene_id"]] <- 0
    dfIDTable <- unique(dfIDTable[,c("gene", "gene_id")])
    dfIDTable <- dfIDTable[order(dfIDTable$gene, decreasing = F), ]
    dfIDTable[["gene_id"]] <- 1:nrow(dfIDTable)


    ###############################################################################
    ## Add percentage expressed genes                                            ##


    top30Var <- head(
      x = VariableFeatures(object = testObj),
      30
    )

    my_genes <- rownames(x = testObj@assays$RNA)
    exp <- FetchData(testObj, my_genes)
    ExprMatrix <- round(as.matrix(colMeans(exp  > 0)) *100,1)
    colnames(ExprMatrix)[1] <- "count_cut_off"
    dfExprMatrix <- data.frame(ExprMatrix)
    dfExprMatrix[["gene"]] <- row.names(dfExprMatrix)
    dfExprMatrix <- dfExprMatrix[dfExprMatrix$gene %in% top30Var, ]
    dfExprMatrix <- dfExprMatrix[order(dfExprMatrix$count_cut_off, decreasing = T),]
    geneDefault = as.vector(dfExprMatrix[1,"gene"])

    ############
    ## Create database
    projectDir <- paste0(projectName)
    dataDir <- paste0(projectDir, "/data")
    paramDir <- paste0(projectDir, "/parameters")

    connectDir <- paste0(projectDir, "/connect")
    projectDB <- paste0(projectName, "_DB")

    coordTb <- paste0(projectName, "_meta_data")
    exprTb <- paste0(projectName, "_gene_expr_tb")
    geneTb <- paste0(projectName, "_geneID_tb")

    if (!dir.exists(projectDir)){
      dir.create(projectDir)
    }

    if (!dir.exists(dataDir)){
      dir.create(dataDir)
    }

    if (!dir.exists(connectDir)){
      dir.create(connectDir)
    }

    if (!dir.exists(paramDir)){
      dir.create(paramDir)
    }

    setwd(dataDir)
    conn <- RSQLite::dbConnect(RSQLite::SQLite(), projectDB)
    RSQLite::dbWriteTable(conn, coordTb, dfCoord, overwrite =T)
    RSQLite::dbWriteTable(conn, exprTb, dfExpr, overwrite =T)
    RSQLite::dbWriteTable(conn, geneTb, dfIDTable, overwrite =T)
    RSQLite::dbDisconnect(conn)


    dfID <- data.frame(
      type = "RSQLite",
      url = "",
      id = "",
      id2 = "",
      db = projectDB,
      coordTb,
      exprTb,
      geneTb,
      default = geneDefault
    )
    setwd("../connect")
    write.table(
      dfID,
      "db.txt",
      sep = "\t",
      row.names = F
    )

    setwd("../parameters")
    yamlList <- list(
      "XYsel" = params[["x_axis"]],
      "allColorOptions" = params[["colorPlotsBy"]],
      "splitOptions" = params[["splitPlotsBy"]],
      "sampleColorList" = params[["sampleColorList"]]
    )

    FN <- paste0("parameters.yaml")
    yaml::write_yaml(yaml::as.yaml(yamlList), FN, fileEncoding = "UTF-8")
    setwd("..")

    file.copy(system.file("ui.r", package = "biologicSC"), ".")
    file.copy(system.file("server.r", package = "biologicSC"), ".")

  }
)

###############################################################################


