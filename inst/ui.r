###############################################################################
##                                                                           ##

##                                                                           ##
###############################################################################

###############################################################################
##                                                                           ##
dfkey <- read.delim("connect/db.txt", header = T, sep="\t", stringsAsFactors = F)

library(shiny)
library(ggplot2)
library(yaml)

pos <- grep("type", names(dfkey))
if (length(pos) > 0 & dfkey$type == "RSQLite"){
  library(RSQLite)
  mode = "SQLite"
} else {
  library(RMySQL)
  host <- as.character(dfkey$url)
  user <- as.character(dfkey$id)
  DBpd <- as.character(dfkey$id2)
  mode = "SQL"
}

library(DT)

geneDefault = as.character(dfkey$default)
dbname <- as.character(dfkey$db)
coordinateTbName <- as.character(dfkey$coordTb)
exprTbName <- as.character(dfkey$exprTb)
geneID_TbName <- as.character(dfkey$geneTb)

##                                                                           ##
###############################################################################


###############################################################################
##                                                                           ##


##                                                                           ##
###############################################################################


###############################################################################
##                                                                           ##
oldw <- getOption("warn")
options(warn = -1)
if (mode == "SQLite"){
  setwd("data")
  dbDB <- dbConnect(RSQLite::SQLite(), dbname=dbname)
  setwd("..")
} else {
  dbDB <- dbConnect(MySQL(), user = user, password = DBpd, host = host, dbname=dbname)
}
query <- paste0("SELECT DISTINCT * FROM ", coordinateTbName)
dfCoordSel <- dbGetQuery(dbDB, query)
dbDisconnect(dbDB)

dfCoordSel[["all"]] <- "all"
##                                                                           ##
###############################################################################


###############################################################################
## Get XYsel from yaml if possible                                           ##

if (file.exists("parameters/parameters.yaml")){
  params <- yaml.load(
    read_yaml(
      "parameters/parameters.yaml",
      fileEncoding = "UTF-8"
    )
  )
} else {
  params <- list(
    "XYsel" = c(
      names(dfCoordSel)[grep("UMAP", names(dfCoordSel))],
      names(dfCoordSel)[grep("tSNE", names(dfCoordSel))],
      names(dfCoordSel)[grep("seurat_clusters", names(dfCoordSel))],
      names(dfCoordSel)[grep("sub_clusters", names(dfCoordSel))],
      names(dfCoordSel)[grep("PC", names(dfCoordSel))],
      names(dfCoordSel)[grep("DM_Pseudotime", names(dfCoordSel))],
      names(dfCoordSel)[grep("DF_Classification", names(dfCoordSel))],
      names(dfCoordSel)[grep("nCount", names(dfCoordSel))],
      names(dfCoordSel)[grep("nFeatures", names(dfCoordSel))],
      names(dfCoordSel)[grep("percent", names(dfCoordSel))],
      names(dfCoordSel)[grep("Con_Prad_AZ", names(dfCoordSel))],
      names(dfCoordSel)[grep("Norm_Hyp", names(dfCoordSel))],
      names(dfCoordSel)[grep("Gender", names(dfCoordSel))],
      names(dfCoordSel)[grep("CellFromTumor", names(dfCoordSel))],
      names(dfCoordSel)[grep("Patient", names(dfCoordSel))],
      names(dfCoordSel)[grep("Region", names(dfCoordSel))],
      names(dfCoordSel)[grep("Article_Cell_Type", names(dfCoordSel))]
    ),
    "allColorOptions" = c(
      #"Log10 Expresson" = "lg10Expr",
      "DM Pseudotime"  = names(dfCoordSel)[grep("DM_Pseudotime", names(dfCoordSel))],
      "SampleID" = names(dfCoordSel)[grep("sampleID", names(dfCoordSel))],
      "Gender" = names(dfCoordSel)[grep("Gender", names(dfCoordSel))],
      "Patient" = names(dfCoordSel)[grep("Patient", names(dfCoordSel))],
      "Doublet Classification" = names(dfCoordSel)[grep("DF_Classification", names(dfCoordSel))] ,
      "Cluster" = names(dfCoordSel)[grep("seurat_clusters", names(dfCoordSel))] ,
      "nCount_RNA" = names(dfCoordSel)[grep("nCount_RNA", names(dfCoordSel))],
      "nFeature_RNA" = names(dfCoordSel)[grep("nFeature_RNA", names(dfCoordSel))],
      "percent_mt" = names(dfCoordSel)[grep("percent_mt", names(dfCoordSel))],
      "S Phase Score" = names(dfCoordSel)[grep("S_Score", names(dfCoordSel))],
      "G2M Score" = names(dfCoordSel)[grep("G2M_Score", names(dfCoordSel))],
      "Cell Cycle Phase" = names(dfCoordSel)[grep("Phase", names(dfCoordSel))],
      "Uniform" = names(dfCoordSel)[grep("all", names(dfCoordSel))]
    ),
    "splitOptions" = c(
      #"Log10 Expresson" = "lg10Expr",
      "DM Pseudotime"  = names(dfCoordSel)[grep("DM_Pseudotime", names(dfCoordSel))],
      "SampleID" = names(dfCoordSel)[grep("sampleID", names(dfCoordSel))],
      "Gender" = names(dfCoordSel)[grep("Gender", names(dfCoordSel))],
      "Patient" = names(dfCoordSel)[grep("Patient", names(dfCoordSel))],
      "Doublet Classification" = names(dfCoordSel)[grep("DF_Classification", names(dfCoordSel))] ,
      "Cluster" = names(dfCoordSel)[grep("seurat_clusters", names(dfCoordSel))] ,
      "nCount_RNA" = names(dfCoordSel)[grep("nCount_RNA", names(dfCoordSel))],
      "nFeature_RNA" = names(dfCoordSel)[grep("nFeature_RNA", names(dfCoordSel))],
      "percent_mt" = names(dfCoordSel)[grep("percent_mt", names(dfCoordSel))],
      "S Phase Score" = names(dfCoordSel)[grep("S_Score", names(dfCoordSel))],
      "G2M Score" = names(dfCoordSel)[grep("G2M_Score", names(dfCoordSel))],
      "Cell Cycle Phase" = names(dfCoordSel)[grep("Phase", names(dfCoordSel))],
      "Uniform" = names(dfCoordSel)[grep("all", names(dfCoordSel))]
    )

  )
}

##                                                                           ##
###############################################################################


Xchoices <- c("Log10 Expression" = "lg10Expr",
    params[["allColorOptions"]],
    params[["XYsel"]]
  )


if (length(grep("UMAP_1", Xchoices)) == 1){
  Xsel <- Xchoices[grep("UMAP_1", Xchoices)]
} else {
  Xsel <- Xchoices[1]
}

Ychoices <- c("Log10 Expression" = "lg10Expr",
    params[["XYsel"]]
)


if (length(grep("UMAP_2", Ychoices)) == 1){
  Ysel <- Ychoices[grep("UMAP_2", Ychoices)]
} else {
  Ysel <- Ychoices[2]
}

splitChoices <- params[["splitOptions"]]
if (length(grep("sampleID", splitChoices)) == 1){
  splitSel <- splitChoices[grep("sampleID", splitChoices)]
} else {
  splitSel <- splitChoices[1]
}

colorChoices <- c("Log10 Expression" = "lg10Expr",
    params[["allColorOptions"]]
  )

if (length(grep("^lg10Expr$", colorChoices)) == 1){
  colSel <- colorChoices[grep("^lg10Expr$", colorChoices)]
} else {
  colSel <- colorChoices[1]
}

spectralCols <- c(Darkred = "#D53E4F",
                  Red = "#F46D43",
                  Orange = "#FDAE61" ,
                  Lightorange = "#FEE08B",
                  Yellow = "#FFFFBF",
                  Lightgreen = "#E6F598",
                  Green = "#ABDDA4",
                  Darkgreen = "#66C2A5",
                  Blue =  "#3288BD")

conditionVec <- unique(sort(dfCoordSel$sampleID))

Nsamples <- length(conditionVec)



##                                                                           ##
###############################################################################


###############################################################################
##                                                                           ##

oldw <- getOption("wafrn")
options(warn = -1)

if (mode == "SQLite"){
  setwd("data")
  dbDB <- dbConnect(RSQLite::SQLite(), dbname=dbname)
  setwd("..")
} else {
  dbDB <- dbConnect(MySQL(), user = user, password = DBpd, host = host, dbname=dbname)
}
query <- paste0("SELECT DISTINCT gene FROM ", geneID_TbName)
allGenes <- as.vector(dbGetQuery(dbDB, query)[,"gene"])
dbDisconnect(dbDB)

##                                                                           ##
###############################################################################



###############################################################################
##                                                                           ##

shinyUI(fluidPage(
    navbarPage(


        "bioLOGIC SC",

        tabPanel("FeatureView"),
        tags$head(

          tags$style(type = 'text/css',
                     HTML('.navbar { background-color: #42972050;}
                          .navbar-default .navbar-brand{color: white;}
                          .tab-panel{ background-color: red; color: white}
                          .navbar-default .navbar-nav > .active > a,
                           .navbar-default .navbar-nav > .active > a:focus,
                           .navbar-default .navbar-nav > .active > a:hover {
                                color: #555;
                                background-color: #42972050;
                            }')
          ),
          tags$script(HTML("var header = $('.navbar > .container-fluid');
header.append('<div style=\"float:left\"><ahref=\"URL\"><img src=\"assets/images/logo.ico\" alt=\"alt\" style=\"float:right;width:33px;height:41px;padding-top:10px;\"> </a>`</div>');
    console.log(header)")
          ),
          tags$link(rel="shortcut icon", href="assets/images/logo.ico")
        )

    ),
    #titlePanel("FeatureView"),
    sidebarLayout(
        sidebarPanel(
            tags$style(".well {background-color:#42972050;}"),
            helpText("This application puts more than 100 million data points at your fingertips. Please be patient at the beginning when the application starts."),

            selectizeInput("gene",
                           label = "Gene or Category Selection",
                           choices = NULL, #c(as.vector(sort(unique(allGenes)))),
                           selected = geneDefault,
                           options = list(maxOptions = 50)) ,

            selectInput("x_axis",
                        label = "Choose an X-axis",
                        choices =Xchoices,
                        selected = Xsel),
            selectInput("y_axis",
                        label = "Choose an Y-axis",
                        choices =Ychoices,
                        selected = Ysel),

            selectInput("splitByColumn",
                        label = "Split Plots By",
                        choices = splitChoices,
                        selected = splitSel),

            selectInput("colorBy",
                        label = "Color Plots By",
                        choices = colorChoices,
                        selected = colSel),


            selectInput("dotcolor",
                        label = "Choose dot colorscale",
                        choices =c("Darkblue" = "darkblue", spectralCols),
                        selected = "darkblue"),

            selectInput("lowColor",
                        label = "Choose low colorscale",
                        choices =c("Grey" = "#D3D3D3", "White" = "white", spectralCols),
                        selected = "#D3D3D3"),

            selectInput("background",
                        label = "Select Background",
                        choices =c("Grey" = "grey", "White" = "white","Minimal" = "minimal", "Plain" =  "plain"),
                        selected = "white"),




            radioButtons("dotsize", label = "Choose a Dotsize", choices = c("0.1","0.5","1","2"), selected = "1",
                         inline = FALSE, width = NULL, choiceNames = c("0.1","0.5","1","2"),
                         choiceValues = c("0.1","0.5","1","2")),

            # selectInput("dotsize",
            #             label = "Choose an Dotsize",
            #             choices =c("0.1","0.5","1","2"),
            #             selected = "1"),

            #downloadButton('downloadPlot', "Download Plot"),



            #downloadButton("downloadData", "Download Data"),

            # checkboxGroupInput("selected_sample",
            #                    label = "Select Column",
            #                    choices = seq_along(mtcars))


        ),
        mainPanel(

                            fluidRow(
                                column(12,
                                       uiOutput("multi_plot_ui")
                                )
                            ),
                            fluidRow(
                              column(12,
                                     textOutput("dev_text")
                              )
                            )


        )


    )
))

##                                                                           ##
###############################################################################


