#-----------------------------------------------------------------------------------------------
print("Read in the data and load relevant libraries")
outputFolder  = "output/"
scriptsFolder = "R/"
inputFolder   = "data/"

## Load these libraries
library(dplyr);
library(feather)
library(ggplot2)

## Source the violin plot and heatmap functions
source(paste0(scriptsFolder,"violin_functions.R"))
source(paste0(scriptsFolder,"heatmap_functions.R"))

## Read in the table of genes and clusters for each plot
hgncTable = read.csv(paste0(inputFolder,"HGNC_input.csv"))

## Make the plots and save them as pdf files.
for (i in unique(hgncTable$plotNumber)){
  plotTable   = hgncTable[hgncTable$plotNumber==i,]
  eval(parse(text=paste0("clusters = c(",plotTable$clusters[1],")"))) # Read clusters
  genes       = as.character(plotTable$Gene)
  shinyFolder = paste0(inputFolder, plotTable$species[1], "/")
  fileName    = paste0(outputFolder,plotTable$plot,"_",plotTable$species,".pdf")[1]
  width       = length(clusters)/2+2
  height      = length(genes)/2+5
  plotType    = plotTable$plotType[1]
  groupBy     = plotTable$groupBy[1]
  if(plotType=="heatmap"){
    sample_heatmap_plot(shinyFolder, group_by = groupBy, groups=clusters, genes=genes)
	ggsave(fileName, width = width, height = height/3)
    build_legend_plot(shinyFolder, genes=genes)
	ggsave(gsub(".pdf","_legend.pdf",fileName), width = 5, height = 2)
  }
  if(plotType=="violin"){
    genes = genes[length(genes):1]
    group_violin_plot(data_source = shinyFolder, group_by = groupBy, clusters=clusters, genes=genes, sort=TRUE)
	ggsave(fileName, width = width/2, height = height/2)
  }  
}
