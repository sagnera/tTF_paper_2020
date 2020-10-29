# Single cell RNAseq analysis described in "Temporal patterning of the central nervous system by a shared transcription factor code" by Sagner et al, 2020 ----

# Script by Andreas Sagner (andreas.sagner@manchester.ac.uk)
# 29 Oct 2020

# The following datasets need to be separately downloaded or generated and placed in the input folder:

# For analysis of spinal cord scRNAseq data, the m_neural.rds file from Delile et al. 2019 was used. 
# The file was generated as described in https://github.com/juliendelile/MouseSpinalCordAtlas.

# For analysis of the temporal TFs in the fore-, mid- and hindbrain, please download the dev_all.loom file 
# (app. 19 GB) from Manno et al. 2020 from http://mousebrain.org/downloads.html.

# For comparison of temporal TF expression between in-vivo and in-vitro, RNAseq data from Rayon et al. 2020 was used. 
# The file GSE140748_expression_matrix.mouse.abundance.tsv can be obtained from the GEO database 
# (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140748).

# For analysis of the temporal TFs in the retina, please download 10X scRNAseq data of the developing retina 
# (count matrix, cellular phenotype data and feature data) from https://github.com/gofflab/developing_mouse_retina_scRNASeq (Clark et al. 2019).

# Analysis of temporal genes in the spinal cord ----

rm(list=ls()) # clears environment 

# check if all necessary packages are installed, or install them if not

# * Packages ----
dir <- dirname(rstudioapi::getSourceEditorContext()$path)

if (!dir.exists(paste0(dir, 'output'))) {dir.create(paste(dir, 'output'))}

setwd(dir = paste0(dir, '/output/'))

packages <- c("Biobase", "dplyr", "plyr", "scater", "heatmap.plus", "Seurat", "dendextend", "tibble", "ggplot2", "heatmap.plus")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

library(Biobase)
library(plyr)
library(ggplot2)
library(scater)
library(dplyr)
library(Seurat)
library(dendextend)
library(tibble)

# * Functions ----

# functions for converting ensemblIDs into real gene names and vice versa

convert.to.ensemblID <- function(genes) {
  return(unlist(lapply(genes, function(x) {return(rownames(Biobase::fData(eset))[which(Biobase::fData(eset)[, "external_gene_name"] == x)])})))
}

convert.to.realname <- function(ensemblIDs, eset) {
  return(unlist(lapply(ensemblIDs, function(x) {return(Biobase::fData(eset)$external_gene_name[which(rownames(Biobase::fData(eset)) == x)])})))
}

## This function plots the expression dynamics of a gene from spinal cord scRNAseq data (Delile et al. 2019, Development)

## threshold defines how many cells in a domain must express a gene to be considered
## print.file = TRUE means a pdf will be generated, FALSE just sends the output to the RStudio console
## include.FP = TRUE includes floor plate and roof plate
## celltype can be "both", "Neuron" or "Progenitor"

plot.mean.dynamics <- function(gene, threshold = 3, print.file = FALSE, include.FP = TRUE, celltype = "both") {
  geneID <- rownames(Biobase::fData(eset))[which(Biobase::fData(eset)[, "external_gene_name"] == gene)]
  
  if(celltype == "both") {
    exprs.values <- Biobase::exprs(eset)[geneID, ]
    df <- cbind(Biobase::pData(eset)[, c("timepoint", "DV", "Type_step1", "Type_step2")], exprs.values) 
  } else {
    celllist <- rownames(Biobase::pData(eset))[which(Biobase::pData(eset)$Type_step1 == celltype)]
    exprs.values <- Biobase::exprs(eset)[geneID, celllist]  
    df <- cbind(Biobase::pData(eset)[celllist, c("timepoint", "DV", "Type_step1", "Type_step2")], exprs.values) 
  }
  
  ### cellIDs with less then threshold counts get removed
  
  df.count <- plyr::count(df, c("timepoint", "DV", "Type_step2")) %>%
    dplyr::mutate(num_cells_expressing = df %>% 
                    dplyr::group_by(timepoint, DV, Type_step2) %>%
                    dplyr::summarise(count2 = length(exprs.values[exprs.values > 0])) %>%
                    dplyr::pull(count2)) %>%
    dplyr::mutate(percentage.expressing = num_cells_expressing / freq)
  
  df.agg <- ddply(df, .(timepoint, DV, Type_step2), summarize, mean=mean(exprs.values), "DV" = unique(DV), "class" = unique(Type_step2)) %>%
    dplyr::mutate(norm = mean / max(mean)) %>%
    tibble::add_column(num_cells = df.count$freq, num_cells_expressing = df.count$num_cells_expressing, perc = as.numeric(df.count$percentage.expressing)) %>%
    dplyr::filter(num_cells > threshold) 
  
  if(include.FP == FALSE){
    df.agg <- df.agg %>%
      dplyr::filter(DV != 1) %>%
      dplyr::filter(DV != 13)
  }
  
  if(celltype %in% c("both", "Progenitor")) {
    order.cell.types <- c("FP", "p3", "pMN", "p2", "p1", "p0", "dp6", "dp5", "dp4", "dp3", "dp2", "dp1", "RP")
    df.agg$cell.types <- unlist(lapply(df.agg$DV, function(x) {return(order.cell.types[x])}))
  } else {
    order.cell.types <- c("FP", "V3", "MN", "V2a", "V2b", "V1", "V0", "dl6", "dl5", "dl4", "dl3", "dl2", "dl1", "RP")
    df.agg$cell.types <- df.agg$Type_step2
    df.agg <- dplyr::filter(df.agg, cell.types != "dl6")
  }
  
  df.agg$cell.types <- factor(df.agg$cell.types, levels = order.cell.types) 
  
  df.agg$class <- unlist(lapply(df.agg$class, function(x) {
    if (x %in% c("FP", "p3", "pMN", "p2", "p1", "p0", "dp6", "dp5", "dp4", "dp3", "dp2", "dp1", "RP")) {
      return("Progenitor")
    } else {
      return("Neuron")  
    }
  }))
  
  gg <- ggplot(data = df.agg, aes(x = factor(timepoint), y = cell.types)) +
    geom_count(aes(size = norm, color = perc)) +
    scale_size_area(max_size = 20) +
    scale_color_gradientn(colours = c("blue", "yellow"), limits = c(0,1)) +
    theme_classic() +
    ylab("DV position") +
    xlab("embryonic days") +
    ggtitle(gene) +
    theme(axis.text.y = element_text(size = 24),
          axis.title.y = element_text(size = 28),
          axis.text.x = element_text(size = 24, angle = 45, hjust = 1, vjust = 1),
          axis.title.x = element_text(size = 28),
          legend.title = element_text(size = 20),
          legend.title.align=0.5,
          legend.text = element_text(size = 20),
          strip.text = element_text(size = 24, face = "bold"),
          plot.title = element_text(face = "bold", size = 48))
  
  if(celltype == "both") {
    gg <- gg + facet_wrap(~ unlist(list(df.agg$class)))
  }
  
  if(print.file == TRUE) {
    pdf(paste0("Gene_expression_dynamics_", gene, ".pdf"), width =10, height = 10, useDingbats = FALSE)
    print(gg)
    dev.off()
  }
  
  return(gg)
}

# generate multiple plots in the same PDF

plot.dynamics <- function(genes, title, celltype = "both") {
  celltype = celltype
  gene.plots <- lapply(genes[which(genes %in% Biobase::fData(eset)$external_gene_name)], FUN = plot.mean.dynamics, celltype = celltype)
  
  n <- ceiling(sqrt(length(genes)))
  
  if(celltype == "both") {
    k = 10
  } else {
    k = 6
  }
  
  pdf(paste0(title), width = n*k, height = ceiling(length(genes) / n) * 10)
  print(multiplot(plotlist = gene.plots, cols = n))
  dev.off()
}

# * Load 10X in-vivo data and TF list ----

# list of mouse TFs downloaded from AnimalTFDB3.0 (http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/tf_summary?species=Mus_musculus)
TF.lst <- read.delim(paste0(dir, '/input/Mus_musculus_TF.txt'), header = TRUE)$Symbol 

eset = readRDS(paste0(dir, '/input/m_neural.rds'))
eset = eset$expressionSet

rownames(Biobase::pData(eset)) <- gsub("-", ".", rownames(Biobase::pData(eset)))
colnames(Biobase::exprs(eset)) <- gsub("-", ".", colnames(Biobase::exprs(eset)))

# * Load data into the Seurat package ----
mat <- Biobase::exprs(eset)
rownames(mat) <- Biobase::fData(eset)[,'external_gene_name']

seurat <- CreateSeuratObject(counts = mat, meta.data = Biobase::pData(eset), project = "MouseSpinalCordAtlas")

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")

seurat <- seurat %>% 
  subset(subset = nFeature_RNA > 600 & nFeature_RNA < 6000 & percent.mt < 6) %>%
  NormalizeData(verbose=FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30)

# * Differential gene expression analysis of progenitor domains ----

## define progenitor domains, we exclude dI6 due to low number of cells
domains.p <- c("dp1", "dp2", "dp3", "dp4", "dp5", "p0", "p1", "p2", "pMN", "p3") 

# runs subclustering on progenitor domains
markers.age <- lapply(domains.p, function(x, threshold =  0.001) {
  
  print(paste0("Subclustering domain ", x))
  
  celllist <- rownames(Biobase::pData(eset))[which(Biobase::pData(eset)$Type_step2 == x)]
  
  seurat.sub <- subset(seurat, cells = celllist) %>%
    FindVariableFeatures(selection.method = "vst", verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 30, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:30)
  
  Idents(seurat.sub) <- "timepoint"
  age.markers <- FindAllMarkers(seurat.sub, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
  
  return(age.markers)
})

# limit to to genes identified as differentially expressed in multiple domains

num.domains = 7 ### min 8 domains

age.markers <- table(unlist(lapply(markers.age, function(x) {return(unique(x[, "gene"]))}))) %>%
  as.data.frame() %>%
  dplyr::filter(Freq > num.domains) %>% 
  dplyr::arrange(desc(Freq)) 

age.TFs <- intersect(age.markers$Var1, TF.lst)

print(paste0(length(age.markers$Var1), " genes detected!"))
print(paste0(length(age.TFs), " TFs detected!"))

# * Calculate heatmaps ----

seurat.pt <- subset(seurat, subset = Type_step1 == "Progenitor")
Idents(seurat.pt) <- "timepoint"

# fit gene expression dynamics for heatmaps
pt.mat <- as.matrix(GetAssayData(object = seurat.pt, slot = 'counts'))
colnames(pt.mat) <- gsub("-", ".", colnames(pt.mat))

# log and recenter dataset
pt.mat_log =  log(0.000001+pt.mat)
pt.mat_zscored =  t(scale(t(pt.mat_log), center=T, scale=T))

# generate dataframe of in-vivo genes for plotting with ggplot
pt.mat.TFs <- data.frame(pt.mat_zscored[as.character(age.markers$Var1), ]) %>%
  tibble::rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  tidyr::gather(key = "cellname", value = "readcount", -gene) %>%
  dplyr::mutate("timepoint" = Biobase::pData(eset)[.$cellname, "timepoint"]) %>%
  dplyr::group_by(gene, timepoint) %>%
  dplyr::summarise(mean = mean(readcount)) %>%
  dplyr::ungroup() %>%
  tidyr::spread(gene, value = "mean") %>%
  tibble::column_to_rownames("timepoint") %>%
  t() %>% as.matrix()

# generate hierarchical clustering

correlation_dist =  as.dist((1 - cor(t(pt.mat.TFs), method="pearson")))
hc = stats::hclust(correlation_dist, method="ward.D2")

clusters = cutree(hc, k=2)
table(clusters)

order.TFs.clusters <- intersect(rownames(pt.mat.TFs)[dendextend::order.hclust(hc)], TF.lst)

# * Plot heatmaps ----

# plot in-vivo heatmap (change "height" and "cexRow" to improve readability)
pdf("Figure_5A_TFs.pdf", height=5, width=3)
heatmap.plus::heatmap.plus(
  pt.mat.TFs[order.TFs.clusters, ],
  scale='row',
  Colv = NA,
  Rowv = NA, 
  col=colorRampPalette(c("#191d73", "white", "#ed7901"))(n = 1000),
  labCol = colnames(pt.mat.TFs),
  cexRow=.8# font size
) 
graphics.off()

lab.row <- unlist(lapply(rownames(pt.mat.TFs), function(x) {
  if(x %in% TF.lst) {
    print(x)
    return(x)
  } else {
    return("")  
  }  
})) 

pdf("Figure_5A_all_genes.pdf", height=10, width=3)
heatmap.plus::heatmap.plus(
  pt.mat.TFs,
  scale='row',
  Colv = NA,
  Rowv = as.dendrogram(hc), 
  col=colorRampPalette(c("#191d73", "white", "#ed7901"))(n = 1000),
  labCol = colnames(pt.mat.TFs),
  labRow = lab.row,
  cexRow=.8# font size
) 
graphics.off()

# * Plot dynamics of temporal TFs in neural progenitors ----

# we use the list of all genes that came out of the differential gene expression analysis
plot.dynamics(genes = sort(age.TFs), title = "Figure_S5A.pdf", celltype = "Progenitor")

sink("sessionInfo.txt")
print(sessionInfo())
sink()

# Analysis of temporal TF expression in the fore-, mid- and hindbrain ----

# Analysis of scRNAseq data from developing fore-, mid- and hindbrain from the Linnarson lab (LaManno et al. 2020, bioRxiv)

# * Packages ----

# Use devtools to install hdf5r and loomR from GitHub

if("loomR" %in% rownames(installed.packages())) {
  print("loomR has already been installed!") 
} else {
  print("loomR needs to be installed! This requires devtools! Please install devtools if not present!")
  devtools::install_github(repo = "hhoeflin/hdf5r")
  devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
}

packages <- c("dplyr", "tidyr", "ggplot2")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(paste0(dir, '/output/'))

library(loomR)
library(dplyr)
library(ggplot2)

# * Functions ----

# Parameters are tissue ("Forebrain", "Midbrain", "Hindbrain" or "all"), genes (vector of genenames), 
# celltype ("Neuron" or "Radial glia"), latest.stage (integer up to which embryonic day values shall be plotted) 
# and specific (specific = TRUE means only tissues precisely called as the tissue argument will be considered, 
# specific = FALSE means all values that include the string passed to tissue will be considered 
# [e.g. tissue = Forebrain, specific = FALSE will consider "Forebrain" and "ForebrainVentral"] )

# generates circle plots in which expression level is color coded by size and percentage of expressing cells by color

plot.mean.exp <- function(tissue = "Midbrain", 
                          genes = genes, 
                          celltype = "Neuron", 
                          latest.stage = 15, 
                          specific = TRUE, 
                          normalize = TRUE){
  
  # generate expression matrix consisting of gene levels, timepoint, tissue and celltype 
  
  if(tissue == "all"){
    tissue = unique(sc.loom$col.attrs$Tissue[])
    name = "all"
  } else {
    name = tissue
  }
  
  if(specific == FALSE){
    tissue.id <- which(grepl(tissue, unique(sc.loom$col.attrs$Tissue[])) == TRUE)
    cell.id <- which(sc.meta$tissue %in% unique(sc.loom$col.attrs$Tissue[])[tissue.id])
  } else {
    if(length(tissue) > 1){
      cell.id <- which(sc.meta$tissue %in% tissue)
    } else {
      cell.id <- which(sc.meta$tissue == tissue)  
    }
  }
  
  gene.id <- which(sc.loom$row.attrs$Gene[] %in% genes)
  
  if(normalize == TRUE) {
    exp.mat <- sc.loom[["matrix"]][cell.id, gene.id] * sc.meta$normalization_factor[cell.id]
    name <- paste0(name, "_nomalized")
  } else {
    exp.mat <- sc.loom[["matrix"]][cell.id, gene.id]
  }
  
  colnames(exp.mat) <- sc.loom$row.attrs$Gene[gene.id]
  
  exp.mat <- cbind(exp.mat, sc.meta[cell.id, ])
  
  ###normalization function
  normalit<-function(x){
    x/(max(x))
  }
  
  # calculate number of cells expressing a gene for each timepoint
  num.cells.expressing <- do.call(cbind, lapply(genes, function(x) {
    gene.id <- which(sc.loom$row.attrs$Gene[] == x)
    exprs.values <- cbind(sc.loom[["matrix"]][cell.id, gene.id], sc.meta[cell.id, ]) %>%
      dplyr::filter(class == celltype) %>%
      dplyr::group_by(age)  
    
    colnames(exprs.values)[1] <- "expression"
    
    exprs.values <- exprs.values %>%
      dplyr::group_by(age) %>%
      dplyr::summarise(count2 = length(expression[expression > 0]))
    
    return(exprs.values$count2)
  }))
  
  colnames(num.cells.expressing) <- genes
  
  # calculate percentage of expressing cells for color map
  perc.expressing <- plyr::count(exp.mat, c("age", "class")) %>%
    dplyr::filter(class == celltype) %>% 
    dplyr::bind_cols(. ,as_tibble(num.cells.expressing)) %>%
    dplyr::mutate_each(function(v){v / .[,"freq"]}, all_of(genes)) %>%
    tidyr::gather(genename, percent, 4:(length(genes)+3), factor_key = TRUE) %>%
    dplyr::mutate(age = gsub("e", "", age)) %>%
    dplyr::mutate(age = as.numeric(age)) %>%
    dplyr::filter(age < latest.stage) %>%
    dplyr::pull(percent)
  
  # calculate average level of expression for the genes at each timepoint
  mat <- exp.mat %>%
    dplyr::filter(class == celltype) %>%
    dplyr::group_by(age) %>%
    dplyr::select(c(age, all_of(genes))) %>%
    dplyr::summarise_each(funs(mean)) %>%
    dplyr::mutate(age = gsub("e", "", age)) %>%
    dplyr::mutate(age = as.numeric(age)) %>%
    dplyr::filter(age < latest.stage) %>%
    dplyr::mutate_each(funs(normalit), all_of(genes)) %>%
    tidyr::gather(genename, measurement, 2:(length(genes)+1), factor_key = TRUE) %>%
    dplyr::mutate(percent = perc.expressing) ### merge with percent expressing cells
  
  # plot using ggplot
  
  gg <- ggplot(data = as.data.frame(mat), aes(x = factor(age), y = genename)) +
    geom_count(aes(size = measurement, color = percent)) +
    scale_size_area(max_size = 10) +
    scale_color_gradientn(colours = c("blue", "yellow"), limits = c(0,1)) +
    theme_classic() +
    ylab("genes") +
    xlab("embryonic days") +
    ggtitle(name) +
    theme(axis.text.y = element_text(size = 24),
          axis.title.y = element_text(size = 28),
          axis.text.x = element_text(size = 24, angle = 45, hjust = 1, vjust = 1),
          axis.title.x = element_text(size = 28),
          legend.title = element_text(size = 20),
          legend.title.align=0.5,
          legend.text = element_text(size = 20),
          strip.text = element_text(size = 24, face = "bold"),
          plot.title = element_text(face = "bold", size = 48))
  
  return(gg)
} 

# generates log-scaled z=scored heatmaps 

plot.heatmap.tTFs = function(tissue = "Midbrain", 
                             genes = genes, 
                             celltype = "Radial glia", 
                             latest.stage = 15, 
                             specific = TRUE, 
                             normalize = TRUE) {
  
  if(tissue == "all"){
    tissue = unique(sc.loom$col.attrs$Tissue[])
    name = "all"
  } else {
    name = tissue
  }
  
  if(specific == FALSE){
    tissue.id <- which(grepl(tissue, unique(sc.loom$col.attrs$Tissue[])) == TRUE)
    cell.id <- which(sc.meta$tissue %in% unique(sc.loom$col.attrs$Tissue[])[tissue.id])
  } else {
    if(length(tissue) > 1){
      cell.id <- which(sc.meta$tissue %in% tissue)
    } else {
      cell.id <- which(sc.meta$tissue == tissue)  
    }
  }
  
  gene.id <- which(sc.loom$row.attrs$Gene[] %in% genes)
  
  if(normalize == TRUE) {
    exp.mat <- t(sc.loom[["matrix"]][cell.id, gene.id] * sc.meta$normalization_factor[cell.id])
    name <- paste0(name, "_nomalized")
  } else {
    exp.mat <- t(sc.loom[["matrix"]][cell.id, gene.id])
  }
  
  rownames(exp.mat) <- sc.loom$row.attrs$Gene[gene.id]
  colnames(exp.mat) <- cell.id
  
  exp.mat.zscored <- t(exp.mat + 0.000001 %>% log()) %>% 
    scale(. , center=T, scale=T) %>%
    t()
  
  exp.mat.TFs <- data.frame(exp.mat.zscored) %>%
    tibble::rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    tidyr::gather(key = "cellname", value = "readcount", -gene) %>%
    dplyr::mutate(cellname = as.numeric(gsub("X", "", cellname))) %>%
    dplyr::mutate("timepoint" = sc.meta$age[.$cellname]) %>%
    dplyr::mutate(timepoint = as.numeric(gsub("e", "", timepoint)))  %>%
    dplyr::filter(timepoint < latest.stage) %>%
    dplyr::group_by(gene, timepoint) %>%
    dplyr::summarise(mean = mean(readcount)) %>%
    dplyr::ungroup() %>%
    tidyr::spread(gene, value = "mean") %>%
    tibble::column_to_rownames("timepoint") %>%
    t() %>% as.matrix()
  
  pdf(paste0("Figure_5B_",name,".pdf"), height=5, width=3)
  heatmap.plus::heatmap.plus(
    exp.mat.TFs[order.TFs.clusters, ],
    scale='row',
    Colv = NA,
    Rowv = NA, 
    col=colorRampPalette(c("#191d73", "white", "#ed7901"))(n = 1000),
    labCol = colnames(exp.mat.TFs),
    cexRow=.8# font size
  ) 
  graphics.off()
}

# * Load data ----
sc.loom <- connect(filename = paste0(dir, "/input/dev_all.loom"), mode = 'r+', skip.validate = TRUE)

# Generate sc.meta file by extracting parameters from connected sc.loom file
sc.meta <- data.frame(sc.loom$col.attrs$Age[], 
                      sc.loom$col.attrs$PseudoAge[], 
                      sc.loom$col.attrs$Tissue[], 
                      sc.loom$col.attrs$PseudoTissue[],
                      sc.loom$col.attrs$Class[],
                      sc.loom$col.attrs$Clusters[],
                      10000 / sc.loom$col.attrs$TotalUMI[])

colnames(sc.meta) <- c("age", "pseudoage", "tissue", "pseudotissue", "class", "clusters", "normalization_factor")

# * Plot expression ----

# Plot temporal TF dynamics in neurons 

genes <- c("Onecut1", "Onecut2", "Onecut3", "Pou2f2", "Zfhx3", "Zfhx4", 
           "Nfia", "Nfib", "Neurod2", "Neurod6")


pdf("Figure_2A_midbrain.pdf", useDingbats = FALSE)
print(plot.mean.exp(tissue = "Midbrain", genes = genes, celltype = "Neuron", latest.stage = 15, normalize = TRUE))
graphics.off()

pdf("Figure_2A_Hindbrain.pdf", useDingbats = FALSE)
print(plot.mean.exp(tissue = "Hindbrain", genes = genes, celltype = "Neuron", latest.stage = 15, normalize = TRUE))
graphics.off()

pdf("Figure_2A_Forebrain.pdf", useDingbats = FALSE)
print(plot.mean.exp(tissue = "Forebrain", genes = genes, celltype = "Neuron", latest.stage = 15, specific = FALSE, normalize = TRUE))
graphics.off()

pdf("Figure_2A_all.pdf", useDingbats = FALSE)
print(plot.mean.exp(tissue = "all", genes = genes, celltype = "Neuron", normalize = TRUE))
graphics.off()

# * Plot heatmaps for spinal cord TFs in other regions of the nervous system ----

plot.heatmap.tTFs(tissue = "Midbrain", genes = age.TFs, celltype = "Radial glia", normalize = TRUE)
plot.heatmap.tTFs(tissue = "Hindbrain", genes = age.TFs, celltype = "Radial glia", normalize = TRUE)
plot.heatmap.tTFs(tissue = "Forebrain", genes = age.TFs, celltype = "Radial glia", specific = FALSE, normalize = TRUE)

sink("sessionInfo.txt")
print(sessionInfo())
sink()

# Comparison with in-vitro RNAseq data ----

# Comparison of scRNAseq data with RNAseq data of ventral spinal cord differentiations from Rayon et al 2020
# data was downloaded from GEO database accession number GSE140748

# * Load data ----

dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(paste0(dir, '/output/'))

mouse <- read.table(paste0(dir, '/input/GSE140748_expression_matrix.mouse.abundance.tsv'), sep = '\t', row.names = 1, header = TRUE)

age.TFs.symbols.mouse <- fData(eset)$ensembl_gene_id[which(fData(eset)$current_gene_names %in% age.markers$Var1)]

# * Calculate Heatmaps ----

mouse.df <- mouse[age.TFs.symbols.mouse, ] %>% t() %>% data.frame() %>%
  tibble::rownames_to_column(var = "sample_name") %>%
  tidyr::separate(sample_name, sep = "_", c("species", "day", "repeat_number")) %>%
  dplyr::filter(day %in% c(0,1,2,3,4,5,6,7)) %>%
  dplyr::select(-c(species, repeat_number)) %>%
  dplyr::group_by(day) %>%
  dplyr::summarise_each(funs(mean)) %>%
  dplyr::mutate(day = as.character(day)) %>%
  mutate_if(is.numeric, function(x) x+1) %>%
  mutate_if(is.numeric, function(x) log(x)) %>%
  tibble::column_to_rownames(var = "day") %>%
  scale(. , center=T, scale=T) %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "day") %>%
  tidyr::gather(key = "genename", value = "readcounts", -day) %>%
  dplyr::mutate(genename = fData(eset)[genename, "current_gene_names"]) %>%
  tidyr::spread(genename, readcounts) %>%
  tibble::column_to_rownames(var = "day") %>%
  as.matrix() %>% t()

mouse.df <- mouse.df[, as.character(sort(as.numeric(colnames(mouse.df))))]

order.mouse.genes <- intersect(rownames(pt.mat.TFs)[dendextend::order.hclust(hc)], rownames(mouse.df))

# * Plot Data ----

pdf("Figure_5C_all_genes.pdf", height=5, width=3)
heatmap.plus::heatmap.plus(
  mouse.df[order.mouse.genes, ],
  scale='row',
  Colv = NA,
  Rowv = NA, 
  col=colorRampPalette(c("#191d73", "white", "#ed7901"))(n = 1000),
  cexRow=.8# font size
) 
graphics.off()

pdf("Figure_5C_TFs.pdf", height=5, width=3)
heatmap.plus::heatmap.plus(
  mouse.df[order.TFs.clusters, ],
  scale='row',
  Colv = NA,
  Rowv = NA, 
  col=colorRampPalette(c("#191d73", "white", "#ed7901"))(n = 1000),
  cexRow=.8# font size
) 
graphics.off()


# Analysis of temporal TFs in the mouse retina ----

## Analysis of scRNAseq data from the mouse retina (Clark et al. 2019, Neuron)
## dowloaded from https://github.com/gofflab/developing_mouse_retina_scRNASeq

rm(list=ls()) ## clears environment 

# * Packages ----

packages <- c("monocle3", "dplyr", "Seurat", "grid", "Matrix")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

library(Seurat)
library(monocle3)
library(dplyr)
library(grid)

# * Functions ----

# function for MultiBarHeatmap to sort cells in neuronal subtypes by time 
# this function was written by user arkal and shared via the Seurat Github forum (https://github.com/satijalab/seurat/issues/2201)

suppressPackageStartupMessages({
  library(rlang)
})

DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               additional.group.sort.by = NULL, 
                               cols.use = NULL,
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = TRUE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE, 
                               lines.width = NULL, 
                               group.bar.height = 0.02, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  
  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ", 
                paste(bad.sorts, collapse = ", "))
      }
    }
  }
  
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]  
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]  
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }
    
    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }
      
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        
        # Default
        cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))  
        
        #Overwrite if better value is provided
        if (!is_null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }
        
        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
        
        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}

# * Load 10X retina data and convert to Seurat object ----
dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(paste0(dir, '/input/'))

expression.mtx <- Matrix::readMM(file = "10x_mouse_retina_development.mtx")
pd <- read.table("10x_mouse_retina_development_phenotype.csv", sep = ",", header = TRUE, row.names = 1)
fd <- read.table("10x_mouse_retina_development_feature.csv", sep = ",", header = TRUE, row.names = 1)
rownames(fd) <- make.unique(as.character(fd$gene_short_name))

setwd(paste0(dir, '/output/'))

rownames(expression.mtx) <- rownames(fd)
colnames(expression.mtx) <- rownames(pd)

#Construct monocle cds and convert to Seurat object
monocle_cds <- new_cell_data_set(expression.mtx,
                                 cell_metadata = pd,
                                 gene_metadata = fd)

seurat.retina <- as.Seurat(monocle_cds, counts = "counts", data = NULL, assay = "RNA")

# * Filter cell types and stages ----
Idents(seurat.retina) <- "age"
seurat.df <- subset(seurat.retina, idents = c("E14", "E16", "E18", "P0"))

Idents(seurat.df) <- "CellType"
seurat.df <- subset(seurat.df, idents = c("RPCs", "Neurogenic Cells", "Photoreceptor Precursors", "Cones", "Rods", "Retinal Ganglion Cells", "Amacrine Cells", "Horizontal Cells"))

seurat.df[["percent.mt"]] <- PercentageFeatureSet(seurat.df, pattern = "^mt-")

dim(seurat.df)

seurat.df <- seurat.df %>% 
  subset(subset = nFeature_RNA > 800 & nFeature_RNA < 6000 & percent.mt < 6) %>%
  NormalizeData(verbose=FALSE) %>%
  FindVariableFeatures(selection.method = "vst", verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30)

# * Plot temporal TFs on UMAPs ----

neuronal.tTFs <- c("Onecut2", "Pou2f2", "Zfhx3", "Nfib")

png("Figure_S3C.png", width = 1024, height = 1024)
print(FeaturePlot(seurat.df, features = neuronal.tTFs))
dev.off()

png("Figure_S3B.png", width = 512, height = 512)
print(DimPlot(seurat.df, reduction = "umap", group.by = "CellType", label = FALSE, repel = TRUE) + NoLegend()) 
dev.off()

png("Figure_S3A.png", width = 512, height = 512)
print(DimPlot(seurat.df, reduction = "umap", group.by = "age", label = FALSE, repel = TRUE) + NoLegend())
dev.off()

# * Plot heatmap of temporal TF expression in neuronal subtypes ----
Idents(seurat.df) <- "CellType"
seurat.neurons <- subset(seurat.df, idents = c("Horizontal Cells", "Amacrine Cells", "Retinal Ganglion Cells", "Rods", "Cones")) 

seurat.neurons$CellType <- factor(x = seurat.neurons$CellType, levels = c("Horizontal Cells", "Amacrine Cells", "Retinal Ganglion Cells", "Cones", "Rods"))
seurat.neurons$age <- factor(x = seurat.neurons$age, levels = c("E14", "E16", "E18", "P0"))

genes <- c("Onecut2", "Pou2f2", "Zfhx3", "Nfib",  "Lhx1", "Pax6", "Pou4f2", "Thrb", "Nrl")

cluster.heatmap.bar <- DoMultiBarHeatmap(seurat.neurons, 
                                         features = genes, 
                                         group.by='CellType', 
                                         additional.group.by = 'age') 

pdf(paste0("Figure_S3D.pdf"), height = 5, width = 8, useDingbats=FALSE)
print(cluster.heatmap.bar)
dev.off()

sink("sessionInfo.txt")
print(sessionInfo())
sink()