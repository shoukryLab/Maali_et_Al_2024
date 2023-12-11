library(dplyr)
library(magrittr)
library(limma)
library(data.table)
library(factoextra)
library(tidyverse)
library(pheatmap)
library(DESeq2)
library(fgsea)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(igraph)
library(DOSE)
library(VennDiagram)
library(ggbiplot)

###Function dgeWay::
dgeWay <- function(count_mat){
  library(edgeR)
  
  #This is the limma pipeline!@!@!
  dgee <- DGEList(count_mat)
  keep <- rowSums(cpm(dgee) > 1) >= 3 #Discuss w Youssef @!@
  #keep <- rowSums(count_mat) > 50
  dgee <- dgee[keep, ,keep.lib.sizes = FALSE]
  
  dgee$samples$lib.size <- colSums(dgee$counts)
  dgee <- calcNormFactors(dgee)
  
  return(dgee)
}
#####

#load data:
setwd("C:/Users/Omar/Desktop/counts") #Change to directory w files

data <- fread("countsAll.txt", sep = " ",
              col.names = T)

#remove genes with expression == 0 accross all samples::
which(rowSums(data[,-1]) == 0) -> zeroExpress

data <- data[-zeroExpress,]

#now begin with normalization and everything else::
colnames(data) <- c("Ensembl", "24-9",
                    "L6", "L9","L10","S1B","S1C","72-6","72-7",
                    "24-6","24-7","24-8", "S3D")

rownames(data) <- data$Ensembl
data <- data[,-1]

#remove also last 4 rows: __ambiguous, etc...
data <- data[-c(nrow(data):(nrow(data)-3)),]

#save::
# fwrite(data, file = "counts_notNormalized.txt", sep = " ",
#        row.names = T, col.names = T)

#factors for samples::
tissue <- c(rep("L",4), rep("S", 2), rep("L", 5), "S")
time <- c(24, rep(72, 3), rep(0, 2), rep(72, 2), rep(24,3), 0)

#MOUSEID
#calculate normalisation factors for expression matrix:
data_dge <- dgeWay(count_mat = data)

#put tissue and time factors into DGE object::
#Put it into the dge object for design matrix and downstream analysis
data_dge$samples$tissue <- as.factor(tissue)
data_dge$samples$time <- as.factor(time)
data_dge$samples$tissue_time <- paste(tissue, time, sep = "_")
designo <- model.matrix(~ 0 + tissue_time, data = data_dge$samples)

#Create voom object to normalize and visualize PCA::
voomo <- voom(data_dge, design = designo, plot = T)

# fit <- lmFit(voomo, designo)
# fit <- eBayes(fit)
# topTable(fit, adjust.method = "BY", p.value = 0.05, 
#          number = nrow(data)) -> topi_005
# 
# topi_005$gene <- topi_005 %>% rownames()
#get names of genes::
geneNames <- fread("counts_merged_NS_l6.bam_sam.txt") %>% 
  dplyr::select(V1, V2)

#Visualize with PCA::
pcaObj <- prcomp(t(voomo$E), scale. = TRUE, center = TRUE)

#try using ggplot::
p.pca <- ggbiplot(pcaObj,
                  groups=as.factor(time),
                  choices=c(1,2),
                  circle=FALSE,
                  ellipse=TRUE,
                  var.axes=FALSE,
                  point.size = 0.5
)

#use colors as defined by rgb intensity::
col_24 <- rgb(red = 255, green = 129, blue = 129, maxColorValue = 255)
col_72 <- rgb(red = 144, green = 191, blue = 249, maxColorValue = 255)
col_0 <- rgb(red = 147, green = 247, blue = 164, maxColorValue = 255)

p.pca + 
  geom_point(aes(color = factor(as.character(time))), 
             size = 2) +
  scale_color_manual(name="Time Points", 
                     values=c(col_0,  col_24, col_72)) +
  ggtitle(label="PCA") +theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank()
  ) + theme(
    axis.title.x = element_text(size = 7),
    axis.title.y = element_text(size = 7)
  ) -> pca_f

