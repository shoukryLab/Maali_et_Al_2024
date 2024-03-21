#This script will show how we analysed our data along with publicly available data from GSE180824
library(dplyr)
library(magrittr)
library(limma)
library(data.table)
#library(factoextra)
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

#load files and cbind::
data <- fread("countsAll.txt")

#integrate the public GSE data::
GSE <- fread("GSE180824_Mackey_Neutrophil_maturity_raw_counts.csv") 

#Select only samples of interest - Liver Neutrophils + Spleen for balanced design::
cols_GSE <- colnames(GSE) %>% 
  grep(pattern = "PBS_SP_Hi|PBS_LV|PBS_BL_Neut", 
       ignore.case = T, value = T)

#We have to include some other technical replicates (from the batch named GSE), 
#so that we can properly include technical replicates in our linear model created
#below. Otherwise, we will have a design matrix that is not of full rank.
GSE %<>%
  dplyr::select(Geneid, cols_GSE) 

#remove genes with expression == 0 across all samples in OUR data::
which(rowSums(data[,-1]) == 0) -> zeroExpress
data <- data[-zeroExpress,]

#Do same w/ GSE data::
which(rowSums(GSE[,-1]) == 0) -> zeroExpress
GSE <- GSE[-zeroExpress,]

###### VOOM VOOM VOOM
#now begin with normalization and everything else::
colnames(data) <- c("Gene", "24-9",
                    "L6", "L9","L10","S1B","S1C","72-6","72-7",
                    "24-6","24-7","24-8", "S3D")

#remove last 4 rows from data::
data <- data[-c(nrow(data):(nrow(data)-3)),]

#convert Ensembl -> gene IDs using this file::
forConversion <- fread("forConversion.txt")

#identify dupGenes for exclusion from merge. Will convert to Gene names 
#here for visualization of Volcano plots::
data %>%
  mutate(Gene = sub("\\.\\d+$", "", Gene)) %>%
  inner_join(GSE, by = c("Gene" = "Geneid")) %>%
  inner_join(forConversion, by = c("Gene"="V1")) %>%
  .$V2 %>% table %>% as.data.frame %>% filter(Freq > 1) %>% 
  .$. %>% as.character() -> dupGenes

data %<>%
  mutate(Gene = sub("\\.\\d+$", "", Gene)) %>%
  inner_join(GSE, by = c("Gene" = "Geneid")) %>%
  inner_join(forConversion, by = c("Gene"="V1")) %>%
  filter(!V2 %in% dupGenes) %>% 
  column_to_rownames("V2") %>%
  dplyr::select(-c("Gene")) 

#MOUSEID
#calculate normalisation factors for expression matrix:
data_dge <- dgeWay(count_mat = data)

#create labels for samples, concatenating GSE data in there as well::
tissue <- c(rep("L",4), rep("S", 2), rep("L", 5), "S") 

#Add GSE samples to tisse::
sapply(colnames(data) %>% grep(pattern= "PBS", value = T) %>% 
         strsplit("_"), function(i){
           i %>% extract2(2)
         }) -> tiss_GSE

#replace "LV" with "L" as with our other samples::
tiss_GSE %<>% gsub("LV","L",.)
tiss_GSE %<>% gsub("SP", "S", .) #Same w SP

tissue <- c(tissue, tiss_GSE)

#Now for time variable
time <- c(24, rep(72, 3), rep(0, 2), rep(72, 2), rep(24,3), 0) 

#Add GSE samples to time variable as well::
time_GSE <- rep(0, length(tiss_GSE))

time <- c(time, time_GSE)

#create tissue_time var::
tissue_time <- paste(tissue, time, sep = "_") %>% as.factor()

#Have to model with batch given that public data is being used
batcho <- c(rep("Us",4), rep("Us", 2), rep("Us", 5), "Us",
            rep("GSE",length(tiss_GSE)))

#put tissue and time factors into DGE object::
#Put it into the dge object for design matrix and downstream analysis
data_dge$samples$tissue <- as.factor(tissue)
data_dge$samples$time <- as.factor(time)
data_dge$samples$tissue_time <- paste(tissue, time, sep = "_")

data_dge$samples$batcho <- as.factor(batcho)

#Batch is included in the linear model::
#removeBatchEffect: "This function is not intended to be used prior 
#to linear modelling. For linear modelling, it is better to include #
#the batch factors in the linear model."
designo <- model.matrix(~0 + tissue_time + batcho)

#Create voom object to normalize and visualize PCA::
#Have to restore matrix to full rank in order to correctly estimate
#coefficients.
voomo <- voom(data_dge, design = designo, plot = T)

#Visualize with PCA::
pcaObj <- prcomp(t(voomo$E), scale. = TRUE, center = TRUE)

#try using ggplot::
tissue_time_batch <- as.factor(c(paste0(tissue_time,"_",batcho)))
p.pca <- ggbiplot(pcaObj,
                  groups=tissue_time_batch,
                  choices=c(1,2),
                  circle=FALSE,
                  ellipse=TRUE,
                  var.axes=FALSE,
                  point.size = 0.5
)

#There seems to be a very strong batch effect, fix this::
dataBC <- sva::ComBat(voomo$E, 
                      batch = as.character(batcho))
# dataBC <- removeBatchEffect(voomo$E,
#                             batch = batcho) #Alternative method

#Now run the PCA on the batch corrected data::
pcaBC <- prcomp(t(dataBC), scale. = T, center = T)

p.pca <- ggbiplot(pcaBC,
                  groups=tissue_time_batch,
                  choices=c(1,2),
                  circle=FALSE,
                  ellipse=TRUE,
                  var.axes=FALSE,
                  point.size = 0.5
)

#complete VOOM/LIMMA to get contrasts + DEG's::
contrastMat <- makeContrasts("72v0" = tissue_timeL_72 - tissue_timeL_0,
                             "24v0" = tissue_timeL_24 - tissue_timeL_0,
                             levels = designo)

fito <- lmFit(voomo, design = designo, weights = voomo$weights)

fitContrast <- contrasts.fit(fito, contrasts = contrastMat)
fitContrast <- eBayes(fit = fitContrast)

fitContrast$tops <- sapply(colnames(fitContrast), function(cont){
  topTable(fitContrast, coef = cont, number = Inf,
           adjust.method = "BH", p.value = 0.05)
}, simplify = FALSE)

#####
#GSEA
#was done using fgsea library and a ranked list based on the 
#coefficients in the fitContrast objcet. 
#Gene set Using pre-converted Biological 
#processes gene sets. File Name: "m5.go.bp.v2022.1.Mm.symbols.gmt"
#####
