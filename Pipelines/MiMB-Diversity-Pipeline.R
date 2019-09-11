library(future)
library(Matrix)
library(dplyr)
library(readr)
library(rdetools)
library(data.table)
library(ggplot2)
library(iterators)
library(bigmemory)
library(Seurat)
sessionInfo()

# Set the working directory to where the dataset was downloaded and un-zipped
setwd(" ")
plan("multiprocess", workers=2)

##### Step 1: Load all the datasets
# Navigate to the directory containing the three files that describe the gene/cell matrix (barcodes.tsv, genes.tsv, matrix.mtx)
bmmc027pre.data <- Read10X(data.dir = "./AML027pre/outs/filtered_gene_bc_matrices/GRCh38/")
bmmc027pre <- CreateSeuratObject(counts=bmmc027pre.data, project="AML027pre")
bmmc027pre <- RenameCells(object=bmmc027pre, add.cell.id="AML027pre")

bmmc027post.data <- Read10X(data.dir="./AML027post/outs/filtered_gene_bc_matrices/GRCh38/")
bmmc027post <- CreateSeuratObject(counts=bmmc027post.data, project="AML027post")
bmmc027post <- RenameCells(object=bmmc027post, add.cell.id="AML027post")

bmmc035pre.data <- Read10X(data.dir="./AML035pre/outs/filtered_gene_bc_matrices/GRCh38/")
bmmc035pre <- CreateSeuratObject(counts=bmmc035pre.data, project = "AML035pre")
bmmc035pre <- RenameCells(object=bmmc035pre, add.cell.id="AML035pre")

bmmc035post.data <- Read10X(data.dir="./AML035post/outs/filtered_gene_bc_matrices/GRCh38/")
bmmc035post <- CreateSeuratObject(counts=bmmc035post.data, project="AML035post")
bmmc035post <- RenameCells(object=bmmc035post, add.cell.id="AML035post")

bmmcHealthy1.data <- Read10X(data.dir="./healthy1/outs/filtered_gene_bc_matrices/GRCh38/")
bmmcHealthy1 <- CreateSeuratObject(counts=bmmcHealthy1.data, project="Healthy1")
bmmcHealthy1 <- RenameCells(object=bmmcHealthy1, add.cell.id="Healthy1")

bmmcHealthy2.data <- Read10X(data.dir="./healthy2/outs/filtered_gene_bc_matrices/GRCh38/")
bmmcHealthy2 <- CreateSeuratObject(counts=bmmcHealthy2.data, project="Healthy2")
bmmcHealthy2 <- RenameCells(object=bmmcHealthy2, add.cell.id="Healthy2")


# see number of cells per dataset
# dim(x.data) gives [1] number of genes and  [2] number of single cells
dim(bmmc027pre.data)
dim(bmmc027post.data)
dim(bmmc035pre.data)
dim(bmmc035post.data)
dim(bmmcHealthy1.data)
dim(bmmcHealthy2.data)

##### Step 2: Aggregate Datafiles
# aggregate iteratively using merge()
bmmcAll.combined1 <- merge(x=bmmc027pre, y=bmmc035pre, add.cell.id1 = "027pre", add.cell.id2 = "035pre", project = "allAML")
bmmcAll.combined2 <- merge(x=bmmcAll.combined1, y=bmmcHealthy1, add.cell.id2 = "H1", project = "allAML")
bmmcAll.combined3 <- merge(x=bmmcAll.combined2, y=bmmcHealthy2, add.cell.id2 = "H2", project = "allAML")
bmmcAll.combined4 <- merge(x=bmmcAll.combined3, y=bmmc027post, add.cell.id2 = "027post", project = "allAML")
bmmcAll.combined <- merge(x=bmmcAll.combined4, y=bmmc035post, add.cell.id2 = "035post", project = "allAML")

##### Step 3: Quality Control Step
# identify mitochondrial genes
mito.genes <- grep(pattern = "^MT-", bmmcAll.combined@assays$RNA@counts@Dimnames[[1]], value = TRUE)
percent.mito <- Matrix::colSums(x=GetAssayData(object=bmmcAll.combined, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=bmmcAll.combined, slot='counts'))
bmmcAll.combined[['percent.mito']] <- percent.mito

# plotting the nFeatures, nCount, & % mito across all samples
VlnPlot(object=bmmcAll.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

# calculating cut-offs for nFeature content (mean +/- 2*std)
nFeatUpper <- mean(bmmcAll.combined@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(bmmcAll.combined@meta.data$nFeature_RNA, na.rm=TRUE)
nFeatLower <- mean(bmmcAll.combined@meta.data$nFeature_RNA, na.rm=TRUE) - 2*sd(bmmcAll.combined@meta.data$nFeature_RNA, na.rm=TRUE) 
perMitoUpper <- 0.05

# plotting data with cutoffs
VlnPlot(object=bmmcAll.combined, features = c("nFeature_RNA"))+geom_hline(yintercept=nFeatUpper, linetype="dashed", color="red", size=1)

VlnPlot(object=bmmcAll.combined, features = c("percent.mito"))+geom_hline(yintercept=perMitoUpper, linetype="dashed", color="red", size=1)

# subset the data to use only cells within the cutoffs for nFeature & percent.mito
bmmcAll.combined <- subset(x=bmmcAll.combined, subset=nFeature_RNA > nFeatLower & nFeature_RNA < nFeatUpper & percent.mito < 0.05)
# can check dim(bmmcAll.combined) before and after subsetting dataset, to confirm cells were removed (13878 --> 12662)

##### Step 4: Normalize the data
bmmcAll.combined <- NormalizeData(object=bmmcAll.combined, normalization.method = "LogNormalize", scale.factor = 10000)

##### Step 5: Determine variable genes across cells (Optional)
bmmcAll.combined <- FindVariableFeatures(object = bmmcAll.combined, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf), do.plot=TRUE)
length(x = VariableFeatures(object=bmmcAll.combined)) 

##### Step 6: Scale data and remove unwanted variation
bmmcAll.combined <- ScaleData(object=bmmcAll.combined, vars.to.regress = c("nCount_RNA", "percent.mito"))


##### Step 7: Perform linear dimensional reduction
bmmcAll.combined <- RunPCA(object=bmmcAll.combined, features=VariableFeatures(object=bmmcAll.combined), verbose=TRUE, npcs=200)
bmmcAll.combined <- ProjectDim(object=bmmcAll.combined, verbose=FALSE)
#<--- here 
##### Step 8: Cluster Data - Graph-Based Cluster (like Loupe Browser)
# following along: https://satijalab.org/seurat/pbmc3k_tutorial.html 
bmmcAll.combined <- FindNeighbors(object=bmmcAll.combined, dims = 1:50)
bmmcAll.combined <- FindClusters(object=bmmcAll.combined, resolution = 0.6)

# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
#
# Number of nodes: 12662
# Number of edges: 834022
#
# Running Louvain algorithm...
# 0%   10   20   30   40   50   60   70   80   90   100%
# [----|----|----|----|----|----|----|----|----|----|
# **************************************************|
# Maximum modularity in 10 random starts: 0.9153
# Number of communities: 18
# Elapsed time: 2 seconds

# counts the total number of clusters identified
nlevels(bmmcAll.combined@meta.data$RNA_snn_res.0.6) # 18

# visualize clustering with UMAP
bmmcAll.combined <- RunUMAP(object=bmmcAll.combined, reduction = "pca", dims = 1:50)
pdf("./BMMC-Pre+Post+Healthy_UMAP_ByClusters.pdf")
dPlot <- DimPlot(bmmcAll.combined, reduction.use="umap")
print(dPlot)
dev.off()
pdf("./BMMC-Pre+Post+Healthy_UMAP_ByIdentity.pdf")
dPlot <- DimPlot(bmmcAll.combined, reduction.use="umap", group.by='orig.ident')
print(dPlot)
dev.off()


###### Step 9: Diversity Scoring
# collect number of cells per cluster per type
file1 <- "./BMMCoutputDataNames.csv";
file2 <- "./BMMCoutputData.csv"
cat(bmmcAll.combined@meta.data$RNA_snn_res.0.6, file=file1, sep=",\n")
cat(bmmcAll.combined@meta.data$orig.ident, file=file2, sep=",\n")
mydat1 <- read.csv(file2)
mydat2 <- read.csv(file1)
fulldat <- cbind(mydat2[1],mydat1[1])
# count number of unique barcodes per cluster
fulltab <- as.data.table(fulldat)
names(fulltab)[1] <- paste("cluster")
names(fulltab)[2] <- paste("type")
# group cells by cluster
group_by(fulltab, cluster)
tabPerClus <- fulltab %>% group_by(cluster) %>% count()
# then group by type
type <- sub("\\_.*","",fulltab$UMI)
fulltab <- cbind(fulltab, type)
tabPerClusType <- fulltab %>% group_by(cluster, type) %>% count()
# covert the cell counts to frequencies used to calculate the diveristy score
tabPerType <- fulltab %>% group_by(type) %>% count()
tmp <- matrix(ncol=1, nrow=dim(tabPerClusType)[1])
for(i in 1:dim(tabPerClusType)[1]){
  if(tabPerClusType$type[i]==tabPerType$type[1]) {
    ## if AML027pre
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[1]
  } else if (tabPerClusType$type[i]==tabPerType$type[2]) {
    ## if AML035pre
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[2]
  } else if (tabPerClusType$type[i]==tabPerType$type[3]) {
    ## if Healthy 1
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[3]
  } else if (tabPerClusType$type[i]==tabPerType$type[4]) {
    ## if Healthy 2
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[4]
  } else if (tabPerClusType$type[i]==tabPerType$type[5]) {
    ## if AML027post
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[5]
  } else if (tabPerClusType$type[i]==tabPerType$type[6]) {
    ## if AML035post
    tmp[i] <- tabPerClusType$n[i]/tabPerType$n[6]
  }else {
    tmp[i] <- 0.0
  }
}
tmp <- as.data.frame(tmp)
names(tmp)[1] <- paste("freq")
tabPerClusType <- cbind(as.data.frame(tabPerClusType), tmp)

# filter/collect per cell type (names must match names you set up in Seurat object)
AML1 <- filter(tabPerClusType,type=="AML027pre")
AML2 <- filter(tabPerClusType,type=="AML035pre")
H1 <- filter(tabPerClusType,type=="Healthy1")
H2 <- filter(tabPerClusType,type=="Healthy2")
postAML1 <- filter(tabPerClusType,type=="AML027post")
postAML2 <- filter(tabPerClusType,type=="AML035post")

# Function to calculate the divesrity index
calcqD <- function(dataset, q){
  diversity <- 0.0;
  for(row in 1:dim(dataset)[1]){
    diversity <- diversity + (dataset$freq[row])^q
  }
  diversity <- diversity^(1/(1-q))
}

# Calculate the diversity index over a range of q (set by qRange)
qRange <- logspace(-2, 2, n = 1000)
qDAML1 <- calcqD(AML1,qRange)
qDAML1tab <- as.data.frame(cbind(qRange,qDAML1))
qDAML2 <- calcqD(AML2,qRange)
qDAML2tab <- as.data.frame(cbind(qRange,qDAML2))

qDH1 <- calcqD(H1,qRange)
qDH1tab <- as.data.frame(cbind(qRange,qDH1))
qDH2 <- calcqD(H2,qRange)
qDH2tab <- as.data.frame(cbind(qRange,qDH2))

qDpostAML1 <- calcqD(postAML1,qRange)
qDpostAML1tab <- as.data.frame(cbind(qRange,qDpostAML1))
qDpostAML2 <- calcqD(postAML2,qRange)
qDpostAML2tab <- as.data.frame(cbind(qRange,qDpostAML2))

### Plots the spectrum of diversity scores
# Plotting all data
pdf("./BMMC_Diversity-Curves.pdf")
qDPlot <- ggplot(data=qDAML1tab, aes(x=qRange, y=qDAML1),color="blue") +
  geom_smooth(size=4) + scale_x_log10() +
  geom_smooth(data=qDAML2tab,aes(x=qRange,y=qDAML2),size=4,color="navy") +
  geom_smooth(data=qDH1tab,aes(x=qRange,y=qDH1),size=4,color="dark gray") + 
  geom_smooth(data=qDH2tab,aes(x=qRange,y=qDH2),size=4,color="black") +
  geom_smooth(data=qDpostAML2tab,aes(x=qRange,y=qDpostAML1),size=4,color="purple") + 
  geom_smooth(data=qDpostAML2tab,aes(x=qRange,y=qDpostAML2),size=4,color="blueviolet") + 
  labs(title="Diversity Score",x="q", y = "qD")
print(qDPlot)
dev.off()
