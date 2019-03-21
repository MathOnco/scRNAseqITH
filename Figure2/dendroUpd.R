
install.packages("tripack")
install.packages("colorRamps")
install.packages("qgraph")
install.packages("plotrix")


library("colorRamps")
library(tripack)
library(RColorBrewer)



library(qgraph)
library("plotrix")


mat <- read.csv("DistanceMatrix_MeanUMI.csv",header=T)
graph <- hclust(as.dist(mat),  method = "cen")
plot(graph,hang = -1, cex = 0.8, ylab = "distance")

dist_mi <- 1/as.matrix(mat)

size <- as.matrix(read.csv("clustersize.csv",header=F))

colscale <- as.matrix(read.csv("colorListHEXupd.csv",header=F))

pdf("dendrogram_out_upd.pdf",width=16,height=12)
qgraph(dist_mi, layout='spring', vsize=size , color=colscale , labels=TRUE, label.prop = 1.5, label.color="white") 
dev.off()

