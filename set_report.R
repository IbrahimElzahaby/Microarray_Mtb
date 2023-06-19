---
Title: "MicroArray gene expression network analysis of Mycobacterium tuberculosis (Mtb)"
Authers: "Ibrahim ElZahaby & Luka Lmelias"
Team: "Group_5"
  ---
# Libraries required
library(impute)
library(dendextend)
library(minet)
library(RColorBrewer) 

# Read in the file
data <- read.csv('set_5.csv', row.names = 1)

# explore the data

# count the number of genes and contions
dim(data)
# get the condition names
colnames(data)

# Check the size of missing data
missing_data = which(is.na(data), arr.ind = TRUE) # get row and col with NA

# Impute missing data
data_impute = impute.knn(as.matrix(data))
data = data_impute$data
after_impute = which(is.na(data), arr.ind = TRUE)

# Check for data normality. 
boxplot(data) # turns out to be too large to interpret, we try QQplot

qqnorm(data)

## Explore the PCA spaces
pca_gene = prcomp(as.matrix(data), center = TRUE, scale. = TRUE)

pca_summary = summary(pca_gene) # we need the importance of each PC

plot(abs(pca_gene$rotation[,1]),abs(pca_gene$rotation[,2]))

barplot(abs(pca_gene$x)) # check the distribution of PCs


# check the % variance explained by each PC
var_expl = (pca_gene$sdev^2)/sum(pca_gene$sdev^2)*100
barplot(var_expl, main = 'Variance explained',
        xlab = 'componet',
        ylab = '% variance explained')


dim(pca_gene$x)

scores = pca_gene$x[,1:4] # we added 4 PCs cause each PC explains only little % variance

# The fisrt and second PCA
plot(scores[,1:2], main = 'PCA scores of the 1st and 2nd PCs',
     xlab = paste0("PC1 (",format(100* pca_summary$importance[2,1],digits=2),"%)"),
     ylab = paste0("PC2 (",format(100* pca_summary$importance[2,2],digits=2),"%)"),
     pch=21
)

# The third and fourth PCA
plot(scores[,3:4], main = 'PCA scores of The 3rd and 4th PCs',
     xlab = paste0("PC3 (",format(100* pca_summary$importance[2,3],digits=2),"%)"),
     ylab = paste0("PC4 (",format(100* pca_summary$importance[2,4],digits=2),"%)"),
     pch=21
)

# Heatmap of gene correlation
colors<-brewer.pal(7,"RdBu")

heatmap(cor(t(data)),symm=TRUE, col=colors, 
        main = "Heatmap of gene correlation")
legend(x = "topright", legend = c("negative", "medium", "high"), fill = colors,
       cex = 0.7) # we could not exactly get the legend colors to show blue when there is positive correlation

# Hierarchical  clustering
distance_gene = dist(data, method = 'euclidean')
hclust_gene = hclust(as.dist(distance_gene), method = 'average') # We could also compare 'average' vs 'complete' if time allows

# Plot the dendrogram
plot(hclust_gene, main = 'clustering in data space',
     xlab = '',  sub = 'average hclust distance') # data space

# Network Inference
data = t(as.matrix(data))
nbins <- sqrt(nrow(data)) # compute mutual information

# method 1 clr
clr.net <- minet(data, method="clr", 
                 estimator="mi.empirical",
                 disc="equalfreq", nbins=nbins)

clr.net[1:5,1:5]

# creating network table
clr.th <- clr.net

diag(clr.th) <- 0

clr.th[lower.tri(clr.th )] <- 0

clr.th[which(clr.th>0)] <- 1
xx <- which(clr.th==1, arr.ind=TRUE)
write.table(cbind(rownames(clr.th)[xx[,1]],
                  rownames(clr.th)[xx[,2]]),
            file="clr.csv", col.names=FALSE, 
            row.names=FALSE, sep=",", quote=FALSE)
length(xx)

# set fixed threshold
thrsh <- 0.25
clr.th <- clr.net
clr.th[which(clr.net<= thrsh)] <- 0
clr.th[which(clr.net> thrsh)] <- 1
diag(clr.th) <- 0
length(clr.th)

#to write to a file (that can be imported in cytoscape)
clr.th[lower.tri(clr.th )] <- 0
xx <- which(clr.th==1, arr.ind=TRUE)
write.table(cbind(rownames(clr.th)[xx[,1]],
                  rownames(clr.th)[xx[,2]]),
            file="clr_th.csv", col.names=FALSE, 
            row.names=FALSE, sep="\t", quote=FALSE)

length(xx)

# method 2 mrnet 
mrnet.net <- minet(data, method="mrnet", 
                   estimator="mi.empirical", 
                   disc="equalfreq", nbins=nbins)
mrnet.net[1:5,1:5]#have a look at the adjacency matrix

# create the table
mrnet.th <- mrnet.net
diag(mrnet.th) <- 0 
mrnet.th[lower.tri(mrnet.th )] <- 0
mrnet.th[which(mrnet.th>0)] <- 1
xx <- which(mrnet.th==1, arr.ind=TRUE)
write.table(cbind(rownames(mrnet.th)[xx[,1]],
                  rownames(mrnet.th)[xx[,2]]),
            file="mrnet.csv", col.names=FALSE, 
            row.names=FALSE, sep=",", quote=FALSE)
length(xx)

# fixing the the threshold
thrsh <- 0.1
mrnet.th <- mrnet.net
mrnet.th[which(mrnet.net<= thrsh)] <- 0
mrnet.th[which(mrnet.net> thrsh)] <- 1
diag(mrnet.th) <- 0
length(mrnet.th)

#to write to a file (that can be imported in cytoscape)
mrnet.th[lower.tri(mrnet.th )] <- 0
xx <- which(mrnet.th==1, arr.ind=TRUE)
write.table(cbind(rownames(mrnet.th)[xx[,1]],
                  rownames(mrnet.th)[xx[,2]]),
            file="mrnet_th.csv", col.names=FALSE, 
            row.names=FALSE, sep=",", quote=FALSE)
length(xx)

# method 3 Aracne
aracne.net <- minet(data, method="aracne",
                    estimator="mi.empirical",
                    disc="equalfreq", nbins=nbins)
aracne.net[1:5,1:5]

# write table without threshold
aracne.th <- aracne.net
diag(aracne.th) <- 0
aracne.th[lower.tri(aracne.th )] <- 0
aracne.th[which(aracne.th>0)] <- 1
xx <- which(aracne.th==1, arr.ind=TRUE)
write.table(cbind(rownames(aracne.th)[xx[,1]],
                  rownames(aracne.th)[xx[,2]]),
            file="aracne.csv", col.names=FALSE, 
            row.names=FALSE, sep=",", quote=FALSE)
length(xx)

# fix threshold
thrsh <- 0.1
aracne.th <- aracne.net
aracne.th[which(aracne.net<= thrsh)] <- 0
aracne.th[which(aracne.net> thrsh)] <- 1
diag(aracne.th) <- 0
length(aracne.th)

#to write to a file (that can be imported in cytoscape)
aracne.th[lower.tri(aracne.th )] <- 0
xx <- which(aracne.th==1, arr.ind=TRUE)
write.table(cbind(rownames(aracne.th)[xx[,1]],
                  rownames(aracne.th)[xx[,2]]),
            file="aracne_th.csv", col.names=FALSE, 
            row.names=FALSE, sep=",", quote=FALSE)
length(xx)


