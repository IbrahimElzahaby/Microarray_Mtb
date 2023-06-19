# Microarray gene expression analysis of Mycobacterium tuberculosis (Mtb)



## Project description:

Our project investigated 177 genes across 287 conditions using different treatments, 
antibiotics, and chemical agents to identify the genes involved in a biological process that 
could be inferred the similar expression levels and the regulatory networks between the genes 
among the conditions.



## Analysis workflow:

#### 1. Initial analysis:

We imported data into R for analysis. We determined the number for genes and conditions in our 
data. We made a Q-Qplot to test for normality in our data. We identified and imputed missing 
data points using impute R package.

#### 2. PCA and Gene Correlation analysis:

We analysed variance-covariance structure and determined percent of variance explained by each 
Principal component. We also analysed the degree of association between genes using a heatmap 
of gene correlation.

#### 3. Gene distance Clustering Analysis:

We performed hierarchical Euclidean distance clustering of genes and made a dendrogram to 
visualize the grouping of genes in our data.

#### 4. Network Inference Analysis:

We used minet package for constructing networks. We used three methods for comparison; clr, 
mrnet, and aracne. We wrote out the resulting networks to a table for further analysis in 
CytoScape.

#### 5. Network Property Analysis with CytoScape:

We imported the network table into CytoScape, we optimised the threshold for each network 
until we obtained an easily interpretable network. We settled on a threshold of 0.25 for clr, 
0.1 for mrnet, for some reason changing threshold did not seem to affect aracne network. For 
this reason, we did not proceed with network from aracne.

#### 6. Network Cluster Analysis:

We explored network properties using analyse network tool in CytoScape. We used clusterONE app 
to determine clusters in our network.

#### 7. GO Biological Processes:

We analysed the GO biological processes and molecular functions associated with the clusters 
using Bingo.
