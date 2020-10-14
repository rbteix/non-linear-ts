#0. Dimensionality Reduction methods.
#0.1 This code has the objective of applying dimensionality reduction methods for embedding of phase spaces.

setwd("C://Users//Guilherme//Documents//MSc. Unifesp - ITA//6. Séries Temporais Não-Lineares//Artigo//Código R")

library(dimRed)
library(Rdimtools)
#1. Use the series of failure from the transformed data set.

#1.1 Define function to perform embedding
embedding<- function(x,m,d){
  n<-length(x)
  ne<- n-(m-1)*d # number of vectors in the reconstructed phase space
  y <- matrix(0,ne,m) # ne rows, m column
  for(i in 1:m) {y[,i] <- x[((i-1)*d+1):(ne+(i-1)*d)]}
  y
}

trajectory_matrix <- as.data.frame(embedding(unlist(data_fail), 50, 1))

#1.2 Transforms trajectory matrix into a dimRed data object.

trajectory_matrix <- dimRedData(trajectory_matrix)

tentative_dimensions <- 1:50

#2. Creates code that optimized dimension

#a. kPCA
optimim_kPCA <- function(dim){ #Builds a wrapper for optimizing kPCA
  
  #Defines trajectore matrix inside function.
  
  trajectory_matrix <- as.data.frame(embedding(unlist(data_fail), 50, 1))
  
  #1.2 Transforms trajectory matrix into a dimRed data object.
  
  trajectory_matrix <- dimRedData(trajectory_matrix)
  
  #2. 
  real_dim <- round(dim,0)
  kpCA_result <- embed(trajectory_matrix, .method = 'kPCA', ndim = real_dim)
  
  AUC <- (-1)*quality(kpCA_result, .method = "mean_R_NX")
  
  #Returns the area under the curve.
  
  return(AUC)
  
}

kPCA_results <- lapply(tentative_dimensions, optimim_kPCA)
windows()
plot(tentative_dimensions, (-1)*unlist(kPCA_results))


#b. nMDS

optimim_nMDS <- function(dim){ #Builds a wrapper for optimizing kPCA
  
  #Defines trajectorz matrix inside function.
  
  trajectory_matrix <- as.data.frame(embedding(unlist(data_fail), 50, 1))
  
  #1.2 Transforms trajectory matrix into a dimRed data object.
  
  trajectory_matrix <- dimRedData(trajectory_matrix)
  
  #2. 
  real_dim <- round(dim,0)
  nMDS_result <- embed(trajectory_matrix, .method = 'nMDS', ndim = real_dim)
  
  AUC <- (-1)*quality(nMDS_result, .method = "mean_R_NX")
  
  #Returns the area under the curve.
  
  return(AUC)
  
}

nMDS_results <- lapply(tentative_dimensions, optimim_nMDS)
plot(tentative_dimensions, (-1)*unlist(nMDS_results))


#c. LLE Method

optimim_LLE <- function(dim){ #Builds a wrapper for optimizing LLE
  
  #Defines trajectorz matrix inside function.
  
  trajectory_matrix <- as.data.frame(embedding(unlist(data_fail), 50, 1))
  
  #1.2 Transforms trajectory matrix into a dimRed data object.
  
  trajectory_matrix <- dimRedData(trajectory_matrix)
  
  #2. 
  real_dim <- round(dim,0)
  LLE_result <- embed(trajectory_matrix, .method = 'LLE', ndim = real_dim)
  
  AUC <- (-1)*quality(LLE_result, .method = "mean_R_NX")
  
  #Returns the area under the curve.
  
  return(AUC)
  
}
LLE_results <- lapply(tentative_dimensions, optimim_LLE)
plot(tentative_dimensions, (-1)*unlist(LLE_results))

#d.  Laplacian Eigenmaps

optimim_LaplacEigen <- function(dim){ #Builds a wrapper for optimizing LLE
  
  #Defines trajectorz matrix inside function.
  
  trajectory_matrix <- as.data.frame(embedding(unlist(data_fail), 50, 1))
  
  #1.2 Transforms trajectory matrix into a dimRed data object.
  
  trajectory_matrix <- dimRedData(trajectory_matrix)
  
  #2. 
  real_dim <- round(dim,0)
  Laplac_result <- embed(trajectory_matrix, .method = "LaplacianEigenmaps", ndim = real_dim)
  
  AUC <- (-1)*quality(Laplac_result, .method = "mean_R_NX")
  
  #Returns the area under the curve.
  
  return(AUC)
  
}
Laplac_results <- lapply(tentative_dimensions, optimim_LaplacEigen)
plot(tentative_dimensions, (-1)*unlist(Laplac_results))

#f. Isomap

optimim_isomap <- function(dim){ #Builds a wrapper for optimizing LLE
  
  #Defines trajectorz matrix inside function.
  
  trajectory_matrix <- as.data.frame(embedding(unlist(data_fail), 50, 1))
  
  #1.2 Transforms trajectory matrix into a dimRed data object.
  
  trajectory_matrix <- dimRedData(trajectory_matrix)
  
  #2. 
  real_dim <- round(dim,0)
  isomap_result <- embed(trajectory_matrix, .method = "Isomap", ndim = real_dim)
  
  AUC <- (-1)*quality(isomap_result, .method = "mean_R_NX")
  
  #Returns the area under the curve.
  
  return(AUC)
  
}
isomap_results <- lapply(tentative_dimensions, optimim_isomap)
plot(tentative_dimensions, (-1)*unlist(isomap_results))

#g. tSNE


#optimim_tsne <- function(dim){ #Builds a wrapper for optimizing LLE
  
  #Defines trajectorz matrix inside function.
  
  trajectory_matrix <- as.data.frame(embedding(unlist(data_fail), 50, 1))
  
  #1.2 Transforms trajectory matrix into a dimRed data object.
  
  trajectory_matrix <- dimRedData(trajectory_matrix)
  
  #2. 
  real_dim <- round(dim,0)
  tsne_result <- embed(trajectory_matrix, .method = "tSNE", ndim = real_dim)
  
  AUC <- (-1)*quality(tsne_result, .method = "mean_R_NX")
  
  #Returns the area under the curve.
  
  return(AUC)
  

#tsne_results <- lapply(tentative_dimensions, optimim_tsne)
#plot(tentative_dimensions, (-1)*unlist(tsne_results))

tSNE_red <- embed(trajectory_matrix, .method = 'tSNE', ndim = 3)
Rnx_tSNE <- quality(tSNE_red, .method = 'mean_R_NX')
Rnx_tSNE

#h. PCA

optimim_pca <- function(dim){ #Builds a wrapper for optimizing LLE
  
  #Defines trajectorz matrix inside function.
  
  trajectory_matrix <- as.data.frame(embedding(unlist(data_fail), 50, 1))
  
  #1.2 Transforms trajectory matrix into a dimRed data object.
  
  trajectory_matrix <- dimRedData(trajectory_matrix)
  
  #2. 
  real_dim <- round(dim,0)
  pca_result <- embed(trajectory_matrix, .method = "PCA", ndim = real_dim)
  
  AUC <- (-1)*quality(pca_result, .method = "mean_R_NX")
  
  #Returns the area under the curve.
  
  return(AUC)
  
}
pca_results <- lapply(tentative_dimensions, optimim_pca)

#i. Standard Method = False Nearest Neighboors

standard_embedding <- embed_udf(data_fail)
standard_trajectory_matrix <- standard_embedding[[4]]

tSNE_original <- embed(dimRedData(standard_trajectory_matrix), .method = "tSNE", ndim = 3)
Rnx_tSNE_original <- quality(tSNE_original, .method = 'mean_R_NX')
Rnx_tSNE_original

#2. Prepares data frame and plots results to compare.

library(ggplot2)
library(reshape2)
RNX_dataframe <- data.frame(kPCA = (-1)*unlist(kPCA_results),
                            nMDS = (-1)*unlist(nMDS_results),
                            LLE = (-1)*unlist(LLE_results),
                            Laplac.Eigen = (-1)*unlist(Laplac_results),
                            Isomap = (-1)*unlist(isomap_results),
                            PCA = (-1)*unlist(pca_results) ,
                            DIM = tentative_dimensions)

melted_RNX <- melt(RNX_dataframe, id = "DIM")

ggplot(data = melted_RNX, aes(x= DIM,y = value, colour = variable)) + geom_line(size = 0.75)


#3. Creates 3D Plot for each method

#a. Isomap

Isomap_red <- embed(trajectory_matrix, .method = 'Isomap', ndim = 3)
x1 <- as.data.frame(unlist(Isomap_red@data))
x1_new <- x1[2:nrow(x1), ]
x1 <- x1[(1:nrow(x1) - 1), ]

isomap_attractor <- plot3D::border3D(x1$iso1, x1$iso2, x1$iso3, x1_new$iso1, x1_new$iso2, x1_new$iso3)


#b. PCA

PCA_red <- embed(trajectory_matrix, .method = 'PCA', ndim = 3)
x2 <- as.data.frame(unlist(PCA_red@data))
x2_new <- x2[2:nrow(x2), ]
x2 <- x2[(1:nrow(x2) - 1), ]

pca_attractor <- plot3D::border3D(x2$PC1, x2$PC2, x2$PC3, x2_new$PC1, x2_new$PC2, x2_new$PC3)

#c. nMDS

nMDS_red <- embed(trajectory_matrix, .method = 'nMDS', ndim = 3)
x3 <- as.data.frame(unlist(nMDS_red@data))
x3_new <- x3[2:nrow(x3), ]
x3 <- x3[(1:nrow(x3) - 1), ]

nMDS_attractor <- plot3D::border3D(x3$NMDS1, x3$NMDS2, x3$NMDS3, x3_new$NMDS1, x3_new$NMDS2, x3_new$NMDS3)

#d. KPCA

kPCA_red <- embed(trajectory_matrix, .method = 'kPCA', ndim = 3)
x4 <- as.data.frame(unlist(kPCA_red@data))
x4_new <- x4[2:nrow(x4), ]
x4 <- x4[(1:nrow(x4) - 1), ]

kpca_attractor <- plot3D::border3D(x4$kPCA1, x4$kPCA2, x4$kPCA3, x4_new$kPCA1, x4_new$kPCA2, x4_new$kPCA3)

#e. Plot Standard tSNE

x5 <- as.data.frame(unlist(tSNE_original@data))
x5_new <- x5[2:nrow(x5), ]
x5 <- x5[(1:nrow(x5) - 1), ]

original_attractor <- plot3D::border3D(x5$tSNE1, x5$tSNE2, x5$tSNE3, x5_new$tSNE1, x5_new$tSNE2, x5_new$tSNE3)


#4. Plots all atractors in same window.

par(mfrow = c(1,1))


isomap_attractor <- plot3D::border3D(x1$iso1, x1$iso2, x1$iso3, x1_new$iso1, x1_new$iso2, x1_new$iso3, main = 'Isomap')
pca_attractor <- plot3D::box3D(x2$PC1, x2$PC2, x2$PC3, x2_new$PC1, x2_new$PC2, x2_new$PC3, main = ' PCA')
nMDS_attractor <- plot3D::border3D(x3$NMDS1, x3$NMDS2, x3$NMDS3, x3_new$NMDS1, x3_new$NMDS2, x3_new$NMDS3, main = 'nMDS')
kpca_attractor <- plot3D::border3D(x4$kPCA1, x4$kPCA2, x4$kPCA3, x4_new$kPCA1, x4_new$kPCA2, x4_new$kPCA3, main = 'kPCA')
original_attractor <- plot3D::border3D(x5$tSNE1, x5$tSNE2, x5$tSNE3, x5_new$tSNE1, x5_new$tSNE2, x5_new$tSNE3)

#5. Trajectory Matrix


official_trajectory_matrix <-  embed(trajectory_matrix, .method = 'PCA', ndim = 10)
official_trajectory_matrix <- as.data.frame(unlist(official_trajectory_matrix@data))

official_trajectory_matrix # Matriz a ser utilizada em análises futuras.
