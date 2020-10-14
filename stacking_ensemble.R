
library(proxy)
library(colorspace)

########################################ensemble with native functions #################################
#Modelos utilizados (já rodados anteriormente)
# model.hc <- hclust(distance_trajectory)
# model.kmeans <- kmeans(official_trajectory_matrix, 3, nstart = 20)#5 centróides, o r irá tentar  20 começos diferentes
# som.model <- som(data.som, grid=som.grid, rlen=7000, alpha=c(0.05,0.01), keep.data =TRUE)

#Stacking, using the hard-voting strategy
master.cluster <-model.kmeans$cluster
slave.hc <- hc.cut
slave.som.kmeans <- decoded.kmeans
slave.som.hc <- decoded.hc

# Preparing the stacked clustering
stacked.clustering <- rep(NA, length(master.cluster))
names(stacked.clustering) <- 1:length(master.cluster)

for (cluster in unique (master.cluster)) {
  indexes = which(master.cluster == cluster, arr.ind =T )
  slave1.votes <- table(slave.hc[indexes])
  slave1.maxcount <- names(slave1.votes)[which.max(slave1.votes)]
 
  
   slave1.indexes = which(slave.hc == slave1.maxcount,
arr.ind = T)
  slave2.votes <- table(slave.som.hc[indexes])
  slave2.maxcount <- names(slave2.votes)[which.max(slave2.votes)]
 
  
   slave2.indexes = which(slave.som.hc == slave2.maxcount,
arr.ind = T)
  slave3.votes <- table(slave.som.kmeans[indexes])
  slave3.maxcount <- names(slave3.votes)[which.max(slave3.votes)]

 
  stacked.clustering[indexes] <- slave3.maxcount
}


#Ploting results
points <- cmdscale(distance_trajectory) # Running the PCA, reducing to 2 dimensions
palette <- colorspace::diverge_hcl(12) # Creating a color
palette

windows()
plot(points, main="K-Means Clustering",
    col= as.factor(master.cluster),
    mai= c(0,0,0,0), mar = c(0,0,0,0),
    xaxt = "n", yaxt = "n", xlab="", ylab="")

windows()
par(mfrow=c(2,2), mar = rep(1.5,4))
plot(points, main="Hierarchical Clustering", col=
       as.factor(slave.hc),
     mai= c(0,0,0,0), mar = c(0,0,0,0),
     xaxt = "n", yaxt = "n", xlab="", ylab="")


plot(points, main="Som.k Clustering", col=
       as.factor(slave.som.kmeans),
     mai= c(0,0,0,0), mar = c(0,0,0,0),
     xaxt = "n", yaxt = "n", xlab="", ylab="")


plot(points, main="Som.hc Clustering", col=
       as.factor(slave.som.hc),
     mai= c(0,0,0,0), mar = c(0,0,0,0),
     xaxt = "n", yaxt = "n", xlab="", ylab="")

plot(points, main="Stacked Clustering", col=
       as.factor(stacked.clustering),
     mai= c(0,0,0,0), mar = c(0,0,0,0), # mar (for margin!), specifying the margins in inches using the mai argument
     xaxt = "n", yaxt = "n", xlab="", ylab="") # xaxt="n" and yaxt="n" suppress the x and y axis respectively



#Internal validation
library(cluster)
D<- daisy(official_trajectory_matrix)
stacked <- as.integer(stacked.clustering)
unname(stacked, force=F)
windows()
plot(silhouette(stacked, D), col =1, border = NA)

#External Validation
library(fpc)
ext.val <- as.numeric(stacked)
clus.stats <- cluster.stats(d= dist(official_trajectory_matrix),
                            ext.val, stacked)
clus.stats$corrected.rand
##[1] 1
clus.stats$vi
##[1] 0