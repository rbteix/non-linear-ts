#Clustering using ML tequiniches using the embedding trajectory matrix
data <- official_trajectory_matrix
windows()
plot(official_trajectory_matrix)
#Finding the ideal number of clusters
#The NbClust package provides 30 indices to determine the number of clusters in a dataset
library(NbClust)
nb_data <- NbClust(data, distance = "euclidean", min.nc = 2,
                   max.nc = 10, method = "kmeans")

##******************************************************************* 
#  * Among all indices:                                                
#  * 10 proposed 2 as the best number of clusters 
#  * 9 proposed 3 as the best number of clusters 
#  * 2 proposed 5 as the best number of clusters 
#  * 2 proposed 10 as the best number of clusters 

##***** Conclusion *****                            

##  * According to the majority rule, the best number of clusters is  2 
##*******************************************************************

windows()
hist(nb_data$Best.nc[1,], breaks = max(na.omit(nb_data$Best.nc[1,])))

##############################################AGRUPAMENTO HIER�RQUICO####################################
library(stats)
library(cluster)

distance_trajectory <- dist((official_trajectory_matrix), diag = TRUE)
model.hc <- hclust(distance_trajectory)
windows()
plot(model.hc, cex = 0.6,  hang=-1)
# Split in 5 clusters
hc.cut <- cutree(model.hc, 2)
rect.hclust(model.hc, k=2, border=2:3)
table(hc.cut)
############################
##hc.cut
##1    2 
##1392  972 
#############################


#Valida��o
library("factoextra")

# can change the clustering algorithm in FUN (ex kmeans)
windows()
fviz_nbclust(as.matrix(unlist(model.hc$height)), FUN = hcut, method = "silhouette", ylim =0.8)+geom_vline(xintercept = F,linetype = 2)

windows()
fviz_nbclust(as.matrix(unlist(model.hc$height)), FUN = kmeans, method = "silhouette")+geom_vline(xintercept = F, linetype = 2)

library("cluster")
sil.hc <- silhouette((hc.cut), dist(official_trajectory_matrix))
windows()
plot(sil.hc)

library(factoextra)
windows()
fviz_silhouette(sil.hc)


clus_hc <- hc.cut
############################################### MAPAS AUTO-ORGANIZ�VEIS (SOM) ##############################
library(kohonen)
library(ggplot2)
library(reshape2)


set.seed(382)

#Defini��o cores paleta
pretty_palette <- c("#1f77b4", '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2')

#Costumizando a paleta para o pacote kohonen
coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6, alpha=alpha)[n:1]}
#Em kohonen as linhas dos dados s�o amostras e as colunas vari�veis
data.som <- as.matrix(official_trajectory_matrix)

#Tornando o data set uma matriz
#Centralizar e dimensionar todas as vari�veis para dar-lhes igual import�ncia durante
# o processo de treinamento do SOM (scale).

#Criar a Grade SOM  geralmente se especifica o tamanho da grade de treinamento antes de treinar.
#Hexagonal e Circular,s�o tpologias poss�veis.

som.grid <- somgrid(xdim = 17, ydim=15, topo = "hexagonal")

#Modelo: op��es para o n�mero de itera��es (rlen), as taxas de aprendizagem(alpha)
som.model <- som(data.som, grid=som.grid, rlen=7000, alpha=c(0.05,0.01), keep.data =TRUE)

#Gr�ficos, visualiza��o
#1 Progresso de Treinamento:
#� medida que as itera��es de treinamento do SOM progridem, a dist�ncia dos pesos de cada n� �s amostras
#representadas por esse n� � reduzida. Idealmente, essa dist�ncia deve atingir um patamar m�nimo.
#Esta op��o de plotagem mostra o progresso ao longo do tempo. Se a curva estiver diminuindo continuamente,
#mais itera��es ser�o necess�rias.
windows()
plot(som.model, type="changes")               

#Permite visualizar a contagem de quantas amostras s�o mapeadas para cada n�. Essa m�trica pode ser usada 
#como uma medida da qualidade do mapa - idealmente, a distribui��o da amostra � relativamente uniforme.
#N�s vazios indicam que o tamanho do seu mapa � muito grande para o n�mero de amostras. 
windows()
plot(som.model, type= "count", main="Node Counts", palette.name=coolBlueHotRed) 

#Mapa de qualidade, um bom resultado o mostra quase uniforme
windows()
plot(som.model, type = "quality", main="Node Quality/Distance", palette.name=coolBlueHotRed)

#Permite visualizar a dist�ncia entre cada n� e seus vizinhos. As �reas de baixa dist�ncia vizinha indicam grupos de n�s semelhantes. 
#�reas com grandes dist�ncias indicam que os n�s s�o muito mais dissimilares                
#windows()
#plot(som.model, type="dist.neighbours", main = "SOM neighbour distances", palette.name=grey.colors)

#Os vetores de peso do n�, ou "c�digos", s�o compostos de valores normalizados das vari�veis originais usadas para gerar o SOM. 
#O vetor de pesos de cada n� � representativo,semelhante das amostras mapeadas para esse n�. Ao visualizar os vetores de peso no mapa, 
#pode-se ver padr�es na distribui��o de amostras e vari�veis
#windows()
#plot(som.model, type="codes")


# Plotar heatmap 
windows()
plot(som.model, type = "property", property = som.model$codes[[1]], main= "Fails", palette.name=coolBlueHotRed)

# Clustering os resultados 

# Mostra o WCSS m�trica para k-means para diferentes tamanhos de clustering, confirmar se entre 3 e 5 tamb�m podem
#ser tamanhos de clusters ideais
mydata <-unlist(som.model$codes)
mydata <- as.matrix(mydata)
wss <- (nrow(mydata)-1)*sum(apply(mydata,2, var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
windows()
par(mar=c(5.1,4.1,4.1,2.1))
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares", main="Within cluster sum of squares (WCSS)")


# Formando grade de clusters
## Usando HC do agrupar codebook vectors
som.cluster <- cutree(hclust(dist(mydata)), 3)

# Mostrar gr�fico com clusters em diferentes cores	
windows()
plot(som.model, type="mapping", bgcol = pretty_palette[som.cluster], main = "Clusters"  )
add.cluster.boundaries(som.model, som.cluster)

#Mostrar o mesmo gr�fico adicionado os c�digos
windows()
plot(som.model, type="codes", bgcol = pretty_palette[som.cluster], main = "Clusters")
add.cluster.boundaries(som.model, som.cluster)

#Valida��o
library("factoextra")

# can change the clustering algorithm in FUN (ex kmeans)
windows()
fviz_nbclust(as.matrix(unlist(som.model$codes)), FUN = hcut, method = "silhouette", ylim =0.8)+geom_vline(xintercept = 3,linetype = 2)

windows()
fviz_nbclust(as.matrix(unlist(som.model$codes)), FUN = kmeans, method = "silhouette")+geom_vline(xintercept = 3, linetype = 2)

library("cluster")
sil.som1 <- silhouette((decoded_hc_clus), dist(official_trajectory_matrix))
windows()
plot(sil.som1)

library(factoextra)
windows()
fviz_silhouette(sil.som1)

# Stacking Hierarchical CLustering e SOM
#som.cluster <- cutree(hclust(dist(as.vector(unlist(som.model$codes)))), 3)

# Stacking K-Means and SOM

similarity.vec <- as.vector(unlist(som.model$codes))
kmeans.model.som <- kmeans(similarity.vec, 3, iter.max = 100)


#Decodifica��o dados reais
decoded.hc<- som.cluster[som.model$unit.classif]
decoded_hc_clus <- decoded.hc
decoded.kmeans <- kmeans.model.som$cluster[som.model$unit.classif]
decoded_kmeans_clus <- decoded.kmeans

#############################################Fuzzy c-means clustering################
set.seed(333)


# Compute fuzzy clustering
library(e1071)
data <- official_trajectory_matrix
data_cm <- cmeans(data, 2)

library(factoextra)
windows()
fviz_cluster(list(data = data, cluster=data_cm$cluster), 
             ellipse.type = "norm",
             ellipse.level = 0.68,
             palette = "jco",
             ggtheme = theme_minimal())

clus_cm <- data_cm$cluster

Valida��o
library("factoextra")

# can change the clustering algorithm in FUN (ex kmeans)
windows()
fviz_nbclust(as.matrix(unlist(data_cm$membership)), FUN = hcut, method = "silhouette", ylim =0.8)+geom_vline(xintercept = F,linetype = 2)

windows()
fviz_nbclust(as.matrix(unlist(data_cm$membership)), FUN = kmeans, method = "silhouette")+geom_vline(xintercept = F, linetype = 2)

library("cluster")
sil.cm <- silhouette((clus_cm), dist(official_trajectory_matrix))
windows()
plot(sil.cm)

library(factoextra)
windows()
fviz_silhouette(sil.cm)
#

##############################################DBSCAN#####################################

library(dbscan) 
data <- official_trajectory_matrix
dbsc.trajectory.matrix <- dbscan(official_trajectory_matrix, eps = 0.14, minPts = 6)
# Plot DBSCAN results
windows()
fviz_cluster(dbsc.trajectory.matrix, data = data, stand = FALSE,
             ellipse = F, show.clust.cent = F,
             geom = "point",palette = "jco", ggtheme = theme_classic())


######################################### Ensemble ######################################

#Juntando clusters em um dataset
data.ens <- official_trajectory_matrix

#HMM <- (hmm.fit@posterior$state)
#data.ens["HMM"] <- c(HMM)


SOM1 <- decoded.kmeans
data.ens["SOM1"] <- c(SOM1)

SOM2 <- decoded.hc 
data.ens["SOM2"] <- c(SOM2)

HC <- hc.cut
data.ens["HC"] <- c(HC)

CM <- clus_cm
data.ens["CM"] <- c(CM)

#Ordenar colunas
#data.ens= data.ens[,c(2,1,3,4, 5,6)]

#Plotar gr�fico de cada modelo e comparar

#HMM
#data.matrix.hmm <- cbind(data.matrix$Values, hmm.fit5@posterior$state)
#data.matrix.hmm = as.data.frame(data.matrix.hmm)
#data.matrix.hmm["Date"] = c(Data.1)
#colnames(data.matrix.hmm) = c("Data", "HMM", "Date")
#data.matrix.hmm = data.matrix.hmm[,c(3,1,2)]

#Geometric line
#windows()
#hml <-ggplot( data.matrix.hmm, aes(Date, Data, color =( data.matrix.hmm$HMM))) +geom_line()  
#hml + scale_color_gradientn(colours = rainbow(5))

#Geometric points
#windows()
#hmg <- ggplot( data.matrix.hmm, aes(Date, Data, color = data.matrix.hmm$HMM)) +geom_point()
#hmg + scale_color_gradientn(colours = rainbow(5))

#SOM
#Geometric line
#trajectory.som1 <- cbind(official_trajectory_matrix, decoded.kmeans)
#trajectory.som1  = data.frame(trajectory.som1)
#windows()
#plot(trajectory.som1)
#data.matrix.som1["Date"] = c(Data.1)
#colnames(trajectory.som1)= c("Data", "SOM1")
#data.matrix.som1 = data.matrix.som1[c(3,1,2)]

#Geometric line
#windows()
#soml1 <-ggplot( trajectory.som1, aes(Date, Data, color =(traject$SOM1))) +geom_line()  
#soml1 + scale_color_gradientn(colours = rainbow(5))

#Geometric points
#windows()
#somg1 <- ggplot( data.matrix.som1, aes(Date, Data, color = data.matrix.som1$SOM1)) +geom_point()
#somg1 + scale_color_gradientn(colours = rainbow(5))

#SOM2
#data.matrix.som2 <- cbind(data.matrix$Values, decoded.hc)
#data.matrix.som2  = data.frame(data.matrix.som2)
#data.matrix.som2["Date"] = c(Data.1)
#colnames(data.matrix.som2)= c("Data","SOM2", "Date")
#data.matrix.som2 = data.matrix.som2[c(3,1,2)]

#Geometric line
#windows()
#soml2 <-ggplot( data.matrix.som2, aes(Date, Data, color =( data.matrix.som2$SOM2))) +geom_line()  
#soml2 + scale_color_gradientn(colours = rainbow(5))

#Geometric pints
#windows()
#somg2 <- ggplot( data.matrix.som2, aes(Date, Data, color = data.matrix.som2$SOM2)) +geom_point()
#somg2 + scale_color_gradientn(colours = rainbow(5))


#HC
#data.matrix.hc <- cbind(data.matrix$Values,  hc.cut)
#data.matrix.hc = as.data.frame(data.matrix.hc)
#data.matrix.hc["Date"] = c(Data.1)
#colnames(data.matrix.hc) =  c("Data",  "HC", "Date")
#data.matrix.hc = data.matrix.hc[,c(3,1,2)]

#Geometric line
#windows()
#hcl <-ggplot( data.matrix.hc, aes(Date, Data, color =( data.matrix.hc$HC))) +geom_line()  
#hcl + scale_color_gradientn(colours = rainbow(5))

#Geometric points
#windows()
#hcg <- ggplot( data.matrix.hc, aes(Date, Data, color = data.matrix.hc$HC)) +geom_point()
#hcg + scale_color_gradientn(colours = rainbow(5))


#Matrix de Correla��o dos modelos 2 a 2
#cor.hc.km <- cor(cbind(HC, KMEANS))
#cor.hc.km
#cor.hc.hm <- cor(cbind(HC, HMM))
#cor.hc.hm
#cor.km.hm <- cor(cbind(HMM, KMEANS))
#cor.km.hm
#cor.som1.km <- cor(cbind(SOM1, KMEANS))
#cor.som1.km
#cor.som1.hmm <- cor(cbind(SOM1, HMM))
#cor.som1.hmm
#cor.som2.km <- cor(cbind(SOM2, KMEANS))
#cor.som2.km
#cor.som2.hmm <- cor(cbind(SOM2, HMM))
#cor.som2.hmm

#cor.som1.hc <- cor(cbind(SOM1, HC))
#cor.som1.hc

#cor.som2.hc <- cor(cbind(SOM2, HC))
#cor.som2.hc

cor_mat <- cor(cbind(SOM1, SOM2, HC, CM))
##############################################ENSEMBLE - KM�DIAS##########################
#Utilizando k-m�dias para fazer um unico agrupamento dos 4 modelos

#Elbow m�todo, para identificar o n�mero ideal de clusters para o dataset
k.max <- 15
data <- data.ens[,-12]
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=20,iter.max = 15 )$tot.withinss})
wss
windows()
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

# Elbow method for kmeans
fviz_nbclust(data.ens[,-12], kmeans, method = "wss") +
  geom_vline(xintercept = F, linetype = 2)

library(NbClust)
nb_data <- NbClust(data.ens[,-12], distance = "euclidean", min.nc = 2,
                   max.nc = 10, method = "kmeans")
#******************************************************************* 
#  * Among all indices:                                                
#  * 4 proposed 2 as the best number of clusters 
#* 13 proposed 3 as the best number of clusters 
#* 2 proposed 4 as the best number of clusters 
#* 1 proposed 5 as the best number of clusters 
#* 1 proposed 7 as the best number of clusters 
#* 1 proposed 9 as the best number of clusters 
#* 1 proposed 10 as the best number of clusters 
#
#***** Conclusion *****                            
  
#  * According to the majority rule, the best number of clusters is  3 

#******************************************************************* 


#kmeans ensemble model
library(stats)
set.seed(5555)

#Criar modelo
ens.model.kmeans3 <- kmeans(data.ens[,-12], 3, nstart = 20)#7 centr�ides, o r ir� tentar  20 come�os diferentes

windows()
fviz_cluster(list(data = data.ens[,-12], cluster=ens.model.kmeans3$cluster), 
             ellipse.type = "norm",
             ellipse.level = 0.68,
             palette = "jco",
             ggtheme = theme_minimal(),
             geom = "point")


library(cluster) 
windows()
clusplot(data.ens[,-12], ens.model.kmeans3$cluster, color=TRUE, shade=TRUE, 
         labels=1, lines=0)


#Plotar gr�ficos
#library(cluster)
#library(fpc)

#Geometric line
#windows()
#ens.plot <- ggplot( data.ens, aes(Date, Values, color = ens.model.kmeans$cluster)) +geom_line()
#ens.plot + scale_color_gradientn(colours = rainbow(5))

#Geometric points
#windows()
#ens.plot <- ggplot( data.ens, aes(Date, Values, color = ens.model.kmeans$cluster)) +geom_point()
#ens.plot + scale_color_gradientn(colours = rainbow(5))

data.final1 <- as.data.frame (ens.model.kmeans1$cluster)
data.final1["PC1"] <- c(official_trajectory_matrix[,1])
data.final1["PC2"] <- c(official_trajectory_matrix[,2])
data.final1["PC3"] <- c(official_trajectory_matrix[,3])
data.final1["PC4"] <- c(official_trajectory_matrix[,4])
data.final1["PC5"] <- c(official_trajectory_matrix[,5])
data.final1["PC6"] <- c(official_trajectory_matrix[,6])
data.final1["PC7"] <- c(official_trajectory_matrix[,7])
data.final1["PC8"] <- c(official_trajectory_matrix[,8])
data.final1["PC9"] <- c(official_trajectory_matrix[,9])
data.final1["PC10"] <- c(official_trajectory_matrix[,10])

#data.final["Date"]= c(Date)
#colnames(data.final) = c("Clusters", "Data")
#data.final = data.final[,c(2,1)]

###################################### Valida��o ###############################
library("clValid")
# comparing hierarchical and kmeans with internal validation indexes
clmethods <- c("hierarchical","kmeans")
intern <- clValid(data.ens[,-12], nClust = 2:13, maxitems = 2364, 
                  clMethods = clmethods, validation = "internal", method = "complete")
summary(intern)
################################################################################
#Validation Measures:
#                                2        3        4        5        6        7        8        9

#hierarchical Connectivity    0.0000   1.5357   1.5357   1.5357   1.5357   1.5357   3.8738 466.2706
#Dunn                         0.5245   0.6162   0.7171   0.7974   0.8219   0.9009   1.1868   0.0721
#Silhouette                   0.5956   0.5279   0.5635   0.6064   0.6218   0.6432   0.6458   0.3765
#kmeans       Connectivity    2.3381   1.5357   1.5357   1.5357   1.5357   1.5357   3.8738 183.9591
#Dunn                         0.5838   0.5833   0.7171   0.7974   0.8219   0.9009   1.1868   0.0462
#Silhouette                   0.6362   0.5396   0.5635   0.6064   0.6218   0.6432   0.6458   0.4274

#Optimal Scores:
  
#  Score  Method       Clusters
#Connectivity 0.0000 hierarchical 2       
#Dunn         1.1868 hierarchical 8       
#Silhouette   0.6458 hierarchical 8       
######################################################################################


library("factoextra")
windows()
fviz_nbclust(data.ens[,-12], FUN = hcut, method = "silhouette", ylim =0.8)+geom_vline(xintercept = F,linetype = 2)


windows()
fviz_nbclust(data.ens[,-12], FUN = kmeans, method = "silhouette")+geom_vline(xintercept = F, linetype = 2)


library("cluster")
sil.kmeans3 <- silhouette((ens.model.kmeans3$cluster), dist(data.ens[,-12]))
windows()
plot(sil.kmeans3)

library(factoextra)
windows()
fviz_silhouette(sil.kmeans3)




##################################################### mean and sd of each cluster  #########################
clus1 <- as.matrix(official_trajectory_matrix)[which(ens.model.kmeans$cluster ==1),]
mean(clus1)
sd(clus1)

clus2 <- as.matrix(official_trajectory_matrix)[which(ens.model.kmeans$cluster ==2),]
mean(clus2)
sd(clus2)

clus3 <- as.matrix(official_trajectory_matrix)[which(ens.model.kmeans$cluster ==3),]
mean(clus3)
sd(clus3)

clus4 <- as.matrix(official_trajectory_matrix)[which(ens.model.kmeans$cluster ==4),]
mean(clus4)
sd(clus4)

clus5 <- as.matrix(official_trajectory_matrix)[which(ens.model.kmeans$cluster ==5),]
mean(clus5)
sd(clus5)

###################################Time seires from clusters###################

#Step 3: Time series reconstruction (diagonal averaging)
#Step 3a: User-defined function for averaging of minor diagonals
diag.ave<-function(mat) {
  nrow <- nrow(mat)
  ncol <- ncol(mat)
  hold<-matrix(0,(nrow+(ncol-1)))
  for(i in 1:(nrow+(ncol-1))) {
    if(i==1) {d<-mat[1,1]}
    if(i>1 & i<=ncol) {d<-diag(mat[i:1,1:i])}
    if(i>ncol & i<=nrow) {d<-diag(mat[i:(i-(ncol-1)),1:ncol])}
    if(i>nrow & i<(nrow+(ncol-1))) {
      d<-diag(mat[nrow:(i-(ncol-1)),(i-(nrow-1)):ncol])}
    if(i==(nrow+(ncol-1))) {d<-mat[nrow,ncol]}
    d.ave<-mean(d) 
    hold[i,]<-d.ave
    #average minor diagonals
  } #end loop
  return(hold)
} #end function
#Step 3b: Reconstructed time series
clus1<-diag.ave(clus1)
clus2<-diag.ave(clus2)
clus3<-diag.ave(clus3)
clus4<-diag.ave(clus4)
clus5<-diag.ave(clus5)
traj_matrix <- diag.ave(as.matrix(official_trajectory_matrix))
windows()
plot(traj_matrix, type = "l")