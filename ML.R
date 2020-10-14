#Clustering using ML tequiniches using the embedding trajectory matrix

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

##############################################AGRUPAMENTO HIERÁRQUICO####################################
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


############################################### MAPAS AUTO-ORGANIZÁVEIS (SOM) ##############################
library(kohonen)
library(dummies)
library(ggplot2)
library(sp)
library(reshape2)
library(rgeos)

set.seed(382)

#Definição cores paleta
pretty_palette <- c("#1f77b4", '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2')

#Costumizando a paleta para o pacote kohonen
coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6, alpha=alpha)[n:1]}
#Em kohonen as linhas dos dados são amostras e as colunas variáveis
data.som <- as.matrix(official_trajectory_matrix)

#Tornando o data set uma matriz
#Centralizar e dimensionar todas as variáveis para dar-lhes igual importância durante
# o processo de treinamento do SOM (scale).

#Criar a Grade SOM  geralmente se especifica o tamanho da grade de treinamento antes de treinar.
#Hexagonal e Circular,são tpologias possíveis.

som.grid <- somgrid(xdim = 20, ydim=15, topo = "hexagonal")

#Modelo: opções para o número de iterações (rlen), as taxas de aprendizagem(alpha)
som.model <- som(data.som, grid=som.grid, rlen=7000, alpha=c(0.05,0.01), keep.data =TRUE)

#Gráficos, visualização
#1 Progresso de Treinamento:
#À medida que as iterações de treinamento do SOM progridem, a distância dos pesos de cada nó às amostras
#representadas por esse nó é reduzida. Idealmente, essa distância deve atingir um patamar mínimo.
#Esta opção de plotagem mostra o progresso ao longo do tempo. Se a curva estiver diminuindo continuamente,
#mais iterações serão necessárias.
windows()
plot(som.model, type="changes")               

#Permite visualizar a contagem de quantas amostras são mapeadas para cada nó. Essa métrica pode ser usada 
#como uma medida da qualidade do mapa - idealmente, a distribuição da amostra é relativamente uniforme.
#Nós vazios indicam que o tamanho do seu mapa é muito grande para o número de amostras. 
windows()
plot(som.model, type= "count", main="Node Counts", palette.name=coolBlueHotRed) 

#Mapa de qualidade, um bom resultado o mostra quase uniforme
windows()
plot(som.model, type = "quality", main="Node Quality/Distance", palette.name=coolBlueHotRed)

#Permite visualizar a distância entre cada nó e seus vizinhos. As áreas de baixa distância vizinha indicam grupos de nós semelhantes. 
#Áreas com grandes distâncias indicam que os nós são muito mais dissimilares                
#windows()
#plot(som.model, type="dist.neighbours", main = "SOM neighbour distances", palette.name=grey.colors)

#Os vetores de peso do nó, ou "códigos", são compostos de valores normalizados das variáveis originais usadas para gerar o SOM. 
#O vetor de pesos de cada nó é representativo,semelhante das amostras mapeadas para esse nó. Ao visualizar os vetores de peso no mapa, 
#pode-se ver padrões na distribuição de amostras e variáveis
#windows()
#plot(som.model, type="codes")


# Plotar heatmap com as variáveis scaled/normalizadas
windows()
plot(som.model, type = "property", property = som.model$codes[[1]], main= "Fails", palette.name=coolBlueHotRed)

# Clustering os resultados 

# Mostra o WCSS métrica para k-means para diferentes tamanhos de clustering, confirmar se entre 3 e 5 também podem
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
som.cluster <- cutree(hclust(dist(mydata)), 5)

# Mostrar gráfico com clusters em diferentes cores	
windows()
plot(som.model, type="mapping", bgcol = pretty_palette[som.cluster], main = "Clusters"  )
add.cluster.boundaries(som.model, som.cluster)

#Mostrar o mesmo gráfico adicionado os códigos
windows()
plot(som.model, type="codes", bgcol = pretty_palette[som.cluster], main = "Clusters")
add.cluster.boundaries(som.model, som.cluster)
#Cinco clusters mostrou-se uma boa quantidade de clusters

# Stacking Hierarchical CLustering e SOM
som.cluster <- cutree(hclust(dist(as.vector(unlist(som.model$codes)))), 5)

# Stacking K-Means and SOM

similarity.vec <- as.vector(unlist(som.model$codes))
kmeans.model.som <- kmeans(similarity.vec, 5, iter.max = 100)


#Decodificação dados reais
decoded.hc<- som.cluster[som.model$unit.classif]
decoded.kmeans <- kmeans.model.som$cluster[som.model$unit.classif]

#############################################Fuzzy c-means clustering################
set.seed(333)


# Compute fuzzy clustering
library(e1071)
data_cm <- cmeans(official_trajectory_matrix, 2)
library(factoextra)
windows()
fviz_cluster(list(data = official_trajectory_matrix, cluster=data_cm$cluster), 
             ellipse.type = "norm",
             ellipse.level = 0.68,
             palette = "jco",
             ggtheme = theme_minimal(), geom = "point")

data_clus <- data_cm$cluster

##############################################DBSCAN#####################################

library(dbscan) 
dbsc.trajectory.matrix <- dbscan(official_trajectory_matrix, eps = 0.14, minPts = 6)
# Plot DBSCAN results
windows()
fviz_cluster(dbsc.trajectory.matrix, data = official_trajectory_matrix, stand = FALSE,
             ellipse = T, show.clust.cent = F,
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

CM <- data_clus
data.ens["CM"] <- c(CM)

#Ordenar colunas
#data.ens= data.ens[,c(2,1,3,4, 5,6)]

#Plotar gráfico de cada modelo e comparar

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


#Matrix de Correlação dos modelos 2 a 2
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

cor.ens <- cor(cbind(SOM1, SOM2, HC, CM))
cor.ens

##############################################ENSEMBLE - KMÉDIAS##########################
#Utilizando k-médias para fazer um unico agrupamento dos 4 modelos

#Elbow método, para identificar o número ideal de clusters para o dataset
k.max <- 15
wss <- sapply(1:k.max, 
              function(k){kmeans(data.ens[,-12], k, nstart=20,iter.max = 15 )$tot.withinss})
wss
windows()
plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


library(NbClust)
nb_data_ens <- NbClust(data.ens[,-12], distance = "euclidean", min.nc = 2,
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

#***** Conclusion *****                            
  
#  * According to the majority rule, the best number of clusters is  3 


#******************************************************************* 
  

windows()
hist(nb_data_ens$Best.nc[1,], breaks = max(na.omit(nb_data_ens$Best.nc[1,])))


library(factoextra)
windows()
fviz_nbclust(data.ens[,-12], FUN = hcut, method = "silhouette")+geom_vline(xintercept = 5, linetype = 2)
windows()
fviz_nbclust(data.ens[,-12], FUN = kmeans, method = "silhouette")+geom_vline(xintercept = 5, linetype = 2)



#kmeans
library(stats)
set.seed(5555)



#Criar modelo
ens.model.kmeans <- kmeans(data.ens[,-12], 3, nstart = 20)#5 centróides, o r irá tentar  20 começos diferentes

km <- kmeans(official_trajectory_matrix, 2, nstart = 20)#5 centróides, o r irá tentar  20 começos diferentes
windows()
fviz_cluster(km, geom = "point",  data = official_trajectory_matrix) + ggtitle("k = 2")



library(factoextra)
pp3 <- fviz_cluster(ens.model.kmeans, geom = "point",  data = data.ens) + ggtitle("k = 3")
windows()
plot(pp3)

library(cluster) 
windows()
clusplot(data.ens[,-12], ens.model.kmeans$cluster, color=TRUE, shade=TRUE, 
         labels=1, lines=0)
# Centroid Plot against 1st 2 discriminant functions
#library(fpc)
#windows()
#plotcluster(data.ens[,-12], ens.model.kmeans$cluster)



#Plotar gráficos
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

data.final <- as.data.frame (ens.model.kmeans$cluster)
data.final["PC1"] <- c(official_trajectory_matrix[,1])
data.final["PC2"] <- c(official_trajectory_matrix[,2])
data.final["PC3"] <- c(official_trajectory_matrix[,3])
data.final["PC4"] <- c(official_trajectory_matrix[,4])
data.final["PC5"] <- c(official_trajectory_matrix[,5])
data.final["PC6"] <- c(official_trajectory_matrix[,6])
data.final["PC7"] <- c(official_trajectory_matrix[,7])
data.final["PC8"] <- c(official_trajectory_matrix[,8])
data.final["PC9"] <- c(official_trajectory_matrix[,9])
data.final["PC10"] <- c(official_trajectory_matrix[,10])

#data.final["Date"]= c(Date)
#colnames(data.final) = c("Clusters", "Data")
#data.final = data.final[,c(2,1)]

###################################### Validação ###############################
library("clValid")
# comparing hierarchical and kmeans with internal validation indexes
clmethods <- c("hierarchical","kmeans")
intern <- clValid(data.ens[,-12], nClust = 2:13, maxitems = 2364, 
                  clMethods = clmethods, validation = "internal", method = "complete")

summary(intern)
#_______________________________________________

library("factoextra")

# can change the clustering algorithm in FUN (ex kmeans)
windows()
fviz_nbclust(data.ens[,-12], FUN = hcut, method = "wss")+geom_vline(xintercept = 5, linetype = 2)
windows()
fviz_nbclust(data.ens[,-12], FUN = hcut, method = "silhouette", ylim =0.8)+geom_vline(xintercept = 5, linetype = 2)

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

#clus4 <- as.matrix(official_trajectory_matrix)[which(ens.model.kmeans$cluster ==4),]
#mean(clus4)
#sd(clus4)

#clus5 <- as.matrix(official_trajectory_matrix)[which(ens.model.kmeans$cluster ==5),]
#mean(clus5)
#sd(clus5)

#clus6 <- as.matrix(official_trajectory_matrix)[which(ens.model.kmeans$cluster ==6),]
#mean(clus6)
#sd(clus6)

#clus7 <- as.matrix(official_trajectory_matrix)[which(ens.model.kmeans$cluster ==7),]
#mean(clus7)
#sd(clus7)

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
#clus4<-diag.ave(clus4)
#clus5<-diag.ave(clus5)
#clus6<-diag.ave(clus6)
#clus7<-diag.ave(clus7)
traj_matrix <- diag.ave(as.matrix(official_trajectory_matrix))
windows()
plot(traj_matrix, type = "l")