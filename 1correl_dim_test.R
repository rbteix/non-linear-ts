#Importing data set
setwd("E://Documentos//Mestrado//Projeto")
data_gru <- read.csv("dados_gr.csv", sep = ";")


#Hides a non numeric data
data_gru1 <- na.omit(data_gru)
na.fail(data_gru1) #verifies if it has only numeric observations



#Removes outliers
data_gru2 <- data_gru1
data_gru2<-data_gru2[!(data_gru2$"Realizado">=750 | 
                         data_gru2$"Realizado"<=450 |
                         data_gru2$"Realizado.com.atraso">=350 |
                         data_gru2$"Cancelado" > 205), ]

# Column date, converting factor in date
Date <-as.Date(data_gru2$Data.1, format = "%d/%m/%Y ")


#Planned flights(including all flights:canceled, realized and delayed)
Planejados <- data_gru2$Total.geral 
Planejados = data_gru2$Cancelado + data_gru2$Realizado + data_gru2$Realizado.com.atraso

#Fail percentage (canceled + realized)/ total
data_fail<- (data_gru2$Cancelado+data_gru2$Realizado.com.atraso) / data_gru2$Total.geral

#Plot
windows()
plot(data_fail, type = "l")

#Converts in data frame
data_fail<- matrix(data_fail, length(data_fail), 1)
colnames(data_fail)= c("Fail_Percentage") #Gives column name
data_fail = as.data.frame(data_fail)

#Tests of linearity
################################### BDS TEST ##################################

library(fNonlinear)
##Warning messages:
##1: package 'fNonlinear' was built under R version 3.4.4 
##2: package 'timeSeries' was built under R version 3.4.4 
##3: package 'fBasics' was built under R version 3.4.4
bds <-bdsTest(unlist(data_fail), m = 6, eps = NULL, title = NULL, description = NULL) 

##################################Surrogate test###############
library(nonlinearTseries)
##Warning messages:
##1: package 'nonlinearTseries' was built under R version 3.4.4 
##2: package 'rgl' was built under R version 3.4.4 
##3: package 'TSA' was built under R version 3.4.4 
##4: package 'leaps' was built under R version 3.4.4 
##5: package 'locfit' was built under R version 3.4.4 

windows()
surrogate <- surrogateTest(time.series = unlist(data_fail), significance = 0.05, one.sided= FALSE,
             K=1, FUN=timeAsymmetry, verbose = T, do.plot = T, xlab = "Values of the statistic",
             ylab = "", main = "Surrogate data testing")

###############################Embedding nonlinear package#####################
library(nonlinearTseries)

# Estimates tau-delay (time step) based on Average Mutual Information
tm.ami <- timeLag(unlist(data_fail), technique = "ami", 
                   lag.max = 50, do.plot = T)


#Uses Cao's algorithm to estimate an apropriate embedding dimension  
windows()
emb.dm = estimateEmbeddingDim(unlist(data_fail), time.lag = tm.ami,
                               max.embedding.dim = 15)

#Dimension correlation
library(nonlinearTseries)

data_fail_num<- as.numeric(unlist(data_fail)) # converts it in numeric data

#Calculates the correlation dimension
windows()
cdi = corrDim(data_fail_num,
             min.embedding.dim = emb.dm,
             max.embedding.dim = emb.dm + 5,
             time.lag = tm.ami, 
             min.radius = 0.001, max.radius = 70, #Verificar se há regiões lineares em diferentes dimensões embeddig
             n.points.radius = 40,
             theiler.window = 4,
             do.plot=FALSE)
plot(cdi)

cd.est = estimate(cdi, regression.range=c(0.2,50),
                  use.embeddings = 15:17)

cat("expected: 0  --- estimate: ",cd.est,"\n")
