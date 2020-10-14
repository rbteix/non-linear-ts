library(tseriesChaos)

data_lyp <- data_fail
data_lyp = as.ts(data_lyp)

S_nu <- lyap_k(data_lyp,m=6,d=5,t=4,k=8,ref=1500,s=200, eps = 0.05)

windows()
plot(S_nu,xlab = expression(paste(nu)),ylab=expression(paste("S",(nu))))
lmin<- 50
lmax<- 200
lyap(S_nu, 50, 200)
# add the vertical lines delimiting the linear region:
#abline(v=lmin,lty=4,lwd=2,col="black")
#abline(v=lmax,lty=4,lwd=2,col="black")

#computing the slope of the segment for v to linear line, through regression analysis
step <- 0.01
lr<- S_nu[lmin:lmax] # output of lyap_k() in [lmin, lmax]
l_lr<- length(lr)-1
lt<- seq(lmin*step,lmax*step,by=step)
lm(lr~lt) # regression analysis in R

#abline(lyap(S_nu,lmin,lmax),lty=2,lwd=2,col="black") # add the regression line
#abline()

#Parâmetros
#séries : séries temporais
#m: dimensão embedding
#d: tempo de atraso
#t: Janela de Theiler - pontos como menos que na séries são
#exluídos da busca por pontos vizinhos.
#k: número de vizinhos considerados
#ref: número de pontos considerados da série
#s: iterações
#eps: raio dentro do qual os vizinhos mais próximos são procurados
#definição do que é uma trajetória próxima.

