logit <- function(x){log10(x/(1-x))}
tigol <- function(x){1/(1+10^(-x))}


N = 1000
x1 = tigol(rnorm(N,-2, 0.5)) 
x2 = tigol(rnorm(N, 2, 0.5))
x3 = tigol(rnorm(N, 0, 0.2))
x4 = tigol(rnorm(N, 0, 0.2))

df1 = data.frame(mean=x1, var=x1*(1-x1)*runif(N,0,0.8), col=1)
df2 = data.frame(mean=x2, var=x2*(1-x2)*runif(N,0,0.8), col=2)
df3 = data.frame(mean=x3, var=x3*(1-x3)*tigol(rnorm(N, -0.8,  0.5)),col=3)
df4 = data.frame(mean=x4, var=x4*(1-x4)*tigol(rnorm(N, 0.8,  0.5)),col=4)


Q = rbind(df1,df2,df3,df4)
colnames(Q) = c('mean','var','col')

b = 50

layout(matrix(c(0,5,0,2,1,3,0,4,0),3,3,byrow=T), widths=c(1,2,1),heights=c(1,2,1))

x=seq(0,1,0.02)
plot(x,x*(1-x),'l',lty=2,xlab='mean',ylab='variance') + with(Q,points(mean, var, col=col, pch=19,cex=0.1))

fr <- function(x) {c(x,0,1)}

hist(fr(tigol(rnorm(N,-2, 0.5))), breaks=b, xlab=expression(Psi), ylab='',main='',xlim=c(0,1), freq=T, yaxt='n',border=1,col=1)
hist(fr(tigol(rnorm(N, 2, 0.5))), breaks=b, xlab=expression(Psi), ylab='',main='',xlim=c(0,1), freq=T, yaxt='n',border=2,col=2)

hist(fr(tigol(rnorm(N, 0,  0.1))), breaks=b, xlab=expression(Psi), ylab='',main='',xlim=c(0,1), freq=T, yaxt='n', border=3,col=3)
u = rbinom(N,1,0.6)
hist(fr((1-u)*tigol(rnorm(N,-2, 0.5))+u*tigol(rnorm(N,2, 0.5))), breaks=b, xlab=expression(Psi), ylab='',main='',xlim=c(0,1), freq=T, yaxt='n', border=4,col=4)



