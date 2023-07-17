library(copent)
library(mnormt)

tst<-function(s0,s1){
  n0 = dim(s0)[1]
  n1 = dim(s1)[1]
  x = rbind(s0,s1)
  y1 = c(rep(1,n0),rep(2,n1)) + runif(n0+n1, max = 0.0000001)
  y0 = c(rep(1,n0+n1)) + runif(n0+n1,max = 0.0000001)
  
  copent(cbind(x,y1)) - copent(cbind(x,y0))
}

m0 = c(0,0)
rho = 0.5
v0 = matrix(c(1,rho,rho,1), nrow = 2)
n = 500 # sample size
sample0 = rmnorm(n,m0,v0)
data = matrix(0,n*11,2)
data[1:n,] = sample0
stat1 = rep(0,10)
for(i in 1:10){
  m1 = m0 + i - 1
  sample1 = rmnorm(n,m1,v0)
  data[1:n+i*n,] = sample1

  stat1[i] = tst(sample0,sample1)
}

x11()
plot(stat1, ylab = "estimated statistic");lines(stat1)

x11()
col1 = rep(1:11,each = n)
plot(data, xlab = "x1", ylab = "x2", col = col1)
legend(x = min(data[,1]),y = max(data[,2]), legend = 0:10,col = 1:11, pch = 1)
