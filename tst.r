library(copent)
library(mnormt)
library(maotai)
library(energy)
library(copula)

# copula entropy based test
tst.ce<-function(s0,s1,n = 8){
  n0 = dim(s0)[1]
  n1 = dim(s1)[1]
  x = rbind(s0,s1)
  result = 0
  for(i in 1:n){
    y1 = c(rep(1,n0),rep(2,n1)) + runif(n0+n1, max = 0.0000001)
    y0 = c(rep(1,n0+n1)) + runif(n0+n1,max = 0.0000001)
    result = result + copent(cbind(x,y1)) - copent(cbind(x,y0))
  }
  result/n
}

# mutual information based test
tst.mi<-function(s0,s1,n = 8){
  n0 = dim(s0)[1]
  n1 = dim(s1)[1]
  x = rbind(s0,s1)
  result = 0
  for(i in 1:n){
    y1 = c(rep(1,n0),rep(2,n1)) + runif(n0+n1, max = 0.0000001)
    result = result + copent(cbind(x,y1)) - copent(x)
  }
  result / n
}

# kernel-based test
tst.kernel<-function(s0,s1){
  dmat <- as.matrix(dist(rbind(s0, s1))) 
  kmat <- exp(-(dmat^2))                     
  lab  <- c(rep(1,dim(s0)[1]), rep(2,dim(s1)[1])) 
  
  mmd2test(kmat, lab, method = "u")$statistic
}

m0 = c(0,0)
# rho = 0.5 # simulation 1
rho = 0.0 # simulation 2,3
v0 = matrix(c(1,rho,rho,1), nrow = 2)
n = 500 # sample size
sample0 = rmnorm(n,m0,v0)
# data = matrix(0,n*11,2)
# data[1:n,] = sample0
stat1 = kstat1 = stat2 = estat1 = rep(0,10)
for(i in 1:10){
  # simulation 1
  # m1 = m0 + i - 1
  # sample1 = rmnorm(n,m1,v0)

  # simulation 2
  # rho1 = (i-1) * 0.1
  # v1 = matrix(c(1,rho1,rho1,1),nrow = 2)
  # sample1 = rmnorm(n,m0,v1)
  
  # simulation 3
  mv.cop <- mvdc(normalCopula((i-1)*0.1), c("norm", "exp"), list(list(mean = 0, sd =2), list(rate = 0.5)))
  sample1 <- rMvdc(n, mv.cop)

  # data[1:n+i*n,] = sample1
  
  stat1[i] = tst.ce(sample0,sample1)
  stat2[i] = tst.mi(sample0,sample1)
  kstat1[i] = tst.kernel(sample0,sample1)
  estat1[i] = eqdist.e( rbind(sample0,sample1), c(n,n) ) # energy statistics based test
}#i

x11()
par(mfrow = c(2,2))
plot(stat1, ylab = "estimated statistic", main = "CE");lines(stat1)
plot(stat2, ylab = "estimated statistic", main = "MI");lines(stat2)
plot(kstat1, ylab = "estimated statistic", main = "Kernel");lines(kstat1)
plot(estat1, ylab = "estimated statistic", main = "Energy");lines(estat1)


# x11()
# col1 = rep(1:11,each = n)
# plot(data, xlab = "x1", ylab = "x2", col = col1)
# legend(x = min(data[,1]),y = max(data[,2]), legend = 0:10,col = 1:11, pch = 1)
