library(copent)
library(mnormt)
library(maotai)
library(energy)
library(Ball)
library(copula)
library(hypoRF)
library(HHG)
library(cramer)

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
rho = 0.5 # simulation 1
# rho = 0.0 # simulation 2,3
v0 = matrix(c(1,rho,rho,1), nrow = 2)
n0 = 300 # size of sample0
n1 = 350 # size of sample1
sample0 = rmnorm(n0,m0,v0)
stat1 = kstat1 = stat2 = estat1 = ball1 = rf1 = hhg1 = hhg2 = hhg3 = hhg4 = cramer1 = rep(0,10)
for(i in 1:10){
  # simulation 1
  m1 = m0 + i - 1
  sample1 = rmnorm(n1,m1,v0)
  
  # simulation 2
  # rho1 = (i-1) * 0.1
  # v1 = matrix(c(1,rho1,rho1,1),nrow = 2)
  # sample1 = rmnorm(n1,m0,v1)
  
  # simulation 3
  # mv.cop <- mvdc(normalCopula((i-1)*0.1), c("norm", "exp"), list(list(mean = 0, sd =2), list(rate = 0.5)))
  # sample1 <- rMvdc(n1, mv.cop)

  stat1[i] = tst(sample0,sample1)
  stat2[i] = tst.mi(sample0,sample1)
  kstat1[i] = tst.kernel(sample0,sample1)
  estat1[i] = eqdist.e( rbind(sample0,sample1), c(n0,n1) ) # energy statistics based test
  ball1[i] = bd.test(sample0,sample1)$statistic # Ball divergence based test
  rf1[i] = hypoRF(as.data.frame(sample0),as.data.frame(sample1))$obs
  Dx = as.matrix(dist((rbind(sample0,sample1)), diag = TRUE, upper = TRUE))
  y = c(rep(0,n0),rep(1,n1))
  hhg = hhg.test.2.sample(Dx,y, nr.perm = 1000)
  hhg1[i] = hhg$sum.chisq; hhg2[i] = hhg$sum.lr; hhg3[i] = hhg$max.chisq; hhg4[i] = hhg$max.lr
  cramer1[i] = cramer.test(sample0,sample1)$statistic
}#i

x11(width = 12, height = 9)
par(mfrow = c(3,4))
plot(stat1, ylab = "statistic", main = "CE");lines(stat1)
plot(stat2, ylab = "statistic", main = "MI");lines(stat2)
plot(kstat1, ylab = "statistic", main = "Kernel");lines(kstat1)
plot(estat1, ylab = "statistic", main = "Energy");lines(estat1)
plot(ball1, ylab = "statistic", main = "Ball");lines(ball1)
plot(rf1, ylab = "statistic", main = "RF");lines(rf1)
plot(hhg1, ylab = "statistic", main = "HHG sum.chisq");lines(hhg1)
plot(cramer1, ylab = "statistic", main = "Cramer");lines(cramer1)
plot(hhg2, ylab = "statistic", main = "HHG sum.lr");lines(hhg2)
plot(hhg3, ylab = "statistic", main = "HHG max.chisq");lines(hhg3)
plot(hhg4, ylab = "statistic", main = "HHG max.lr");lines(hhg4)
