library(copent)
library(mnormt)
library(energy)
library(Ball)
library(copula)
library(hypoRF)
library(HHG)
library(cramer)
library(kernlab)
library(TwoSampleTest.HD)
library(fasano.franceschini.test)
library(Peacock.test)
library(diproperm)
library(latex2exp)

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

m0 = c(0,0)
rho = 0.5 # simulation 1
rho = 0.0 # simulation 2,3
v0 = matrix(c(1,rho,rho,1), nrow = 2)
n0 = 400 # size of sample0
n1 = 450 # size of sample1
sample0 = rmnorm(n0,m0,v0)
stat1 = kstat1 = stat2 = estat1 = ball1 = rf1 = hhg1 = hhg2 = hhg3 = hhg4 = cramer1 = tsthd1 = ff1 = peacock1 = dpp1 = rep(0,10)
for(i in 1:10){
  # simulation 1
  # m1 = m0 + i - 1
  # sample1 = rmnorm(n1,m1,v0)
  
  # simulation 2
  rho1 = (i-1) * 0.1
  v1 = matrix(c(1,rho1,rho1,1),nrow = 2)
  sample1 = rmnorm(n1,m0,v1)
  
  # simulation 3
  # mv.cop <- mvdc(normalCopula(0.8), c("norm", "exp"), list(list(mean = 0, sd = 2), list(rate = i)))
  # sample1 <- rMvdc(n1, mv.cop)

  stat1[i] = tst(sample0,sample1)
  stat2[i] = tst.mi(sample0,sample1)
  kstat1[i] = kmmd(sample0,sample1)@mmdstats[1]
  estat1[i] = eqdist.e( rbind(sample0,sample1), c(n0,n1) ) # energy statistics based test
  ball1[i] = bd.test(sample0,sample1)$statistic # Ball divergence based test
  rf1[i] = hypoRF(as.data.frame(sample0),as.data.frame(sample1))$obs
  Dx = as.matrix(dist((rbind(sample0,sample1)), diag = TRUE, upper = TRUE))
  y = c(rep(0,n0),rep(1,n1))
  hhg = hhg.test.2.sample(Dx,y, nr.perm = 1000)
  hhg1[i] = hhg$sum.chisq; hhg2[i] = hhg$sum.lr; hhg3[i] = hhg$max.chisq; hhg4[i] = hhg$max.lr
  cramer1[i] = cramer.test(sample0,sample1)$statistic
  tsthd1[i] = TwoSampleTest.HD(sample0,sample1)$statistic
  ff1[i] = fasano.franceschini.test(sample0,sample1)$statistic
  peacock1[i] = peacock2(sample0,sample1)
  dpp1[i] = DiProPerm(rbind(sample0,sample1),c(rep(-1,n0),rep(1,n1)))$obs_teststat
}#i

x11(width = 12, height = 6)
par(mfrow = c(3,5))
# x1 = seq(0,9); xlab1 = TeX(r'($u_1$)') # simulation 1
x1 = seq(0,0.9,0.1); xlab1 = TeX(r'($\rho_1$)') # simulation 2
# x1 = 1:10; xlab1 = "rate" # simulation 3
plot(x1,stat1, xlab = xlab1, ylab = "statistic", main = "CE");lines(x1,stat1)
plot(x1,stat2, xlab = xlab1, ylab = "statistic", main = "MI");lines(x1,stat2)
plot(x1,kstat1, xlab = xlab1, ylab = "statistic", main = "Kernel");lines(x1,kstat1)
plot(x1,estat1, xlab = xlab1, ylab = "statistic", main = "Energy");lines(x1,estat1)
plot(x1,ball1, xlab = xlab1, ylab = "statistic", main = "Ball");lines(x1,ball1)
plot(x1,rf1, xlab = xlab1, ylab = "statistic", main = "RF");lines(x1,rf1)
plot(x1,hhg1, xlab = xlab1, ylab = "statistic", main = "HHG sum.chisq");lines(x1,hhg1)
plot(x1,hhg2, xlab = xlab1, ylab = "statistic", main = "HHG sum.lr");lines(x1,hhg2)
plot(x1,hhg3, xlab = xlab1, ylab = "statistic", main = "HHG max.chisq");lines(x1,hhg3)
plot(x1,hhg4, xlab = xlab1, ylab = "statistic", main = "HHG max.lr");lines(x1,hhg4)
plot(x1,cramer1, xlab = xlab1, ylab = "statistic", main = "Cramer");lines(x1,cramer1)
plot(x1,tsthd1, xlab = xlab1, ylab = "statistic", main = "TST.HD");lines(x1,tsthd1)
plot(x1,ff1, xlab = xlab1, ylab = "statistic", main = "F-F");lines(x1,ff1)
plot(x1,peacock1, xlab = xlab1, ylab = "statistic", main = "Peacock");lines(x1,peacock1)
plot(x1,dpp1, xlab = xlab1, ylab = "statistic", main = "DiProPerm");lines(x1,dpp1)
