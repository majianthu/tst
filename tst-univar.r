library(copent)
library(twosamples)
library(robusTest)
library(R2sample)
library(tnl.Test)
library(latex2exp)

sample0 = matrix( rnorm(300,0,1), 300, 1)
stat1 = cvm1 = ks1 = kuiper1 = wass1 = dts1 = ad1 = wilcox1 = kruskal1 = wilcox2 = vartest1 = 0
lr1 = za1 = zk1 = zc1 = tnl1 = 0
for(i in 1:10){
  # sample1 = matrix( rnorm(350,(i-1),1), 350, 1) # mean
  sample1 = matrix( rnorm(350, 0, i), 350, 1) # var
  stat1[i] = tst(sample0,sample1)
  wilcox1[i] = wilcox.test(sample0,sample1)$statistic
  kruskal1[i] = kruskal.test(c(sample0,sample1),c(rep(1,300),rep(2,350)))$statistic
  cvm1[i] = cvm_stat(sample0,sample1)
  ks1[i] = ks_stat(sample0,sample1)
  kuiper1[i] = kuiper_stat(sample0,sample1)
  wass1[i] = wass_stat(sample0,sample1)
  dts1[i] = dts_stat(sample0,sample1)
  ad1[i] = ad_stat(sample0,sample1)
  # wilcox2[i] = wilcoxtest(sample0,sample1)$statistic # simulation 1
  vartest1[i] = vartest(sample0,sample1)$statistic # simulation 2
  lr1[i] = twosample_test(sample0,sample1,doMethod = "LR")$statistics
  za1[i] = twosample_test(sample0,sample1,doMethod = "ZA")$statistics
  zk1[i] = twosample_test(sample0,sample1,doMethod = "ZK")$statistics
  zc1[i] = twosample_test(sample0,sample1,doMethod = "ZC")$statistics
  tnl1[i] = tnl.test(sample0,sample1,2)$statistic
}
x11(width = 12, height = 8);
par(mfrow=c(3,5))
# x1 = seq(0,9); xlab1 = TeX(r'($u_1$)') # simulation 1
x1 = 1:10; xlab1 = TeX(r'($\delta_1$)') # simulation 2
plot(x1,stat1, xlab = xlab1, ylab = "statistic", main = "CE"); lines(x1,stat1)
plot(x1,wilcox1, xlab = xlab1, ylab = "statistic", main = "Wilcoxon1"); lines(x1,wilcox1)
plot(x1,kruskal1, xlab = xlab1, ylab = "statistic", main = "Kruskal-Wallis"); lines(x1,kruskal1)
plot(x1,cvm1, xlab = xlab1, ylab = "statistic", main = "CVM");lines(x1,cvm1)
plot(x1,ks1, xlab = xlab1, ylab = "statistic", main = "KS");lines(x1,ks1)
plot(x1,kuiper1, xlab = xlab1, ylab = "statistic", main = "Kuiper");lines(x1,kuiper1)
plot(x1,wass1, xlab = xlab1, ylab = "statistic", main = "WASS");lines(x1,wass1)
plot(x1,dts1, xlab = xlab1, ylab = "statistic", main = "DTS");lines(x1,dts1)
plot(x1,ad1, xlab = xlab1, ylab = "statistic", main = "AD");lines(x1,ad1)
# plot(x1,wilcox2, xlab = xlab1, ylab = "statistic", main = "Wilcoxon2"); lines(x1,wilcox2) # simulation 1
plot(x1,vartest1, xlab = xlab1, ylab = "statistic", main = "Vartest"); lines(x1,vartest1) # simulation 2
plot(x1,lr1, xlab = xlab1, ylab = "statistic", main = "LR");lines(x1,lr1)
plot(x1,za1, xlab = xlab1, ylab = "statistic", main = "ZA");lines(x1,za1)
plot(x1,zk1, xlab = xlab1, ylab = "statistic", main = "ZK");lines(x1,zk1)
plot(x1,zc1, xlab = xlab1, ylab = "statistic", main = "ZC");lines(x1,zc1)
plot(x1,tnl1, xlab = xlab1, ylab = "statistic", main = "TNL");lines(x1,tnl1)
