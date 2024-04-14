library(copent)
library(twosamples)

sample0 = matrix( rnorm(300,0,1), 300, 1)
stat1 = cvm1 = ks1 = kuiper1 = wass1 = dts1 = ad1 = wilcox1 = kruskal1 = 0
for(i in 1:10){
  sample1 = matrix( rnorm(350,(i-1)*0.5,1), 350, 1)
  stat1[i] = tst(sample0,sample1)
  wilcox1[i] = wilcox.test(sample0,sample1)$statistic
  kruskal1[i] = kruskal.test(c(sample0,sample1),c(rep(1,300),rep(2,350)))$statistic
  cvm1[i] = cvm_stat(sample0,sample1)
  ks1[i] = ks_stat(sample0,sample1)
  kuiper1[i] = kuiper_stat(sample0,sample1)
  wass1[i] = wass_stat(sample0,sample1)
  dts1[i] = dts_stat(sample0,sample1)
  ad1[i] = ad_stat(sample0,sample1)
}
x11(); 
par(mfrow=c(3,3))
plot(stat1, ylab = "statistic", main = "CE"); lines(stat1)
plot(wilcox1, ylab = "statistic", main = "Wilcoxon"); lines(wilcox1)
plot(kruskal1, ylab = "statistic", main = "Kruskal-Wallis"); lines(kruskal1)
plot(cvm1, ylab = "statistic", main = "CVM");lines(cvm1)
plot(ks1, ylab = "statistic", main = "KS");lines(ks1)
plot(kuiper1, ylab = "statistic", main = "Kuiper");lines(kuiper1)
plot(wass1, ylab = "statistic", main = "WASS");lines(wass1)
plot(dts1, ylab = "statistic", main = "DTS");lines(dts1)
plot(ad1, ylab = "statistic", main = "AD");lines(ad1)
