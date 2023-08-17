library(copent)

s0 = matrix( rnorm(300,0,1), 300, 1)
stat1 = 0
for(i in 1:10){
  s1 = matrix( rnorm(350,(i-1)*0.5,1), 350, 1)
  stat1[i] = tst(s0,s1)
}
x11(); plot(stat1, ylab = "statistic", main = "univariate two-sample test"); lines(stat1)

