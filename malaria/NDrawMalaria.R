##
##  R code to generate Figure 3 of main text to
##  "Simulation-based inference methods for partially observed Markov model via the R package is2"
##  by Duc Anh Doan, Dao Nguyen and Xin Dang.
##

load("Nmalaria-mif-50-1000.rda")
load("Nmalaria-mif2-50-1000.rda")
load("Nmalaria-is2-50-1000.rda")
load("Nmalaria-is3-aif-50-1000.rda")

m1.lik<-m1.lik[is.finite(m1.lik)]
m2.lik<-m2.lik[is.finite(m2.lik)]
m3.lik<-m3.lik[is.finite(m3.lik)]
m7.lik<-m7.lik[is.finite(m7.lik)]


res1<-m1.lik
res2<-m2.lik
res3<-m3.lik
res7<-m7.lik

MLE<-max(max(res1),max(res2),max(res3),max(res7))


pdf("compareAllmalaria.pdf")

s1<-log(MLE-res1)

k1<- density(s1, bw=1.5)
L1<- MLE-exp(k1$x)
L1
plot(L1,k1$y/exp(k1$x), main="",type="l",lty=3, col=3, lwd=2, axes=T,xlab='',ylab='', xlim=c(-2000,-1800),ylim=c(0,0.06))


s2<-log(MLE-res2)
k2<- density(s2, bw=1.5)
L2<- MLE-exp(k2$x)
lines(L2,k2$y/exp(k2$x), main="",lty=2, col=2, lwd=2,axes=T,xlab='',ylab='', ylim=c(0,0.06), xlim=c(-2000,-1800))


s3<-log(MLE-res3)
k3<- density(s3, bw=1.5)
L3<- MLE-exp(k3$x)
lines(L3,k3$y/exp(k3$x), main="",lty=1, col=1, lwd=2,axes=T,xlab='',ylab='', ylim=c(0,0.06), xlim=c(-2000,-1800))

s7<-log(MLE-res7)
k7<- density(s7, bw=1.5)
L7<- MLE-exp(k7$x)
lines(L7,k7$y/exp(k7$x), main="",lty=5, col=5, lwd=2,axes=T,xlab='',ylab='', ylim=c(0,0.06), xlim=c(-2000,-1800))



abline(v = MLE, col = "darkblue", lwd=2, lty=5)

XLIM <- c(-2000,-1800)
YLIM <- c(0,0.06)
LINE.YAXIS <- 2
LINE.XAXIS <- 2.5
X.LABEL <- 0.87
Y.LABEL <- 0.87  

CEX.TRIANGLE <- CEX.SQUARE <- 1.5
CEX.POINTS <- 1.5
CEX.LAB <- 0.9
CEX.AXIS <- 0.9

CEX.TRIANGLE <- CEX.SQUARE <- 1
CEX.POINTS <- 1
CEX.LAB <- 1.2
CEX.AXIS <- 0.9
CEX.AXIS.NUMBERS <- 1


box()
axis(side=1,cex=CEX.AXIS.NUMBERS)
mtext(side=1,bquote("log likelihood"),line=LINE.XAXIS,cex=CEX.AXIS)

abline(h=0)
legend(-2000,0.06, c("IF1","IF2", "IS2", "AIF"), col = c(3,2,1,5),
        lty = c(3, 2, 1,5),
       merge = TRUE, bg = "gray90")
plot.window(c(0,1),c(0,1))
dev.off()
