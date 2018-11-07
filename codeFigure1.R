##
##  R code to generate Figure 1 of main text to
##  "Simulation-based inference methods for partially observed Markov model via the R package is2"
##  by Duc Anh Doan, Dao Nguyen and Xin Dang.
##

CORES <- 8      ## number of parallel processes
JOBS <- 30      ## number of estimation runs
NMIF <- 20      ## number of IF iterations per estimation run
NP <- 1000      ## number of particles in pfilter operations
seed <- 123
require(doParallel)
registerDoParallel(CORES)
require(is2)

pompExample(ou2)
p.truth <- coef(ou2)
guess1 <- p.truth
th2.range <- c(-1,0)
th3.range <- c(0,1)
X1.range <- c(-5,5)
X2.range <- c(-5,5)
tic <- Sys.time()

binary.file <- "ou2toy.rda"
if (file.exists(binary.file)){
  load(binary.file)
} else {
  ## RUN MIF1 TO ESTIMATE PARAMETERS
  tic <- Sys.time()
  mpar1 <- foreach(i=1:JOBS, 
                   .packages='pomp', 
                   .inorder=FALSE) %dopar% {
    
    jobno<-seed+i
    set.seed(jobno)
    th.draw <-c(alpha.2=runif(1,min=th2.range[1],max=th2.range[2]), 
                alpha.3=runif(1,min=th3.range[1],max=th3.range[2]),
                X1.0=runif(1,min=X1.range[1],max=X1.range[2]),
                X2.0=runif(1,min=X2.range[1],max=X2.range[2]))
    guess1[c('alpha.2','alpha.3','x1.0','x2.0')] <-th.draw
    
    oo<-try(mif(ou2,Nmif=NMIF,start=guess1,
                pars=c('alpha.2','alpha.3'),ivps=c('x1.0','x2.0'),
                rw.sd=c(x1.0=0.2,x2.0=0.2,
                        alpha.2=0.02,alpha.3=0.02),
                Np=NP,
                var.factor=1,
                ic.lag=10,
                method="mif",
                cooling.type="geometric",
                cooling.fraction=0.3,
                max.fail=100
    )
    )
    if(class(oo)[1]=='try-error'){
      c(alpha.2=NA,alpha.3=NA)
    } 
    else{
      coef(oo)[c('alpha.2','alpha.3')]
    }
                   }
  toc <- Sys.time()
  etime1 <- toc-tic
  m1.out <- do.call(rbind,mpar1)
  m1.out<-m1.out[is.finite(m1.out[,1]),]
  m1.out<-t(m1.out)
  m1.out
  etime1
  ## RUN MIF2 TO ESTIMATE PARAMETERS
  tic <- Sys.time()
  mpar2 <- foreach(i=1:JOBS,
                   .packages='is2',
                   .inorder=FALSE) %dopar% {
    
  jobno<-seed+i
  set.seed(jobno)
  th.draw <-c(alpha.2=runif(1,min=th2.range[1],max=th2.range[2]), 
              alpha.3=runif(1,min=th3.range[1],max=th3.range[2]),
              X1.0=runif(1,min=X1.range[1],max=X1.range[2]),
              X2.0=runif(1,min=X2.range[1],max=X2.range[2]))
  guess1[c('alpha.2','alpha.3','x1.0','x2.0')] <-th.draw
  oo<-try(mif3(ou2,Nmif=NMIF,start=guess1,
               ivps=c('x1.0','x2.0'),
               rw.sd=c(x1.0=0.2,x2.0=0.2,
                       alpha.2=0.02,alpha.3=0.02),
               Np=NP,
               var.factor=1,
               ic.lag=10,
               method="mif2",
               cooling.type="geometric",
               cooling.fraction=0.3,
               max.fail=100
  )
  )
  if(class(oo)[1]=='try-error'){
    c(alpha.2=NA,alpha.3=NA)
  } 
  else {
    coef(oo)[c('alpha.2','alpha.3')]
  }
                   }
  
  toc <- Sys.time()
  etime2 <- toc-tic

  m2.out <- do.call(rbind,mpar2)
  m2.out<-m2.out[is.finite(m2.out[,1]),]
  m2.out<-t(m2.out)
  m2.out
  etime2
  ## RUN IS2 lag 1 TO ESTIMATE PARAMETERS
  
  tic <- Sys.time()
  mpar3 <- foreach(i=1:JOBS,
                   .packages='is2',
                   .inorder=FALSE) %dopar% {
    jobno<-seed+i
    set.seed(jobno)
    th.draw <-c(alpha.2=runif(1,min=th2.range[1],max=th2.range[2]), 
                alpha.3=runif(1,min=th3.range[1],max=th3.range[2]),
                X1.0=runif(1,min=X1.range[1],max=X1.range[2]),
                X2.0=runif(1,min=X2.range[1],max=X2.range[2]))
    guess1[c('alpha.2','alpha.3','x1.0','x2.0')] <-th.draw
  
    oo<-try(is2(ou2,Nis=NMIF,start=guess1,
                pars=c('alpha.2','alpha.3'),ivps=c('x1.0','x2.0'),
                rw.sd=c(x1.0=0.2,x2.0=0.2,
                        alpha.2=0.02,alpha.3=0.02),
                Np=NP,
                var.factor=1,
                ic.lag=10,
                method="is2",
                lag=1,
                cooling.type="geometric",
                cooling.fraction=0.3,
                max.fail=100
              )
            )
    if(class(oo)[1]=='try-error'){
      c(alpha.2=NA,alpha.3=NA)
    } 
    else{
      coef(oo)[c('alpha.2','alpha.3')]
    }
  }
  toc <- Sys.time()
  etime3 <- toc-tic

  m3.out <- do.call(rbind,mpar3)
  m3.out<-m3.out[is.finite(m3.out[,1]),]
  m3.out<-t(m3.out)
  m3.out
  etime3
  ## RUN IS2 lag 5 TO ESTIMATE PARAMETERS

  tic <- Sys.time()
  mpar4 <- foreach(i=1:JOBS,
                   .packages='is2',
                   .inorder=FALSE) %dopar% {
    jobno<-seed+i
    set.seed(jobno)
    th.draw <-c(alpha.2=runif(1,min=th2.range[1],max=th2.range[2]), 
                alpha.3=runif(1,min=th3.range[1],max=th3.range[2]),
                X1.0=runif(1,min=X1.range[1],max=X1.range[2]),
                X2.0=runif(1,min=X2.range[1],max=X2.range[2]))
    guess1[c('alpha.2','alpha.3','x1.0','x2.0')] <-th.draw
  
    oo<-try(is2(ou2,Nis=NMIF,start=guess1,
                pars=c('alpha.2','alpha.3'),ivps=c('x1.0','x2.0'),
                rw.sd=c(x1.0=2,x2.0=2,
                        alpha.2=.23,alpha.3=.23),
                Np=NP,
                var.factor=1,
                ic.lag=10,
                method="is2",
                lag=5,
                cooling.type="geometric",
                cooling.fraction=0.3,
                max.fail=100
              )
            )
    if(class(oo)[1]=='try-error'){
      c(alpha.2=NA,alpha.3=NA)
    } 
    else{
      coef(oo)[c('alpha.2','alpha.3')]
    }
  }
  toc <- Sys.time()
  etime4 <- toc-tic

  m4.out <- do.call(rbind,mpar4)
  m4.out<-m4.out[is.finite(m4.out[,1]),]
  m4.out<-t(m4.out)
  m4.out
  etime4
  
  ## RUN avif lag 1 TO ESTIMATE PARAMETERS
  
  tic <- Sys.time()
  mpar5 <- foreach(i=1:JOBS, 
                   .packages='is2', 
                   .inorder=FALSE) %dopar% {
                     jobno<-seed+i
                     set.seed(jobno)
                     th.draw <-c(alpha.2=runif(1,min=th2.range[1],max=th2.range[2]), 
                                 alpha.3=runif(1,min=th3.range[1],max=th3.range[2]),
                                 X1.0=runif(1,min=X1.range[1],max=X1.range[2]),
                                 X2.0=runif(1,min=X2.range[1],max=X2.range[2]))
                     guess1[c('alpha.2','alpha.3','x1.0','x2.0')] <-th.draw
                     oo<-try(avif(ou2,Nis=NMIF,start=guess1,
                                  ivps=c('x1.0','x2.0'),
                                  rw.sd=c(x1.0=0.2,x2.0=0.2,
                                          alpha.2=0.02,alpha.3=0.02),
                                  Np=NP,
                                  var.factor=1,
                                  ic.lag=10,
                                  method="avif",
                                  lag=1,
                                  cooling.type="geometric",
                                  cooling.fraction=0.3,
                                  max.fail=100
                     )
                     )
                     if(class(oo)[1]=='try-error'){
                       c(alpha.2=NA,alpha.3=NA)
                     } 
                     else {
                       coef(oo)[c('alpha.2','alpha.3')]
                     }
                   }
  toc <- Sys.time()
  etime5 <- toc-tic
  
  m5.out <- do.call(rbind,mpar5)
  m5.out<-m5.out[is.finite(m5.out[,1]),]
  m5.out<-t(m5.out)
  m5.out
  etime5
  
  ## RUN aif filtering TO ESTIMATE PARAMETERS
  
  tic <- Sys.time()
  mpar6 <- foreach(i=1:JOBS, 
                   .packages='is2', 
                   .inorder=FALSE) %dopar% {
                     jobno<-seed+i
                     set.seed(jobno)
                     th.draw <-c(alpha.2=runif(1,min=th2.range[1],max=th2.range[2]), 
                                 alpha.3=runif(1,min=th3.range[1],max=th3.range[2]),
                                 X1.0=runif(1,min=X1.range[1],max=X1.range[2]),
                                 X2.0=runif(1,min=X2.range[1],max=X2.range[2]))
                     guess1[c('alpha.2','alpha.3','x1.0','x2.0')] <-th.draw
                     oo<-try(aif(ou2,Nis=NMIF,start=guess1,
                                  ivps=c('x1.0','x2.0'),
                                  rw.sd=c(x1.0=0.2,x2.0=0.2,
                                          alpha.2=0.02,alpha.3=0.02),
                                  Np=NP,
                                  var.factor=1,
                                  ic.lag=10,
                                  method="aif",
                                  lag=0,
                                  cooling.type="geometric",
                                  cooling.fraction=0.3,
                                  max.fail=100
                     )
                     )
                     if(class(oo)[1]=='try-error'){
                       c(alpha.2=NA,alpha.3=NA)
                     } 
                     else {
                       coef(oo)[c('alpha.2','alpha.3')]
                     }
                   }
  toc <- Sys.time()
  etime6 <- toc-tic
  
  m6.out <- do.call(rbind,mpar6)
  m6.out<-m6.out[is.finite(m6.out[,1]),]
  m6.out<-t(m6.out)
  m6.out
  etime6
  
  ## RUN aif lag 5 TO ESTIMATE PARAMETERS
  
  tic <- Sys.time()
  mpar7 <- foreach(i=1:JOBS, 
                   .packages='is2', 
                   .inorder=FALSE) %dopar% {
                     jobno<-seed+i
                     set.seed(jobno)
                     th.draw <-c(alpha.2=runif(1,min=th2.range[1],max=th2.range[2]), 
                                 alpha.3=runif(1,min=th3.range[1],max=th3.range[2]),
                                 X1.0=runif(1,min=X1.range[1],max=X1.range[2]),
                                 X2.0=runif(1,min=X2.range[1],max=X2.range[2]))
                     guess1[c('alpha.2','alpha.3','x1.0','x2.0')] <-th.draw
                     oo<-try(aif(ou2,Nis=NMIF,start=guess1,
                                 ivps=c('x1.0','x2.0'),
                                 rw.sd=c(x1.0=0.2,x2.0=0.2,
                                         alpha.2=0.02,alpha.3=0.02),
                                 Np=NP,
                                 var.factor=1,
                                 ic.lag=10,
                                 method="aif",
                                 lag=5,
                                 cooling.type="geometric",
                                 cooling.fraction=0.3,
                                 max.fail=100
                     )
                     )
                     if(class(oo)[1]=='try-error'){
                       c(alpha.2=NA,alpha.3=NA)
                     } 
                     else {
                       coef(oo)[c('alpha.2','alpha.3')]
                     }
                   }
  toc <- Sys.time()
  etime7 <- toc-tic
  
  m7.out <- do.call(rbind,mpar7)
  m7.out<-m7.out[is.finite(m7.out[,1]),]
  m7.out<-t(m7.out)
  m7.out
  etime7
  
  ## RUN mif_Momentum TO ESTIMATE PARAMETERS
  
  mpar8 <- foreach(i=1:JOBS,
                   .packages='is2',
                   .inorder=FALSE) %dopar% {
                     jobno<-seed+i
                     set.seed(jobno)
                     th.draw <-c(alpha.2=runif(1,min=th2.range[1],max=th2.range[2]), 
                                 alpha.3=runif(1,min=th3.range[1],max=th3.range[2]),
                                 X1.0=runif(1,min=X1.range[1],max=X1.range[2]),
                                 X2.0=runif(1,min=X2.range[1],max=X2.range[2]))
                     guess1[c('alpha.2','alpha.3','x1.0','x2.0')] <-th.draw
                     
                     oo<-try(mifMomentum(ou2,Nis=NMIF,start=guess1,
                                          pars=c('alpha.2','alpha.3'),ivps=c('x1.0','x2.0'),
                                          rw.sd=c(x1.0=0.2,x2.0=0.2,
                                                  alpha.2=0.02,alpha.3=0.02),
                                          Np=NP,
                                          var.factor=1,
                                          ic.lag=10,
                                          cooling.type="geometric",
                                          cooling.fraction=0.3,
                                          max.fail=100
                     )
                     )
                     if(class(oo)[1]=='try-error'){
                       c(alpha.2=NA,alpha.3=NA)
                     } 
                     else{
                       coef(oo)[c('alpha.2','alpha.3')]
                     }
                   }
  toc <- Sys.time()
  etime8 <- toc-tic
  
  m8.out <- do.call(rbind,mpar8)
  m8.out<-m8.out[is.finite(m8.out[,1]),]
  m8.out<-t(m8.out)
  m8.out
  etime8
  
  
  set.seed(577639485L)

  kalman.filter <- function (y, x0, a, b, sigma, tau) {
    n <- nrow(y)
    ntimes <- ncol(y)
    sigma.sq <- sigma%*%t(sigma)
    tau.sq <- tau%*%t(tau)
    inv.tau.sq <- solve(tau.sq)
    cond.dev <- numeric(ntimes)
    filter.mean <- matrix(0,n,ntimes)
    pred.mean <- matrix(0,n,ntimes)
    pred.var <- array(0,dim=c(n,n,ntimes))
    dev <- 0
    m <- x0
    v <- diag(0,n)
    for (k in seq(length=ntimes)) {
      pred.mean[,k] <- M <- a%*%m
      pred.var[,,k] <- V <- a%*%v%*%t(a)+sigma.sq
      q <- b%*%V%*%t(b)+tau.sq
      r <- y[,k]-b%*%M
      cond.dev[k] <- n*log(2*pi)+log(det(q))+t(r)%*%solve(q,r)
      dev <- dev+cond.dev[k]
      q <- t(b)%*%inv.tau.sq%*%b+solve(V)
      v <- solve(q)
      filter.mean[,k] <- m <- v%*%(t(b)%*%inv.tau.sq%*%y[,k]+solve(V,M))
    }
    list(
      pred.mean=pred.mean,
      pred.var=pred.var,
      filter.mean=filter.mean,
      cond.loglik=-0.5*cond.dev,
      loglik=-0.5*dev
    )
  }

  kalman <- function (x, object, params) {
    y <- object@data
    p <- params
    p[names(x)] <- x
    x0 <- init.state(object,params=p)
    a <- matrix(p[c('alpha.1','alpha.2','alpha.3','alpha.4')],2,2)
    b <- diag(1,2)
    sigma <- matrix(p[c('sigma.1','sigma.2','sigma.2','sigma.3')],2,2)
    sigma[1,2] <- 0
    tau <- diag(p['tau'],2,2)
    -kalman.filter(y,x0,a,b,sigma,tau)$loglik
  }

  # true coefficients
  p.truth <- coef(ou2)
  cat("coefficients at `truth'\n")
  print(p.truth[c('alpha.2','alpha.3','x1.0','x2.0')],digits=4)
  cat("Kalman filter log likelihood at truth\n")
  print(loglik.truth <- -kalman(p.truth,ou2,p.truth),digits=4)
  
  # make a wild guess
  p.guess <- p.truth[c('alpha.2','alpha.3','x1.0','x2.0')]*exp(rnorm(n=4,mean=0,sd=0.5))
  cat("coefficients at guess\n")
  print(p.guess,digits=4)
  cat("Kalman filter log likelihood at guess\n")
  print(loglik.guess <- -kalman(p.guess,ou2,p.truth),digits=4)
  
  # find MLE using Kalman filter starting at the guess
  cat("running Kalman filter estimation\n")
  tic <- Sys.time()
  kalm.fit1 <- optim(p.guess,kalman,object=ou2,params=p.truth,hessian=T,control=list(trace=2))
  toc <- Sys.time()
  print(toc-tic)
  tic <- Sys.time()
  print(loglik.mle <- -kalm.fit1$value,digits=4)
  mle.exact<-kalm.fit1$par
  
  ### COMPUTE EXACT LIKELIHOODS
  
  N1 <- 100
  N2 <- 100
  th1.vec <- seq(from=-1,to=0,length=N1)
  th2.vec <- seq(from=0,to=1,length=N2)
  lik.array <- matrix(NA,nr=N1,nc=N2)
  
  for (i1 in 1:N2) {
   for (i2 in 1:N2) {
     th <- c(alpha.2=th1.vec[i1],alpha.3=th2.vec[i2])
     lik.array[i1,i2]<--kalman(th, object=ou2, params=p.truth)
   }
  }
  save(m1.out,etime1,
       m2.out,etime2,
       m3.out,etime3,
       m4.out,etime4,
       m5.out,etime5,
       m6.out,etime6,
       m7.out,etime7,
       m8.out,etime8,
       mle.exact,lik.array,
       th1.vec,th2.vec,
       file="ou2toy.rda",compress='xz')
}
load("ou2toy.rda")

pdf("ou2compareAllSupp3.pdf")
XLIM <- c(-0.95,-0.05)
YLIM <- c(0.05,0.95)
LINE.YAXIS <- 2
LINE.XAXIS <- 2.5
LEVELS <- c(0,2,4,10,10000)
X.LABEL <- 0.87
Y.LABEL <- 0.87	

op <- par(mfrow=c(2,4),mai=c(0,0.2,0.2,0),omi=c(0.7,0.4,0,0.1))
CEX.TRIANGLE <- CEX.SQUARE <- 1.5
CEX.POINTS <- 1
CEX.LAB <- 1.5
CEX.AXIS <- 1.2

CEX.TRIANGLE <- CEX.SQUARE <- 1
CEX.POINTS <- 1
CEX.LAB <- 1.2
CEX.AXIS <- 0.8
CEX.AXIS.NUMBERS <- 1
#col1=heat.colors(12)
col1=terrain.colors(20,alpha = 1)
########### if1 ###############

require(graphics)
#require(colorspace)
#require(RColorBrewer)
plot(x=m1.out["alpha.2",],y=m1.out["alpha.3",],xlim=XLIM,ylim=YLIM,
     xlab='',ylab='',axes=F,type='n')
box()
axis(side=2,cex.axis=CEX.AXIS.NUMBERS)
axis(side=1,labels=F,cex.axis=CEX.AXIS.NUMBERS)
mtext(side=2,bquote(alpha[2]),line=LINE.YAXIS,cex=CEX.AXIS)
.filled.contour(x=th1.vec,y=th2.vec,z=(max(lik.array)-lik.array),levels=LEVELS,col=col1[c(1,5,8,11)])
points(x=m1.out["alpha.2",],y=m1.out["alpha.3",],pch=3,cex=CEX.POINTS)
points(x=mle.exact["alpha.2"],y=mle.exact["alpha.3"],pch=3,col="green",cex=CEX.TRIANGLE)

plot.window(c(0,1),c(0,1))
text(x=X.LABEL,y=Y.LABEL,"A",cex=CEX.LAB)

########### if2 ###############

plot(x=m2.out["alpha.2",],y=m2.out["alpha.3",],xlim=XLIM,ylim=YLIM,
     xlab='',ylab='',axes=F,type='n')
box()
axis(side=2,labels=F,cex.axis=CEX.AXIS.NUMBERS)
axis(side=1,labels=F,cex.axis=CEX.AXIS.NUMBERS)

.filled.contour(x=th1.vec,y=th2.vec,z=(max(lik.array)-lik.array),levels=LEVELS,col=col1[c(1,5,8,11)])

points(x=m2.out["alpha.2",],y=m2.out["alpha.3",],pch=3,cex=CEX.POINTS)
points(x=mle.exact["alpha.2"],y=mle.exact["alpha.3"],pch=3,col="green",cex=CEX.TRIANGLE)

plot.window(c(0,1),c(0,1))
text(x=X.LABEL,y=Y.LABEL,"B",cex=CEX.LAB)

########### is2 lag 1 ###############

plot(x=m3.out["alpha.2",],y=m3.out["alpha.3",],xlim=XLIM,ylim=YLIM,
     xlab='',ylab='',axes=F,type='n')
box()
axis(side=2,labels=F,cex.axis=CEX.AXIS.NUMBERS)
axis(side=1,labels=F,cex.axis=CEX.AXIS.NUMBERS)

.filled.contour(x=th1.vec,y=th2.vec,z=(max(lik.array)-lik.array),levels=LEVELS,col=col1[c(1,5,8,11)])

points(x=m3.out["alpha.2",],y=m3.out["alpha.3",],pch=3,cex=CEX.POINTS)
points(x=mle.exact["alpha.2"],y=mle.exact["alpha.3"],pch=3,col="green",cex=CEX.TRIANGLE)

plot.window(c(0,1),c(0,1))
text(x=X.LABEL,y=Y.LABEL,"C",cex=CEX.LAB)

########### is2 lag 5 ###############

plot(x=m4.out["alpha.2",],y=m4.out["alpha.3",],xlim=XLIM,ylim=YLIM,
     xlab='',ylab='',axes=F,type='n')
box()
axis(side=2,labels=F,cex.axis=CEX.AXIS.NUMBERS)
axis(side=1,labels=F,cex.axis=CEX.AXIS.NUMBERS)

.filled.contour(x=th1.vec,y=th2.vec,z=(max(lik.array)-lik.array),levels=LEVELS,col=col1[c(1,5,8,11)])

points(x=m4.out["alpha.2",],y=m4.out["alpha.3",],pch=3,cex=CEX.POINTS)
points(x=mle.exact["alpha.2"],y=mle.exact["alpha.3"],pch=3,col="green",cex=CEX.TRIANGLE)

plot.window(c(0,1),c(0,1))
text(x=X.LABEL,y=Y.LABEL,"D",cex=CEX.LAB)


########### avif ###############

plot(x=m5.out["alpha.2",],y=m5.out["alpha.3",],xlim=XLIM,ylim=YLIM,
     xlab='',ylab='',axes=F,type='n')
box()
axis(side=2,cex.axis=CEX.AXIS.NUMBERS)
axis(side=1,cex.axis=CEX.AXIS.NUMBERS)
mtext(side=2,bquote(alpha[2]),line=LINE.YAXIS,cex=CEX.AXIS)
mtext(side=1,bquote(alpha[3]),line=LINE.YAXIS,cex=CEX.AXIS)

.filled.contour(x=th1.vec,y=th2.vec,z=(max(lik.array)-lik.array),levels=LEVELS,col=col1[c(1,5,8,11)])

points(x=m5.out["alpha.2",],y=m5.out["alpha.3",],pch=3,cex=CEX.POINTS)
#points(x=m2.out["th1.start",],y=m2.out["th2.start",],pch=3,cex=0.3)
points(x=mle.exact["alpha.2"],y=mle.exact["alpha.3"],pch=3,col="green",cex=CEX.TRIANGLE)

plot.window(c(0,1),c(0,1))
text(x=X.LABEL,y=Y.LABEL,"E",cex=CEX.LAB)

########### aif ###############

plot(x=m6.out["alpha.2",],y=m6.out["alpha.3",],xlim=XLIM,ylim=YLIM,
     xlab='',ylab='',axes=F,type='n')
box()
axis(side=2,labels=F,cex.axis=CEX.AXIS.NUMBERS)
axis(side=1,cex.axis=CEX.AXIS.NUMBERS)
mtext(side=1,bquote(alpha[3]),line=LINE.YAXIS,cex=CEX.AXIS)

.filled.contour(x=th1.vec,y=th2.vec,z=(max(lik.array)-lik.array),levels=LEVELS,col=col1[c(1,5,8,11)])

points(x=m6.out["alpha.2",],y=m6.out["alpha.3",],pch=3,cex=CEX.POINTS)
points(x=mle.exact["alpha.2"],y=mle.exact["alpha.3"],pch=3,col="green",cex=CEX.TRIANGLE)

plot.window(c(0,1),c(0,1))
text(x=X.LABEL,y=Y.LABEL,"F",cex=CEX.LAB)


########### aif lag 5 ###############

plot(x=m7.out["alpha.2",],y=m7.out["alpha.3",],xlim=XLIM,ylim=YLIM,
     xlab='',ylab='',axes=F,type='n')
box()
axis(side=2,labels=F,cex.axis=CEX.AXIS.NUMBERS)
axis(side=1,cex.axis=CEX.AXIS.NUMBERS)
mtext(side=1,bquote(alpha[3]),line=LINE.YAXIS,cex=CEX.AXIS)

.filled.contour(x=th1.vec,y=th2.vec,z=(max(lik.array)-lik.array),levels=LEVELS,col=col1[c(1,5,8,11)])

points(x=m7.out["alpha.2",],y=m7.out["alpha.3",],pch=3,cex=CEX.POINTS)
points(x=mle.exact["alpha.2"],y=mle.exact["alpha.3"],pch=3,col="green",cex=CEX.TRIANGLE)

plot.window(c(0,1),c(0,1))
text(x=X.LABEL,y=Y.LABEL,"G",cex=CEX.LAB)

########### mif.Momentum ###############

plot(x=m8.out["alpha.2",],y=m8.out["alpha.3",],xlim=XLIM,ylim=YLIM,
     xlab='',ylab='',axes=F,type='n')
box()
axis(side=2,labels=F,cex.axis=CEX.AXIS.NUMBERS)
axis(side=1,cex.axis=CEX.AXIS.NUMBERS)
mtext(side=1,bquote(alpha[3]),line=LINE.YAXIS,cex=CEX.AXIS)

.filled.contour(x=th1.vec,y=th2.vec,z=(max(lik.array)-lik.array),levels=LEVELS,col=col1[c(1,5,8,11)])

points(x=m8.out["alpha.2",],y=m8.out["alpha.3",],pch=3,cex=CEX.POINTS)
points(x=mle.exact["alpha.2"],y=mle.exact["alpha.3"],pch=3,col="green",cex=CEX.TRIANGLE)

plot.window(c(0,1),c(0,1))
text(x=X.LABEL,y=Y.LABEL,"H",cex=CEX.LAB)


par(op)
dev.off()


