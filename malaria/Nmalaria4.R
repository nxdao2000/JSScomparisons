##
##  R code to generate Figure 3 of main text to
##  "Simulation-based inference methods for partially observed Markov model via the R package is2"
##  by Duc Anh Doan, Dao Nguyen and Xin Dang.
##

require(is2)

CORES <- 8     ## number of parallel processes
JOBS <- 8      ## number of estimation runs
NMIF <- 2      ## number of IF iterations per estimation run
NP <- 1000     ## number of particles in pfilter operations
NLIK <- 2      ## number of likelihood evaluations
NPLIK <- 2000  ## number of particles in likelihood evaluation
## for mif computations:
IC.LAG <- 100   ## fixed smoothing lag for initial conditions

 require(doParallel)
registerDoParallel(CORES)

#patch relapse.c
if(!is.loaded("relapse")){
  system('R CMD SHLIB relapse.c')
  dyn.load("relapse.so")
}

# define scaling functions
logit <- function(p){log(p/(1-p))}		# maps interval [0,1] to [-infty,infty]
expit <- function(p){1/(1+exp(-p))}		# reverse map of logit()

# import data
dat <- read.csv("City-data.csv")							# case + rainfall data
census  <- read.csv("City-census.csv",comment.char="#")	# population data
time.dat <- dat$time[which(dat$time>1987)]					# t0 -> Jan 1987
pv.dat <- dat$pv[which(dat$time>1987)]
pop.year <- census[which(census$Year>=1987),1]
pop.dat <- census[which(census$Year>=1987),2]
rain.dat <- dat$rainfall[which(dat$time>1986)]				# rain data begins a year before (Jan 1986)
# replace NA values in rainfall data by monthly averages
nn <- length(rain.dat)
for(i in 1:12){
  ix <- seq(i,nn,12)
  avr <- mean(rain.dat[ix],na.rm=TRUE)
  rain.dat[ix][is.na(rain.dat[ix])] <- avr
}
# accumulate rainfall (adding preceding 4 months of rain)
rain <- rain.dat[9:296]										# use data from Sep 1986
nn <- length(pv.dat)
acrain.dat <- rep(0,nn)
for(i in 1:nn) acrain.dat[i] <- sum(rain[i:(i+3)])
# remove Jan-Jun & Nov-Dec rain (for max correlation with Pv)
for(i in seq(1,nn,12)) acrain.dat[c(i:(i+5),(i+10):(i+11))] <- 0
acrain.dat <- acrain.dat[1:nn]		# reset to original length

#compute classical basis as seasonality
tbasis.step <- 1/360
tbasis <- seq(min(time.dat)-0.1,max(time.dat)+0.1,by=tbasis.step)
nbasis <- 6
basis <- periodic.bspline.basis(tbasis,nbasis=nbasis,degree=3,period=1)
colnames(basis) <- paste("season",1:nbasis,sep="")

# smoothen (decadal) population data
pop.spline <- smooth.spline(x=pop.year,y=pop.dat)
pop <- predict(pop.spline,tbasis)$y
dpopdt <- c(0,diff(pop))/tbasis.step

# smoothen rainfall data by linear interpolation (time points left-shifted one month)
acrain <- approx(x=time.dat-1/12,y=acrain.dat,xout=tbasis,rule=2)$y
# define covariates
covariate <- cbind(pop=pop,dpopdt=dpopdt,rainfall=(acrain-mean(acrain))/sd(acrain),basis=basis)
######################################
# build pomp object for SEIH^3QS model
######################################
po <- pomp(
  data = data.frame(cbind(time=time.dat,pv=pv.dat)),
  times = "time",
  t0 = 2*time.dat[1] - time.dat[2],
  tcovar = tbasis,
  covar = covariate,
  statenames = c("cases","err","S","E","I","Q","H1","H2","H3","K","F","W"),
  paramnames = c("log.sigOBS","log.muEI","log.muIH","log.muHI","log.muIQ","log.muIS","log.muQS",
    			 "log.tau","log.sigPRO","logit.rho","logit.q","logit.br","delta","b1","logit.S0",
				 "logit.E0","logit.I0","logit.Q0","logit.H10","logit.H20","logit.H30","log.K0","log.F0"),
  covarnames = c("pop","dpopdt","rainfall","season1"),
  zeronames = c("cases","W","err"),
  obsnames = c("pv"),
  rprocess = euler.sim(step.fun="SEIH3QS",delta.t=2/365),
  rmeasure = "nbinom_rmeasure",
  dmeasure = "nbinom_dmeasure",
  initializer = function(params,t0,statenames,covars,...)
  {
    # untransform needed parameters
    snames <-c("cases","err","S","E","I","Q","H1","H2","H3","K","F","W")
    frac1 <- expit(params[c("logit.S0","logit.E0","logit.I0","logit.Q0","logit.H10","logit.H20","logit.H30")])
    frac2 <- exp(params[c("log.K0","log.F0")])
    x0 <- numeric(length(snames))
    names(x0) <- snames
    x0[c("S","E","I","Q","H1","H2","H3","K","F")] <- c(round(covars["pop"]*frac1/sum(frac1)),frac2)
    x0
  }
)



#########################
# run IF1 
#########################


METHOD <- "avif"
IC.LAG <- 1000
tic <- Sys.time()
mpar4 <- foreach(
  i=1:JOBS,
  .packages='pomp',
  .inorder=FALSE
) %dopar% {
  
  set.seed(345+i)
  


# randomize (and transform where needed) initial parameter values
# ("b1...b6" are untransformed, "delta" and "initpop" not estimated)
par.ini.tr <- c(
  log.muEI=log(24), log.muIH=log(runif(1,1,5)), log.muHI=log(runif(1,1,5)),
  log.muIQ=log(runif(1,1,2)), log.muIS=log(runif(1,0.5,2)), log.muQS=log(runif(1,10,20)),
  log.sigPRO=log(runif(1,0.1,0.5)), log.sigOBS=log(runif(1,0.1,0.5)), log.tau=log(runif(1,0.1,0.5)),
  log.K0=log(runif(1,0.1,0.5)), log.F0=log(runif(1,0.1,0.5)),
  logit.rho=logit(runif(1,0.001,0.01)), logit.q=logit(0.001), logit.br=logit(runif(1,0.5,0.9)),
  logit.S0=0, logit.E0=0, logit.I0=0, logit.Q0=0, logit.H10=0, logit.H20=0, logit.H30=0,
  b1=runif(1,-5,5), b2=runif(1,-5,5), b3=runif(1,-5,5), b4=runif(1,-5,5),
  b5=runif(1,-5,5), b6=runif(1,-5,5), delta=0.02, initpop=1163607
)
# randomize S0,E0,I0,Q0,H0 with constraint S0+E0+I0+Q0+H0=1
S0  <- runif(1,0,1);
E0  <- runif(1,0,1)
I0  <- runif(1,0,1)
Q0  <- runif(1,0,1)
H10 <- runif(1,0,1)
H20 <- runif(1,0,1)
H30 <- runif(1,0,1)
fac <- sum(c(S0,E0,I0,Q0,H10,H20,H30))
S0 <- S0/fac; E0 <- E0/fac; I0 <- I0/fac; Q0 <- Q0/fac; H10 <- H10/fac; H20 <- H20/fac; H30 <- H30/fac
par.ini.tr[c("logit.S0","logit.E0","logit.I0","logit.Q0","logit.H10","logit.H20","logit.H30")] <- logit(c(S0,E0,I0,Q0,H10,H20,H30))

# log-transformed parameters

log.par <- c(
  "log.muEI", "log.muIH", "log.muHI", "log.muIQ", "log.muIS", "log.muQS", "log.tau",
  "log.sigPRO", "log.sigOBS", "log.K0", "log.F0"
)
# logit-transformed parameters
logit.par <- c(
  "logit.rho", "logit.q", "logit.br", "logit.S0", "logit.E0", "logit.I0", "logit.Q0",
  "logit.H10", "logit.H20", "logit.H30"
)
# untransformed parameters
untr.par <- c("b1", "b2", "b3", "b4", "b5", "b6", "delta", "initpop")
# untransformed initial parameters
par.ini <- c(exp(par.ini.tr[log.par]),expit(par.ini.tr[logit.par]),par.ini.tr[untr.par])
names(par.ini) <- c(
  "muEI", "muIH", "muHI", "muIQ", "muIS", "muQS", "tau", "sigPRO", "sigOBS", "K0", "F0", "rho",
  "q", "br", "S0", "E0", "I0", "Q0", "H10", "H20", "H30", "b1", "b2", "b3", "b4", "b5", "b6",
  "delta", "initpop"
)

# parameters to be estimated
par.est <- c(
  "log.muEI", "log.muIH", "log.muHI", "log.muIQ", "log.muIS", "log.muQS", "log.sigPRO", "log.sigOBS",
  "log.tau", "logit.rho", "logit.q", "logit.br", "b1", "b2", "b3", "b4", "b5", "b6"
)
# initial states to be estimated
ivp.est <- c("logit.S0", "logit.E0", "logit.I0", "logit.Q0", "logit.H10", "logit.H20", "logit.H30", "log.K0", "log.F0")

# strength of random walk
sd <- rep(0,length(par.ini.tr))
names(sd) <- names(par.ini.tr)
sd[par.est] <- 0.1
sd[ivp.est] <- 2

# parameters to be held fixed at pre-assigned value
# vivax incubation period ~= 15 days (=15/365 yr); q takes a low value <<1 
par.ini.tr[c("log.muEI","logit.q")] <- c(log(24),logit(0.001))
sd[c("log.muEI","logit.q")] <- 0.000001
par.tr <- par.ini.tr

  m <- try(
    is3(
    po,
    Nmif = NMIF,
    start = par.tr,
    pars = par.est,
    ivps = ivp.est,
    rw.sd = sd,
    Np = NP,
    var.factor = 1,
    ic.lag = 10,
    cooling.type="geometric",
    cooling.fraction=0.1,		
    max.fail = 1000,
    method="avif",
    lag=0
    
    
  )
  )
  if(class(m)[1]=='try-error'){
    list(pomp=po,start=par.tr,ll=NA)
  } else{
    ll <- logLik(pfilter(m,Np=NPLIK))
    ll1 <-m$loglik
    list(pomp=m,start=par.tr,ll=ll, ll1=ll1)
  }
}
toc <- Sys.time()
etime4 <- toc-tic


m4.in <- rbind(sapply(mpar4,function(x)x$start))
m4.out <- rbind(sapply(mpar4,function(x)coef(x$pomp)))
m4.lik <- rbind(sapply(mpar4,function(x)x$ll))

m4.lik

save(m4.in,m4.out,m4.lik,etime4,
		 file="Nmalaria-is3-avif-2-1000.rda",compress='xz')

