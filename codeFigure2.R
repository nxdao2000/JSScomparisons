##
##  R code to generate Figure 2 and 3 of main text to
##  "Simulation-based inference methods for partially observed Markov model via the R package is2"
##  by Duc Anh Doan, Dao Nguyen and Xin Dang.
##

library("is2")
library("coda")
library("foreach")
library("doMC")

##' The 'bake' function caches results in an .rda file.
bake <- function(file, expr) {
  expr <- substitute(expr)
  if (file.exists(file)) {
    load(file, envir = parent.frame())
  } else {
    env <- new.env()
    eval(expr, envir = env)
    save(list = ls(envir = env), envir = env, file = file)
    for (n in ls(envir = env)) assign(n, get(n, envir = env), envir = parent.frame())
  }
  invisible(NULL)
}
pompExample(gompertz)

bigtick <- Sys.time()

##' Edit the following line to obtain better results using package 'doMC'.
options(cores = 8)

##' ## Constructing a pomp object.
##' The following codes construct the basic elements of the Gompertz model
##' and construct the 'gompertz' pomp object.
#' ## Particle MCMC

##' We'll need a prior density function:
hyperparams <- list(min = coef(gompertz)/10, max = coef(gompertz) * 10)

gompertz.dprior <- function(params, ..., log) {
  f <- sum(dunif(params, min = hyperparams$min, max = hyperparams$max, log = TRUE))
  if (log) 
    f else exp(f)
}
bake("pmcmc.rda", {
  
  pompExample(gompertz)
  
  tic <- Sys.time()
  library("doMC")
  library("foreach")
  registerDoMC()
  
  save.seed <- .Random.seed
  set.seed(12345678L, kind = "L'Ecuyer")
  
  pmcmc1 <- foreach(i = 1:4, .inorder = FALSE, .packages = "pomp", .combine = c, 
                    .options.multicore = list(set.seed = TRUE)) %dopar% {
                      pmcmc(pomp(gompertz, dprior = gompertz.dprior), start = coef(gompertz), Nmcmc = 10000, 
                            Np = 100, max.fail = Inf, proposal = mvn.diag.rw(c(r = 0.01, sigma = 0.01, 
                                                                               tau = 0.01)))
                    }
  
  toc <- Sys.time()
  pmcmcTime <- toc - tic
  .Random.seed <<- save.seed
  
  pmcmc.traces <- conv.rec(pmcmc1, c("sigma", "tau"))
  pmcmc.traces <- window(pmcmc.traces, start = 5000, thin = 1)
  ess.pmcmc <- effectiveSize(pmcmc.traces)
  rm(pmcmc1, save.seed, tic, toc)
})
ess.pmcmc

##' Plot the traces and densities.
pdf('pmcmc-gompertz.pdf')
plot(pmcmc.traces)
dev.off()


bake("pmif.rda", {
  
  pompExample(gompertz)
  
  tic <- Sys.time()
  library("doMC")
  library("foreach")
  registerDoMC()
  
  save.seed <- .Random.seed
  set.seed(12345678L, kind = "L'Ecuyer")
  
  pmif1 <- foreach(i = 1:4, .inorder = FALSE, .packages = "is2", .combine = c, 
                   .options.multicore = list(set.seed = TRUE)) %dopar% {
                     pmif(pomp(gompertz, dprior = gompertz.dprior), start = coef(gompertz), 
                          Nmif = 10000, 
                          Np = 100, max.fail = Inf, proposal = mvn.diag.rw(c(r = 0.01, sigma = 0.01, 
                                                                             tau = 0.01)))
                   }
  
  toc <- Sys.time()
  pmifTime <- toc - tic
  .Random.seed <<- save.seed
  
  pmif.traces <- conv.rec(pmif1, c("sigma", "tau"))
  pmif.traces <- window(pmif.traces,start=5000, thin = 1)
  ess.pmif <- effectiveSize(pmif.traces)
  rm(pmif1, save.seed, tic, toc)
})
ess.pmif
##' Plot the traces and densities.
pdf('pmif-gompertz.pdf')
plot(pmif.traces)
dev.off()
  
  
  








