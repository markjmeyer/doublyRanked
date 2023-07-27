#### load libraries ####
library(mvnfast)
library(refund)
library(numbers)
library(xtable)
source('wilcox.fda.R')

dt_check <- try(library(depthTests), silent = TRUE)
if(class(dt_check) == "try-error"){
  require(devtools)
  devtools::install_github("julia-wrobel/depthTests")
}
library(depthTests)

#### specifications ####
n1    <- n2 <- 10
N     <- n1 + n2
S     <- c(40, 80, 160, 320)
sc    <- 1
B     <- 1000
sig   <- 0.05
rho   <- c(0.25, 0.5, 0.75)
dots  <- 1000
up    <- 25000

pvals   <- array(NA, dim = c(length(rho), length(S), B, 2))
stats   <- array(NA, dim = c(length(rho), length(S), B, 2))

##### group mean functions and covariance ####
iter <- 0

for(i in 1:length(rho)){
  for(k in 1:length(S)){
    # i <- k <- 1
    CovStr  <- matrix(0, nrow = S[k], ncol = S[k])
    for(r in 1:S[k]){
      for(c in 1:S[k]){
        CovStr[r,c] <- sig*rho[i]^(abs(r - c))
      }
    }
    for(b in 1:B){
      set.seed(b)
      s     <- seq(0, 1, length = S[k])
      mu1   <- dgamma(s, 10, 30)/max(dgamma(s, 10, 30))
      mu2   <- sc*mu1
      X1    <- rmvn(n = n1, mu = mu1, sigma = CovStr)
      X2    <- rmvn(n = n2, mu = mu2, sigma = CovStr)
      X     <- rbind(X1, X2)
      fX    <- fpca.face(X)
      Y     <- fX$Yhat
      G     <- c(rep(0, n1), rep(1, n2))
      
      ###### doubly ranked test ######
      drt   <- wilcox.fda(Y ~ G, method = 'geo.rank')
      pvals[i, k, b, 1] <- drt$p.value
      stats[i, k, b, 1] <- drt$statistic
      
      ###### MBD rank sum test ######
      mrs   <- rankSumTest(X1, X2, floor(n1/2), max(floor(n2/2)-1, 2))
      pvals[i, k, b, 2] <- mrs$p.value
      stats[i, k, b, 2] <- mrs$statistic

      
      ## simulation controls ##
      iter <- iter + 1
      
      if(mod(iter, dots) == 0){
        cat('.')
      }
      
      if(mod(iter, up) == 0){
        cat(paste("\n",iter,"datasets completed\n"))
      }           
    }
  }
}

#### save output ####
filename <- paste('mww_gamma_g_size_', n1, '_fpca.RData', sep = '')
save.image(file = filename)