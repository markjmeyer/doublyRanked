#### load libraries ####
library(refund)
library(numbers)
library(xtable)
library(compositions)
source('wilcox.fda.R')

#### specifications ####
n1    <- n2   <-  n3 <- 10
N     <- n1 + n2 + n3
one   <- TRUE # TRUE, FALSE
S     <- c(40, 80, 160, 320)
sc    <- seq(0.99, 0.5, by = -0.02)
B     <- 500
sig   <- 0.05
rho   <- c(0.25, 0.5, 0.75)
dots  <- 1000
up    <- 25000

pvals   <- array(NA, dim = c(length(rho), length(S), length(sc), B))
stats   <- array(NA, dim = c(length(rho), length(S), length(sc), B))

##### group mean functions and covariance ####
iter <- 0

for(i in 1:length(rho)){
  for(k in 1:length(S)){
    CovStr  <- matrix(0, nrow = S[k], ncol = S[k])
    for(r in 1:S[k]){
      for(c in 1:S[k]){
        CovStr[r,c] <- sig*rho[i]^(abs(r - c))
      }
    }
    for(l in 1:length(sc)){
      for(b in 1:B){
        set.seed(b)
        s     <- seq(0, 1, length = S[k])
        mu1   <- pnorm(s, mean = 0.5, sd = 0.1)
        mu2   <- sc[l]*mu1
        mu3   <- ifelse(one, 1, sc[l])*mu2
        X1    <- rnorm.aplus(n = n1, mean = mu1, var = CovStr)
        X2    <- rnorm.aplus(n = n2, mean = mu2, var = CovStr)
        X3    <- rnorm.aplus(n = n3, mean = mu3, var = CovStr)
        X     <- rbind(X1, X2, X3)
        fX    <- fpca.face(X)
        Y     <- fX$Yhat
        G     <- c(rep(0, n1), rep(1, n2), rep(2, n3))
        
        ###### doubly ranked test with geometric mean ######
        drt   <- kruskal.fda(Y ~ G, method = 'geo.rank')
        pvals[i, k, l, b] <- drt$p.value
        stats[i, k, l, b] <- drt$statistic
        
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
}

#### save output ####
filename <- paste('kw_pnorm_l_', n1, '_r_fpca.RData', sep = '')
save.image(file = filename)
