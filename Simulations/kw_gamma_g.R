#### load libraries ####
library(refund)
library(numbers)
library(xtable)
library(mvnfast)
source('~/Dropbox/Research/Drafts/Wilcoxon FDA/wilcox.fda.R')

#### specifications ####
n1    <- n2   <-  n3 <- 10
N     <- n1 + n2 + n3
one   <- TRUE # TRUE, FALSE
S     <- c(40, 80, 160, 320)
sc    <- seq(0.99, 0.5, by = -0.02)
B     <- 500
sig   <- 0.05 # 0.5
rho   <- c(0.25, 0.5, 0.75)
# n1 <- 5; sn <- 1.25; S <- 160; sc <- 0.005
dots  <- 1000
up    <- 25000

pvals   <- array(NA, dim = c(length(rho), length(S), length(sc), B))
stats   <- array(NA, dim = c(length(rho), length(S), length(sc), B))

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
    for(l in 1:length(sc)){
      for(b in 1:B){
        # l <- b <- 1
        # sts <- proc.time()
        set.seed(b)
        s     <- seq(0, 1, length = S[k])
        mu1   <- dgamma(s, 10, 30)/max(dgamma(s, 10, 30))
        mu2   <- sc[l]*mu1
        mu3   <- ifelse(one, 1, sc[l])*mu2
        X1    <- rmvn(n = n1, mu = mu1, sigma = CovStr)
        X2    <- rmvn(n = n2, mu = mu2, sigma = CovStr)
        X3    <- rmvn(n = n3, mu = mu3, sigma = CovStr)
        X     <- rbind(X1, X2, X3)
        # med   <- which(fMBD2(t(X)) == max(fMBD2(t(X))))
        # Xs    <- X - matrix(rep(X[med[1],], N), nrow = N, byrow = TRUE)
        # Xs2   <- Xs^2
        fX    <- fpca.face(X)
        Y     <- fX$Yhat
        G     <- c(rep(0, n1), rep(1, n2), rep(2, n3))
        
        ###### doubly ranked test with geometric mean ######
        drt   <- kruskal.fda(Y ~ G, method = 'geo.rank')
        pvals[i, k, l, b] <- drt$wc.test$p.value
        stats[i, k, l, b] <- drt$wc.test$statistic
        
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
fname <- paste('~/Dropbox/Research/Drafts/Wilcoxon FDA/Simulations/Gaussian/gamma/KW/N1_', n1, sep = '')
setwd(fname)
filename <- paste('kw_gamma_g_', n1, '_r_fpca.RData', sep = '')
save.image(file = filename)

perMaxPow   <- array(NA, dim = c(2, length(S), length(rho)))
plotFlag    <- TRUE
for(i in 1:length(rho)){
  for(k in 1:length(S)){
    pow1  <- apply(pvals[i, k, , ] < 0.05, 1, mean)
    
    if(plotFlag){
      # if(one){
      #   figname   <- paste('kw_gam_n_', N, '_S', S[k], '_r', round(100*rho[i]),  '_g.pdf', sep = '')
      # } else{
      #   figname   <- paste('kw_gam_n_', N, '_S', S[k], '_r', round(100*rho[i]),  '_c_2.pdf', sep = '')
      # }
      # pdf(figname)
      plot(1-sc, apply(pvals[i, k, , ] < 0.05, 1, mean), type = 'n', ylim = c(0,1), main = paste('N = ', N, sep = ''), ylab = '',
           xlab = '', cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25)
      abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray')
      title(ylab = paste('S = ', S[k], sep = ''), xlab = expression(Delta), cex.lab = 1.5, line = 2)
      lines(1-sc, pow1, type = 'b', col = 'forestgreen', lty = 1, pch = 15, lwd = 2, cex = 1.25)
      # dev.off()
    }
  }
}

tab <- matrix(apply(apply(perMaxPow, c(1,2), mean), 1, mean), nrow = 1)
xtable(tab, digits = 4)

