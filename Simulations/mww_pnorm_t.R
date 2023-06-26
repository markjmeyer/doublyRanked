#### load libraries ####
library(mvnfast)
library(refund)
library(numbers)
library(xtable)
source('~/Dropbox/Research/Drafts/Wilcoxon FDA/wilcox.fda.R')

#### specifications ####
# sn    <- c(1, 1.25, 1.5)
n1    <- n2 <- 10 # 10, 25, 50 -> 20, 50, 100
N     <- n1 + n2
S     <- c(40, 80, 160, 320)
sc    <- seq(0.99, 0.5, by = -0.02)
B     <- 500
sig   <- 0.05 # 0.5
rho   <- c(0.25, 0.5, 0.75)
# n1 <- 5; sn <- 1.25; S <- 160; sc <- 0.005
dots  <- 1000
up    <- 25000

# combinat2 <- function(n,p){
#   if (n<p){combinat=0}
#   else {combinat=exp(lfactorial(n)-(lfactorial(p)+lfactorial(n-p)))}
# }
# 
# fMBD2 <- function(data){
#   p     <- dim(data)[1]
#   n     <- dim(data)[2]
#   rmat  <- apply(data,1,rank)
#   down  <- rmat-1
#   up    <- n-rmat
#   (rowSums(up*down)/p+n-1)/combinat2(n,2)
# }

pvals   <- array(NA, dim = c(length(rho), length(S), length(sc), B, 3))
stats   <- array(NA, dim = c(length(rho), length(S), length(sc), B, 3))

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
        # n2    <- floor(sn[i]*n1)
        # N     <- n1 + n2
        s     <- seq(0, 1, length = S[k])
        mu1   <- pnorm(s, mean = 0.5, sd = 0.1)
        mu2   <- sc[l]*mu1
        X1    <- rmvt(n = n1, mu = mu1, sigma = CovStr, df = 2)
        X2    <- rmvt(n = n2, mu = mu2, sigma = CovStr, df = 2)
        X     <- rbind(X1, X2)
        # med   <- which(fMBD2(t(X)) == max(fMBD2(t(X))))
        # Xs    <- X - matrix(rep(X[med[1],], N), nrow = N, byrow = TRUE)
        # Xs2   <- Xs^2
        fX    <- fpca.face(X)
        Y     <- fX$Yhat
        G     <- c(rep(0, n1), rep(1, n2))
        
        ###### doubly ranked test ######
        drt   <- wilcox.fda(Y ~ G, method = 'geo.rank')
        pvals[i, k, l, b, 1] <- drt$wc.test$p.value
        stats[i, k, l, b, 1] <- drt$wc.test$statistic
        
        # fdrt  <- wilcox.fda(X ~ G, method = 'avg.rank')
        # pvals[i, k, l, b, 2] <- fdrt$wc.test$p.value
        # stats[i, k, l, b, 2] <- fdrt$wc.test$statistic
        
        ###### MBD rank sum test ######
        mrs   <- wilcox.fda(Y ~ G, method = 'mbd.rank')
        pvals[i, k, l, b, 2] <- mrs$wc.test$p.value
        stats[i, k, l, b, 2] <- mrs$wc.test$statistic
        
        ###### random projections based test ######
        rpt   <- wilcox.fda(X ~ G, method = 'rand.proj')
        pvals[i, k, l, b, 3] <- rpt$wc.test$p.value
        stats[i, k, l, b, 3] <- rpt$wc.test$statistic
        # ets <- proc.time() - sts
        
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
fname <- paste('/Users/mjm556/Dropbox/Research/Drafts/Wilcoxon FDA/Simulations/MVT/pnorm/MWW/N1_', n1, sep = '')
setwd(fname)
filename <- paste('mww_pnorm_t_', n1, '_r_fpca.RData', sep = '')
save.image(file = filename)

perMaxPow   <- array(NA, dim = c(3, length(S), length(rho)))
plotFlag    <- TRUE
for(i in 1:length(rho)){
  for(k in 1:length(S)){
    pow1  <- apply(pvals[i, k, , , 1] < 0.05, 1, mean)
    pow2  <- apply(pvals[i, k, , , 2] < 0.05, 1, mean)
    pow3  <- apply(pvals[i, k, , , 3] < 0.05, 1, mean)
    
    pow   <- rbind(pow1, pow2, pow3)
    powR  <- apply(pow, 2, rank)
    maxR  <- apply(powR, 2, max)
    
    maxP  <- matrix(NA, nrow = 3, ncol = length(sc))
    for(l in 1:length(sc)){
      maxP[,l] = 1*(powR[,l] == maxR[l])
    }
    perMaxPow[, k, i] <- apply(maxP, 1, mean)
    
    if(plotFlag){
      figname   <- paste('mww_phi_n1_', n1,'_n2_', n2, '_S', S[k], '_r', round(100*rho[i]),  '_t_fpca.pdf', sep = '')
      pdf(figname)
      plot(1-sc, apply(pvals[i, k, , , 1] < 0.05, 1, mean), type = 'n', ylim = c(0,1), main = paste('S = ', S[k], sep = ''), ylab = '',
           xlab = '', cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25)
      abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray')
      title(ylab = paste('r = ', rho[i], sep = ''), xlab = expression(xi), cex.lab = 1.5, line = 2)
      lines(1-sc, pow1, type = 'b', col = 'forestgreen', lty = 1, pch = 15, lwd = 2, cex = 1.25)
      lines(1-sc, pow2, type = 'b', col = 'mediumpurple3', lty = 2, pch = 16, lwd = 2, cex = 1.25)
      lines(1-sc, pow3, type = 'b', col = 'gold2', lty = 3, pch = 17, lwd = 2, cex = 1.25)
      # lines(1-sc, pow4, type = 'b', col = 'indianred3', lty = 4, pch = 18, lwd = 2, cex = 1.25)
      # lines(1-sc, pow4, type = 'b', col = 'royalblue3', lty = 5, pch = 3, lwd = 2, cex = 1.25)
      dev.off()
    }
  }
}

tab <- matrix(apply(apply(perMaxPow, c(1,2), mean), 1, mean), nrow = 1)
xtable(tab, digits = 4)

