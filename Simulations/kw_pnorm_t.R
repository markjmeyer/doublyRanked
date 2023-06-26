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

pvals   <- array(NA, dim = c(length(rho), length(S), length(sc), B, 2))
stats   <- array(NA, dim = c(length(rho), length(S), length(sc), B, 2))

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
        mu1   <- pnorm(s, mean = 0.5, sd = 0.1)
        mu2   <- sc[l]*mu1
        mu3   <- ifelse(one, 1, sc[l])*mu2
        X1    <- rmvt(n = n1, mu = mu1, sigma = CovStr, df = 4)
        X2    <- rmvt(n = n2, mu = mu2, sigma = CovStr, df = 4)
        X3    <- rmvt(n = n3, mu = mu3, sigma = CovStr, df = 4)
        X     <- rbind(X1, X2, X3)
        # med   <- which(fMBD2(t(X)) == max(fMBD2(t(X))))
        # Xs    <- X - matrix(rep(X[med[1],], N), nrow = N, byrow = TRUE)
        # Xs2   <- Xs^2
        fX    <- fpca.face(X)
        Y     <- fX$Yhat
        G     <- c(rep(0, n1), rep(1, n2), rep(2, n3))
        
        ###### doubly ranked test with geometric mean ######
        drt   <- kruskal.fda(Y ~ G, method = 'geo.rank')
        pvals[i, k, l, b, 1] <- drt$wc.test$p.value
        stats[i, k, l, b, 1] <- drt$wc.test$statistic
        
        # fdrt  <- kruskal.fda(Xs2 ~ G, method = 'geo.rank')
        # pvals[i, k, l, b, 2] <- fdrt$wc.test$p.value
        # stats[i, k, l, b, 2] <- fdrt$wc.test$statistic
        
        ###### doubly ranked test with arithmetic mean ######
        # ars   <- kruskal.fda(X ~ G, method = 'avg.rank')
        # pvals[i, k, l, b, 3] <- ars$wc.test$p.value
        # stats[i, k, l, b, 3] <- ars$wc.test$statistic
        # 
        # frs   <- kruskal.fda(Xs2 ~ G, method = 'avg.rank')
        # pvals[i, k, l, b, 4] <- frs$wc.test$p.value
        # stats[i, k, l, b, 4] <- frs$wc.test$statistic
        
        ###### random projections based test ######
        rpt   <- kruskal.fda(X ~ G, method = 'rand.proj')
        pvals[i, k, l, b, 2] <- rpt$wc.test$p.value
        stats[i, k, l, b, 2] <- rpt$wc.test$statistic
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
fname <- paste('/Users/mjm556/Dropbox/Research/Drafts/Wilcoxon FDA/Simulations/MVT/pnorm/KW/N1_', n1, sep = '')
setwd(fname)
filename <- paste('kw_pnorm_t_', n1, '_r_fpca.RData', sep = '')
save.image(file = filename)

perMaxPow   <- array(NA, dim = c(2, length(S), length(rho)))
plotFlag    <- TRUE
for(i in 1:length(rho)){
  for(k in 1:length(S)){
    pow1  <- apply(pvals[i, k, , , 1] < 0.05, 1, mean)
    pow2  <- apply(pvals[i, k, , , 2] < 0.05, 1, mean)
    # pow3  <- apply(pvals[i, k, , , 3] < 0.05, 1, mean)
    # pow4  <- apply(pvals[i, k, , , 4] < 0.05, 1, mean)
    # pow5  <- apply(pvals[i, k, , , 5] < 0.05, 1, mean)
    # pow6  <- apply(pvals[i, k, , , 6] < 0.05, 1, mean)
    # pow7  <- apply(pvals[i, k, , , 7] < 0.05, 1, mean)
    
    pow   <- rbind(pow1, pow2)
    powR  <- apply(pow, 2, rank)
    maxR  <- apply(powR, 2, max)
    
    maxP  <- matrix(NA, nrow = 2, ncol = length(sc))
    for(l in 1:length(sc)){
      maxP[,l] = 1*(powR[,l] == maxR[l])
    }
    perMaxPow[, k, i] <- apply(maxP, 1, mean)
    
    if(plotFlag){
      # if(one){
      #   figname   <- paste('kw_phi_n_', N, '_S', S[k], '_r', round(100*rho[i]),  '_t.pdf', sep = '')
      # } else{
      #   figname   <- paste('kw_phi_n_', N, '_S', S[k], '_r', round(100*rho[i]),  '_c_2.pdf', sep = '')
      # }
      # pdf(figname)
      plot(1-sc, apply(pvals[i, k, , , 1] < 0.05, 1, mean), type = 'n', ylim = c(0,1), main = paste('N = ', N, sep = ''), ylab = '',
           xlab = '', cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25)
      abline(v = axTicks(1), h = axTicks(2), lty = 6, col = 'lightgray')
      title(ylab = paste('S = ', S[k], sep = ''), xlab = expression(Delta), cex.lab = 1.5, line = 2)
      lines(1-sc, pow1, type = 'b', col = 'forestgreen', lty = 1, pch = 15, lwd = 2, cex = 1.25)
      # lines(sc, pow2, type = 'b', col = 'mediumpurple3', lty = 2, pch = 16, lwd = 2, cex = 1.25)
      lines(1-sc, pow2, type = 'b', col = 'gold2', lty = 3, pch = 17, lwd = 2, cex = 1.25)
      # lines(sc, pow4, type = 'b', col = 'indianred3', lty = 4, pch = 18, lwd = 2, cex = 1.25)
      # lines(sc, pow5, type = 'b', col = 'royalblue3', lty = 5, pch = 3, lwd = 2, cex = 1.25)
      # lines(sc, pow6, type = 'b', col = 'orange', lty = 6, pch = 4, lwd = 2, cex = 1.25)
      # lines(sc, pow7, type = 'b', col = 'black', lty = 7, pch = 5, lwd = 2, cex = 1.25)
      # dev.off()
    }
  }
}

tab <- matrix(apply(apply(perMaxPow, c(1,2), mean), 1, mean), nrow = 1)
xtable(tab, digits = 4)

