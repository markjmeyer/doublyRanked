#### load libraries ####
library(refund)
library(numbers)
library(xtable)
library(mvnfast)
source('~/Dropbox/Research/Drafts/Wilcoxon FDA/wilcox.fda.R')

#### specifications ####
n1    <- n2   <-  n3 <- 50
N     <- n1 + n2 + n3
S     <- c(40, 80, 160, 320)
sc    <- 1
B     <- 1000
sig   <- 0.05 # 0.5
rho   <- c(0.25, 0.5, 0.75)
# n1 <- 5; sn <- 1.25; S <- 160; sc <- 0.005
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
      # b <- 1
      # sts <- proc.time()
      set.seed(b)
      s     <- seq(0, 1, length = S[k])
      mu1   <- dgamma(s, 10, 30)/max(dgamma(s, 10, 30))
      mu2   <- sc*mu1
      mu3   <- sc*mu1
      X1    <- rmvt(n = n1, mu = mu1, sigma = CovStr, df = 2)
      X2    <- rmvt(n = n2, mu = mu2, sigma = CovStr, df = 2)
      X3    <- rmvt(n = n3, mu = mu3, sigma = CovStr, df = 2)
      X     <- rbind(X1, X2, X3)
      # med   <- which(fMBD2(t(X)) == max(fMBD2(t(X))))
      # Xs    <- X - matrix(rep(X[med[1],], N), nrow = N, byrow = TRUE)
      # Xs2   <- Xs^2
      fX    <- fpca.face(X)
      Y     <- fX$Yhat
      G     <- c(rep(0, n1), rep(1, n2), rep(2, n3))
      
      ###### doubly ranked test with geometric mean ######
      drt   <- kruskal.fda(Y ~ G, method = 'geo.rank')
      pvals[i, k, b, 1] <- drt$wc.test$p.value
      stats[i, k, b, 1] <- drt$wc.test$statistic
      
      # fdrt  <- kruskal.fda(Xs2 ~ G, method = 'geo.rank')
      # pvals[i, k, b, 2] <- fdrt$wc.test$p.value
      # stats[i, k, b, 2] <- fdrt$wc.test$statistic
      
      ###### doubly ranked test with arithmetic mean ######
      # ars   <- kruskal.fda(X ~ G, method = 'avg.rank')
      # pvals[i, k, b, 3] <- ars$wc.test$p.value
      # stats[i, k, b, 3] <- ars$wc.test$statistic
      # 
      # frs   <- kruskal.fda(Xs2 ~ G, method = 'avg.rank')
      # pvals[i, k, b, 4] <- frs$wc.test$p.value
      # stats[i, k, b, 4] <- frs$wc.test$statistic
      
      ###### random projections based test ######
      rpt   <- kruskal.fda(X ~ G, method = 'rand.proj')
      pvals[i, k, b, 2] <- rpt$wc.test$p.value
      stats[i, k, b, 2] <- rpt$wc.test$statistic
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

#### save output ####
setwd('~/Dropbox/Research/Drafts/Wilcoxon FDA/Simulations/MVT/gamma/KW/Size')
filename <- paste('kw_pnorm_g_size_', n1, '_fpca.RData', sep = '')
save.image(file = filename)

tabM  <- apply(pvals < 0.05, c(1,2,4), mean)

# [c(1, 2, 4, 5, 3, 6)]
tab <- matrix(apply(apply(pvals < 0.05, c(1,2,4), mean), 3, mean), nrow = 1)
xtable(tab, digits = 4)

for(m in 1:3){
  print(xtable(tabM[,,m], digits = 4))
}

fname <- paste('kw_pnorm_', n1,'_size_c.pdf', sep = '')
pdf(fname)
boxplot(cbind(c(tabM[,,1]), c(tabM[,,2]) ), col = c('forestgreen', 'gold2'),
        names = c('Doubly Ranked', 'Random Projections'), main = paste('n1 =', n1), ylab = 'size')
dev.off()

