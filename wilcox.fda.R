#### Rank-based FDA Tests ####
# by:   MJ Meyer
# Date: 06/08/23

#### general functions ####
gmean <- function(x, na.rm = TRUE){
  if(na.rm){
    xp <- x[complete.cases(x)]
    exp(mean(log(xp)))
  } else{
    exp(mean(log(x)))
  }
}

#### Wilcoxon tests for FDA ####
wilcox.fda  <- function(x, ...){
  UseMethod('wilcox.fda')
}

wilcox.fda.default  <- function(X, G, method = c('avg.rank', 'geo.rank'),
                                data.names = NULL, ...){
  
  # X is matrix class for functional trajectories
  # G is class vector/numeric, one dim must be 1
  
  
  if(length(method) > 1){
    method  <- 'geo.rank'
  }
  
  if(method != 'avg.rank' & method != 'geo.rank'){
    stop("method must be either avg.rank or geo.rank")
  }
  
  if(length(unique(G)) != 2){
    stop('wlicoxon.fda only works for two sample')
  }
  
  if(method == 'avg.rank'){
    
    ##### doubly ranked test with arithmetic mean #####
    R       <- apply(X, 2, rank, na.last = 'keep')
    Y       <- rowMeans(R, na.rm = TRUE)
    wtest   <- suppressWarnings(wilcox.test(Y ~ G, ...))
    
  } else{
    
    ##### doubly ranked test with geometric mean #####
    R       <- apply(X, 2, rank, na.last = 'keep')
    Y       <- apply(R, 1, gmean, na.rm = TRUE)
    wtest   <- suppressWarnings(wilcox.test(Y ~ G, ...))
    
  }
  
  stat.method <- switch(method,
                        avg.rank = 'Doubly Ranked Mann-Whitney-Wilcoxon Test with Arithmetic Mean',
                        geo.rank = 'Doubly Ranked Mann-Whitney-Wilcoxon Test with Geometric Mean'
  )
  
  dataList  <- list(X = X, G = G)
  if(is.null(data.names)){
    data.names <- names(dataList)
  }
  out       <- list(statistic = wtest$statistic, parameter = NULL, p.value = wtest$p.value,
                    null.value = wtest$null.value, alternative = wtest$alternative, 
                    data.name = paste(data.names, collapse = ' by '),
                    wc.test = wtest, method = stat.method, ranks = list(Y = Y, R = R),
                    data = dataList)    
  
  class(out) <- 'htest'
  
  return(out)
}

wilcox.fda.formula <- function(formula, ...){
  X           <- eval(formula[[2]])
  G           <- eval(formula[[3]])
  data.names  <- c(formula[[2]], formula[[3]])
  
  out     <- wilcox.fda(X, G, data.names = data.names, ...)

  return(out)
}


#### Kruskal-Wallis tests for FDA ####
kruskal.fda  <- function(x, ...){
  UseMethod('kruskal.fda')
}

kruskal.fda.default  <- function(X, G, method = c('avg.rank', 'geo.rank'),
                                 data.names = NULL, ...){
  
  # X is matrix class for functional trajectories
  # G is class vector/numeric, one dim must be 1
  
  
  if(length(method) > 1){
    method  <- 'geo.rank'
  }
  
  if(method != 'avg.rank' & method != 'geo.rank'){
    stop("method must be either avg.rank or geo.rank")
  }
  
  if(length(unique(G)) == 2){
    stop('use wlicoxon.fda for two sample')
  }
  
  N     <- nrow(X)
  S     <- ncol(X)
  s     <- seq(0, 1, length = S)
  
  if(method == 'avg.rank'){
    
    ##### doubly ranked test with arithmetic mean #####
    R       <- apply(X, 2, rank, na.last = 'keep')
    Y       <- rowMeans(R, na.rm = TRUE)
    wtest   <- suppressWarnings(kruskal.test(Y ~ G, ...))
    
  } else{
    
    ##### doubly ranked test with geometric mean #####
    R       <- apply(X, 2, rank, na.last = 'keep')
    Y       <- apply(R, 1, gmean)
    wtest   <- suppressWarnings(kruskal.test(Y ~ G, ...))
    
  }
  
  stat.method <- switch(method,
                        avg.rank = 'Doubly Ranked Kruskal-Wallis Test with Arithmetic Mean',
                        geo.rank = 'Doubly Ranked Kruskal-Wallis Test with Geometric Mean'
  )
  
  dataList  <- list(X = X, G = G)
  if(is.null(data.names)){
    data.names <- names(dataList)
  }
  out       <- list(statistic = wtest$statistic, parameter = NULL, p.value = wtest$p.value,
                    null.value = wtest$null.value, alternative = wtest$alternative, 
                    data.name = paste(data.names, collapse = ' by '),
                    wc.test = wtest, method = stat.method, ranks = list(Y = Y, R = R),
                    data = dataList)
  class(out) <- 'htest'
  
  return(out)
}

kruskal.fda.formula <- function(formula, ...){
  X           <- eval(formula[[2]])
  G           <- eval(formula[[3]])
  data.names  <- c(formula[[2]], formula[[3]])
  
  out     <- kruskal.fda(X, G, data.names = data.names, ...)
  
  return(out)
}
