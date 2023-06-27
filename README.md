# doublyRanked
 Code to perform doubly ranked tests is in the file wilcox.fda.R. A description of the method is available on arXiv, [Meyer (2023)](https://arxiv.org/abs/2306.14761). Functions are described below. Folders contain code and data for the simulations and illustrations decribed in the manuscript.

 ## Doubly ranked MWW
 Perform doubly ranked Mann-Whitney-Wilcoxon test on functional data using
```
wilcox.fda(X ~ G)
```
Takes an $N\times S$ matrix class containing functional trajectories, X, and a grouping vector, G. Defaults to using geometric mean, a sufficient statistic, to summarize the ranks of curves over time as decribed in [Meyer (2023)](https://arxiv.org/abs/2306.14761). Arugments can have the formula class, X ~ G as shown above, or be of the form X = , G = .
```
wilcox.fda(X = X, G = G)
```
Both return objects of class htest.

 ## Doubly ranked KW
 Perform doubly ranked Kruskal-Wallis test on functional data using
```
kruskal.fda(X ~ G)
```
 Takes an $N\times S$ matrix class containing functional trajectories, X, and a grouping vector, G. Defaults to using geometric mean, a sufficient statistic, to summarize the ranks of curves over time as decribed in [Meyer (2023)](https://arxiv.org/abs/2306.14761). Arugments can have the formula class, X ~ G as shown above, or be of the form X = , G = .
```
kruskal.fda(X = X, G = G)
```
Both return objects of class htest.

## Geometric mean
 Take the geometric mean of a vector of data, x:
```
gmean(x)
```
For data with missing values, use 
```
gmean(x, na.rm = TRUE)
```
to remove and calculate the geometric mean.
