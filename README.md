# doublyRanked
 Code to perform doubly ranked tests is in the file wilcox.fda.R. A description of the method is available on arXiv, [Meyer (2023)](https://arxiv.org/abs/2306.14761). Functions are described below. Folders contain code and data for the simulations and illustrations, additional details below and in the manuscript.

## Functions
Functions in the file wilcox.fda.R include wilcox.fda() and kruskal.fda().

 ### Doubly ranked MWW
 Perform doubly ranked Mann-Whitney-Wilcoxon test on functional data using
```
wilcox.fda(X ~ G)
```
Takes an $N\times S$ matrix class containing functional trajectories, X, and a grouping vector, G. Defaults to using geometric mean, a sufficient statistic, to summarize the ranks of curves over time as decribed in [Meyer (2023)](https://arxiv.org/abs/2306.14761). Arugments can have the formula class, as shown above, or be of the form
```
wilcox.fda(X = X, G = G)
```
Both return objects of class htest.

 ### Doubly ranked KW
 Perform doubly ranked Kruskal-Wallis test on functional data using
```
kruskal.fda(X ~ G)
```
 Takes an $N\times S$ matrix class containing functional trajectories, X, and a grouping vector, G. Defaults to using geometric mean, a sufficient statistic, to summarize the ranks of curves over time as decribed in [Meyer (2023)](https://arxiv.org/abs/2306.14761). Arugments can have the formula class, as shown above, or be of the form
```
kruskal.fda(X = X, G = G)
```
Both return objects of class htest.

### Geometric mean
 Take the geometric mean of a vector of data, x:
```
gmean(x)
```
For data with missing values, use 
```
gmean(x, na.rm = TRUE)
```
to remove and calculate the geometric mean.

## Simulations folder
This folder contains files that were used for the power and type I error simulations described in [Meyer (2023)](https://arxiv.org/abs/2306.14761). The files are named with "mww" or "kw" first, depending on the test used. Then the true curve setting, "phi" or "gamma." Next is the the error distribution, "g" for Gaussian, "l" for log-normal, or "t" for t with 2 df. Files ending with "size" are for the type I error simulation.

## Illustrations folder
This folder contains a script and dataset for the data illustrations described in [Meyer (2023)](https://arxiv.org/abs/2306.14761). The script is titled Illustrations.R. The code in this file runs the doubly ranked tests on the data and produces graphics of the functional curves. Of the three data examples, two are available in packages on CRAN. The third is in the file driving_requests.txt, also in this folder.
