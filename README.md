# doublyRanked
 Code to perform doubly ranked tests is in the file wilcox.fda.R. A description of the method is available on arXiv, [Meyer (2023)](https://arxiv.org/abs/2306.14761). Functions are described below. Folders contain code and data for the simulations and illustrations decribed in the manuscript.

 ## Doubly ranked MWW
 Performs doubly ranked Mann-Whitney-Wilcoxon test on functional data.
```
wilcox.fda(X ~ G)
```
Takes an $N\times S$ matrix class containing functional trajectories, X, and a grouping vector, G. Defaults to using geometric mean, a sufficient statistic, to summarize the ranks of curves over time as decribed in Meyer (2023). Arugments can be of the form X = , G = or X ~ G.

 ## kruskal.fda()
 Performs doubly ranked Kruskal-Wallis test on functional data. Takes an $N\times S$ matrix class containing functional trajectories, X, and a grouping vector, G. Defaults to using geometric mean, a sufficient statistic, to summarize the ranks of curves over time as decribed in Meyer (2023). Arugments can be of the form X = , G = or X ~ G.
