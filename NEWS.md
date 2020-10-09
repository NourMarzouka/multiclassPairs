# multiclassPairs v0.2.2 (Release date: 2020-10-09)

### Bug fixes:

* easier access to switchBox disjoint argument in train_one_vs_rest_TSP function
* Bug fix in plot_binary_TSP related to using ExpressionSet as input
* Bug fix related to passing additional arguments to SB training function by the user
* Update examples

# multiclassPairs v0.2.1 (Release date: 2020-09-28)

### Dependencies:

* Dependency issue solved (switchBox and Biobase packages are installed separately)

### Minor changes:
* Update examples.



# multiclassPairs v0.2.0 (Release date: 2020-09-24)

### Additions:

* additional function summary_genes_RF to summarize genes to rules stats.
* additional function optimize_RF to help in train_RF parameters optimization.

### Changes:

* plot_binary_RF now supports when RF model is trained with probability = FALSE.
* plot_binary_RF extracts prediction labels for training data from the classifier object.
* imputation is implemented in predict_RF function.
* NA is not allowed for class and platforms labels.

### Optimizations:

* stats for gene repetition in rules are stored in the sorted rules object to make training process faster

### Minor changes:
* Update examples.
* minor bug fixes.



# multiclassPairs v0.1.6 (Release date: 2020-09-08)

* first release on CRAN servers
