# multiclassPairs v0.3.0 (Release date: 2020-11-11)

### changes
* replace the mode imputation method to kNN method in predict_RF function.
* train_RF now stores the whole binary matrix instead of mode matrix.
* change workflow figures in the tutorial.
* the predict_RF function can predict matrix with one sample without errors


# multiclassPairs v0.3.0 (Release date: 2020-11-02)

### changes:
* proximity_matrix_RF replaced cocluster_RF function and it can return and plot the proximity matrix

### Bug fixes:
* FIXED: plot_binary_RF does not get the predictions and scores when using as_training=TRUE and top_anno="platfrom" or "prediction"

# multiclassPairs v0.2.2 (Release date: 2020-10-09)

### Additions:
* Tutorial is available now.

### Minor changes:
* easier access to switchBox disjoint argument in train_one_vs_rest_TSP function.
* Update examples.

### Bug fixes:
* plot_binary_TSP when using ExpressionSet as input with no ref or platform.
* passing additional arguments to SB training function by the user.
* printing number of rules in the print function for sorted rules.
* border = NA instead of border = FALSE in plotting functions.
* optimize_RF can handle two classes problems without errors
* num.trees = num.trees missed in ranger for featureNo_altogether slots


# multiclassPairs v0.2.1 (Release date: 2020-09-28)

### Dependencies:

* Dependency issue solved (switchBox and Biobase packages are installed separately).

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

* stats for gene repetition in rules are stored in the sorted rules object to make training process faster.

### Minor changes:
* Update examples.
* minor bug fixes.



# multiclassPairs v0.1.6 (Release date: 2020-09-08)

* first release on CRAN servers
