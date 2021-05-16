# multiclassPairs v0.4.3 (Release date: 2021-05-16)
### minor CRAN fixes

# multiclassPairs v0.4.1 (Release date: 2021-01-26)

### minor changes 
* minor change in rule_based_RandomForest print method
* default of k_range in train_one_vs_rest_TSP set to 10:50 instead of 2:50
* default of genes_altogether and genes_one_vs_rest in sort_rules_RF set to 50 instead of 200
* default of rules_altogether and rules_one_vs_rest in train_RF set to 50 instead of 200
* Update the tutorial with time and accuracy comparisons

# multiclassPairs v0.4.0 (Release date: 2020-11-19)

### changes
* train_RF has optimized gene_repetition method

# multiclassPairs v0.3.1 (Release date: 2020-11-16)

### changes
* replace the mode imputation method by kNN method in predict_RF function.
* train_RF now stores the whole binary matrix instead of mode matrix.
* change work-flow figures in the tutorial.
* the predict_RF function can predict matrix with one sample with no error


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
