#Function to read the data and Labels
ReadData <- function(Data,
                     Labels,
                     Platform = NULL,
                     verbose = TRUE) {

  # check the input Data format
  if (!is.data.frame(Data) &
      !is.matrix(Data) &
      class(Data)[1] != "ExpressionSet") {
    stop("Bad format for the input Data...should be:
         matrix, data.frame, or ExpressionSet")
  }

  if (is.data.frame(Data)) {
    Data_tmp <- Data
  }

  if (is.matrix(Data)) {
    Data_tmp <- as.data.frame(Data, stringsAsFactors = FALSE)
  }

  # if the input Data is ExpressionSet object
  if (class(Data)[1] == "ExpressionSet") {

    # extract the expression matrix from the ExpressionSet
    # Data_tmp <- as.matrix(exprs(Data))
    Data_tmp <- as.data.frame(exprs(Data), stringsAsFactors = FALSE)

    # if labels are not provided then give the available variables in Eset
    if (!hasArg(Labels)) {
      message(capture.output(cat("Phenotype data has these variables:",
                                 varLabels(Data),
                                 fill = TRUE)))
      stop("input a vector with same length of samples number
           or select one of these variable for Labels")
    }

    # extract the Labels - in case it is stored in the ExpressionSet
    if (is.character(Labels) & length(Labels) == 1) {

      if (Labels %in% varLabels(Data)) {
        Labels_tmp <- as.character(pData(Data)[, Labels])

      } else {
        message(capture.output(cat("Phenotype data has these variables:",
                                   varLabels(Data),
                                   fill = TRUE)))
        stop("Labels variable is not found in the phenotype data of your ExpressionSet")
      }
    }

    # get the input Labels vector as it is
    if ((is.character(Labels) | is.factor(Labels)) & length(Labels) != 1) {
      Labels_tmp <- as.character(Labels)

      if (length(Labels_tmp) != ncol(Data_tmp)) {
        message("Number of samples: ", ncol(Data_tmp))
        message("Labels length: ", length(Labels_tmp))
        stop("Labels vector length are not equal to samples in data")
      }
    }

    # if user input platform name or vector
    if (!is.null(Platform)) {
      # extract the Labels - in case it is stored in the ExpressionSet
      if (is.character(Platform) & length(Platform) == 1) {

        if (Platform %in% varLabels(Data)) {
          Platform_tmp <- as.character(pData(Data)[, Platform])

        } else {
          message(capture.output(
            cat("Phenotype data has these variables:",
                varLabels(Data), fill = TRUE)))
          stop("Platform variable is not found in the phenotype
               data of your ExpressionSet")
        }
      }

      # get the input Platform vector as it is
      if ((is.character(Platform) |
           is.factor(Platform)) &
          length(Platform) != 1) {

        Platform_tmp <- as.character(Platform)

        if (length(Platform_tmp) != ncol(Data_tmp)) {
          message("Number of samples:", ncol(Data_tmp))
          message("Labels length:", length(Platform_tmp))
          stop("Platform vector length are not equal to samples in data")
        }
      }
    }
  }

  # check if rownames is not NULL to avoid error later
  if (is.null(rownames(Data_tmp))) {
    stop("Provide features/genes names as rownames in the Data matrix!")
  }

  # get the input Labels vector as it is
  if ((is.character(Labels) |
       is.factor(Labels)) &
      class(Data)[1] != "ExpressionSet") {

    Labels_tmp <- as.character(Labels)

    if (length(Labels_tmp) != ncol(Data_tmp)) {
      message(paste("Number of samples: ", ncol(Data_tmp)))
      message(paste("Labels length: ", length(Labels_tmp)))
      stop("Labels vector length are not equal to samples in data")
    }
  }

  # get the input Platform vector as it is
  if (!is.null(Platform)) {
    if ((is.character(Platform) |
         is.factor(Platform)) &
        class(Data)[1] != "ExpressionSet") {

      Platform_tmp <- as.character(Platform)

      if (length(Platform_tmp) != ncol(Data_tmp)) {
        message(paste("Number of samples: ", ncol(Data_tmp)))
        message(paste("Labels length: ", length(Platform_tmp)))
        stop("Platform vector length are not equal to samples in data")
      }
    }
  } else {
    Platform_tmp <- NULL
  }

  ###
  # Remove genes with NAs in all samples
  remove_na <- rowSums(is.na(Data_tmp)) == ncol(Data_tmp)

  if (sum(remove_na) > 0) {
    message(paste("These features will be removed because they have NA values
                  in all samples:"))
    message(paste(rownames(Data_tmp)[remove_na], collapse = " "))
    Data_tmp <- Data_tmp[!remove_na, ]
  }

  # Remove genes with NAs in all genes
  remove_na <- colSums(is.na(Data_tmp)) == nrow(Data_tmp)

  if (sum(remove_na) > 0) {
    message(paste("These samples will be removed because they have NA values
                  for all features:"))
    message(paste(colnames(Data_tmp)[remove_na], collapse = " "))
    Data_tmp   <- Data_tmp[, !remove_na]
    Labels_tmp <- Labels_tmp[!remove_na]

    if (!is.null(Platform)) {
      Platform_tmp <- Platform_tmp[!remove_na]
    }
  }

  # print info about the input data and labels
  if (verbose) {
    message("Creating Data object...")
    message("Number of samples: ", ncol(Data_tmp))
    message("Number of genes/features: ", nrow(Data_tmp))
    message(capture.output(cat("Classes:", unique(Labels_tmp), fill = TRUE)))

    if (!is.null(Platform)) {
      message(capture.output(cat("Platforms/studies:",
                                 unique(Platform_tmp),
                                 fill = TRUE)))
    } else {
      message("Platforms/studies: NULL")
    }
  }

  # if any labels are NAs then stop
  if (any(is.na(as.character(Labels_tmp)))) {
    stop("NAs are not allowed in labels!")
  }

  if (any(is.na(as.character(Platform_tmp)))) {
    stop("NAs are not allowed in Platform!")
  }

  # create the object
  object <- list(
    data = list(Data=Data_tmp,
                Labels=Labels_tmp,
                Platform=Platform_tmp)
  )
  class(object) <- "multiclassPairs_object"

  return(object)
}

##### 1-vs-r TSP functions #####
# function for factorizing the labels
group_TSP <- function(label, my_group) {
  label    <- as.character(label)
  my_group <- as.character(my_group)
  label[ label != my_group] <- "rest"
  label[ label == my_group] <- my_group
  return(factor(label, levels=c(my_group, "rest")))
}

# dunn_test function
do_dunn_test <- function(data,
                         data_labels,
                         correction_method = "none") {

  sink(tempfile())
  on.exit(sink())

  if (is.null(rownames(data))) {
    stop("Provide rownames for Data matrix!")
  }

  p_values <- lapply(1:nrow(data), function(x) {
    dunn.test(
      x = as.numeric(data[x,]),
      g = data_labels,
      method = correction_method)
  })

  names(p_values) <- rownames(data)

  return(p_values)
}

# function to filter gene for SB
filter_genes_TSP <- function(data_object,
                             filter = c("one_vs_one", "one_vs_rest"),
                             platform_wise = FALSE,
                             featureNo = 1000,
                             UpDown = TRUE,
                             verbose = TRUE) {

  # check the data_object class
  if (class(data_object)[1] != "multiclassPairs_object") {
    stop("This function requires multiclassPairs_object!
              Use ReadData function to generate it.")
  }

  # check if platform_wise and UpDown are logical
  if (!is.logical(platform_wise) | !is.logical(UpDown)) {
    stop("platform_wise and UpDown arguments should be logical!")
  }

  # check the filter type
  if (length(filter) != 1 | !filter %in% c("one_vs_one","one_vs_rest")) {
    stop("Filter argument should be either: one_vs_one or one_vs_rest")
  }

  # check featureNo if is positive
  if (featureNo < 0 | !is.numeric(featureNo) | length(featureNo) != 1) {
    stop("filtered_genes argument should numeric with positive number")
  }

  # check if platform vector in the data_object
  if (is.null(data_object$data$Platform) & platform_wise == TRUE) {
    stop("platform_wise=TRUE while there is no platform vector in the object!")
  }

  # check the distribution of classes in the platforms
  if (platform_wise == TRUE) {
    # get the table for classes and platforms
    for_check <- as.data.frame.matrix(table(Platform = data_object$data$Platform,
                                            Classes  = data_object$data$Labels,
                                            useNA = "ifany"))

    # print the table for the user
    if (verbose) {
      message("Classes - Platform table:")
      print(for_check)
    }

    # give warning if some platforms lack some classes
    if (any(for_check == 0)) {
      message("Warning! All platforms/studies should have all classes to perform accurate platforms-wise filtering!")
      message("You can turn platform_wise to FALSE!")
    }

    # stop if there is any class found in any platform
    # but not in other platforms and with other classes
    # any class are not together with any other class in any platform then stop

    x <- c()
    for (y in 1:nrow(for_check)) {
      if (sum(for_check[y,] != 0) > 1) {
        x <- unique(c(x, colnames(for_check)[for_check[y,] != 0]))
      }
    }
    if (!all(colnames(for_check) %in% x)) {
      stop("Class:",colnames(for_check)[!colnames(for_check) %in% x],
           " not involved in any comparison! Turn platform_wise to FALSE!")
    }
  }

  # warning if the user need more genes than what is there in the data
  if (featureNo>nrow(data_object$data$Data)) {
    message("Warning!")
    message("featureNo argument > number of genes in your data")
    message("This could mean that all genes in dataset will be included!")
  }

  # create empty SB (one vs rest scheme) object
  if (verbose) {
    message("Creating new filtered genes object for one-vs-rest Scheme")
  }

  object_tmp <- list(
    OnevsrestScheme = list(filtered_genes=NULL,
                           calls=c()))
  class(object_tmp) <- "OnevsrestScheme_genes_TSP"

  # get the call
  param <- match.call()

  # get the data matrix
  D <- data_object$data$Data

  # get labels
  L <- data_object$data$Labels

  # get the classes
  groups <- unique(L)

  # get platforms vector
  if (platform_wise == TRUE) {
    plat_vector <- data_object$data$Platform
    studies     <- unique(plat_vector)
  }

  # prepare empty list for the genes
  filtered_genes <- vector("list", length(groups))
  names(filtered_genes) <- groups

  if (verbose) {
    message("filtering...")
  }

  # one_vs_rest filtering
  if (filter == "one_vs_rest" & platform_wise == FALSE) {

    for (i in groups) {
      if (verbose) {
        message(paste("Class: ",i))
      }
      # wilcoxon test from SwitchBox package
      filtered_genes[[i]] <- SWAP.Filter.Wilcoxon(
        inputMat = as.matrix(D),
        phenoGroup = group_TSP(label = L,
                               my_group = i),
        featureNo = featureNo,
        UpDown = UpDown)
    }

    # store the filtered genes for classes
    object_tmp$OnevsrestScheme$filtered_genes <- filtered_genes
    # save the call
    object_tmp$OnevsrestScheme$calls <- param

    if (verbose) {
      message("DONE!")
    }
    return(object_tmp)
  }

  # one_vs_rest filtering - Platform-wise
  if (filter == "one_vs_rest" & platform_wise == TRUE) {
    # wilcoxon test from SwitchBox package
    # make list for filtered genes for platform-wise
    plat_genes <- vector("list", length(studies))
    names(plat_genes) <- studies

    # fill each study with list of classes for genes
    for (y in studies) {
      plat_genes[[y]] <- filtered_genes
    }

    # this counter to follow how many platforms for each class
    counter <- filtered_genes

    # sort the genes in each platform and each class alone
    for(y in studies) {
      plat_samples <- plat_vector == y

      if (verbose) {
        message(paste("Platform/study: ",y))
      }

      for (i in groups) {
        if (verbose) {
          message(paste("Class: ",i))
        }
        # check if this class is there in this platform or not
        tmp_check <- as.character(group_TSP(label = L[plat_samples],
                                            my_group = i))

        if (length(unique(tmp_check)) == 1) {
          message("skip class ",i," for platform ",y)
          plat_genes[[y]] <- plat_genes[[y]][names(plat_genes[[y]]) != i]
          next
        }

        # wilcoxon test from SwitchBox package
        # this will return the same number of genes in D
        # this will sort the genes based on this platform
        plat_genes[[y]][[i]] <- SWAP.Filter.Wilcoxon(
          inputMat = as.matrix(D[,plat_samples]),
          phenoGroup = group_TSP(label = L[plat_samples],
                                 my_group = i),
          featureNo = nrow(D),
          UpDown = UpDown)

        # store the name this platform for that class
        counter[[i]] <- c(counter[[i]], y)
      }
    }

    if (verbose) {
      message("combining top filtered genes...")
    }

    # get the top genes
    if (UpDown == FALSE) {

      for (i in groups) {
        tmp <- as.vector(t(rbind(sapply(plat_genes[counter[[i]]], '[[', i))))

        # be sure it mentioned same to platforms number
        tmp2 <- table(tmp) == length(counter[[i]])#length(studies)
        tmp2 <- names(tmp2[tmp2 == TRUE])
        tmp  <- tmp[tmp %in% tmp2]

        # last mentioned in the list is last ranked
        tmp2 <- !duplicated(tmp, fromLast = TRUE)
        tmp <- tmp[tmp2][1:featureNo]

        # store it in the class genes
        filtered_genes[[i]] <- tmp[!is.na(tmp)]
      }
    }

    if (UpDown == TRUE) {
      for (i in groups) {
        # up genes
        tmp <- as.vector(t(rbind(sapply(plat_genes[counter[[i]]],
                                        function(x) {
                                          x[[i]][1:(nrow(D)/2)]
                                        }
        ))))
        # be sure it mentioned same to platforms number
        tmp2 <- table(tmp) == length(counter[[i]])#length(studies)
        tmp2 <- names(tmp2[tmp2 == TRUE])
        tmp  <- tmp[tmp %in% tmp2]

        # last mentioned in the list is last ranked
        tmp2 <- !duplicated(tmp, fromLast = TRUE)
        tmp_up <- tmp[tmp2][1:(featureNo/2)]
        rm(tmp, tmp2)

        # Down genes
        tmp <- as.vector(t(rbind(sapply(plat_genes[counter[[i]]],
                                        function(x) {
                                          x[[i]][nrow(D):(nrow(D)/2)]
                                        }
        ))))

        # be sure it mentioned same to platforms number
        tmp2 <- table(tmp) == length(counter[[i]])#length(studies)
        tmp2 <- names(tmp2[tmp2 == TRUE])
        tmp  <- tmp[tmp %in% tmp2]

        # last mentioned in the list is last ranked
        tmp2 <- !duplicated(tmp, fromLast = TRUE)
        tmp_down <- tmp[tmp2][1:(featureNo/2)]


        # store it in the class genes
        tmp <- c(tmp_up, tmp_down)
        filtered_genes[[i]] <- tmp[!is.na(tmp)]
      }
    }


    # store the filtered genes for classes
    object_tmp$OnevsrestScheme$filtered_genes <- filtered_genes
    # save the call
    object_tmp$OnevsrestScheme$calls <- param

    if (verbose) {
      message("DONE!")
    }
    return(object_tmp)
  }

  # one_vs_one filtering
  if (filter == "one_vs_one" & platform_wise == FALSE) {
    # rank data
    if (verbose) {
      message("ranking data...")
    }
    D_ranked <- apply(D,2,rank, ties.method = "min")

    # remove genes with only Zeros because it produce error with Dunn.test
    # because the zeros will be 1s for all samples (after ranking)
    D_ranked <- D_ranked[rowSums(D_ranked)>ncol(D_ranked),]

    # Perform Dunn test
    if (verbose) {
      message("Performing Dunn.test...")
    }
    tmp <- do_dunn_test(data = D_ranked, data_labels = L)

    # get the length of the pairwise comparison
    len_com <- length(tmp[[1]]$comparisons)

    # create empty df
    df <- data.frame(matrix(NA,
                            nrow = length(tmp)*len_com,
                            ncol = 5,
                            dimnames = list(c(),
                                            c("gene", "p_value", "Z",
                                              "comparison1", "comparison2"))))

    # to create the comparison1 and comparison2 column
    comparison1 <- sapply(strsplit(tmp[[1]]$comparisons, split = " - "),
                          function(x) {x[1]})
    comparison2 <- sapply(strsplit(tmp[[1]]$comparisons, split = " - "),
                          function(x) {x[2]})

    df$comparison1 <- rep(comparison1, times= length(tmp))
    df$comparison2 <- rep(comparison2, times= length(tmp))
    rm(comparison1, comparison2)

    # to create the gene column
    df$gene <- rep(names(tmp), each= len_com)

    # get the p values and the Z.statistics
    p_value <- Z <- c()

    for (g in names(tmp)) {
      p_value <- c(p_value, tmp[[g]]$P)
      Z       <- c(Z,       tmp[[g]]$Z)
    }

    df$p_value <- p_value
    df$Z       <- Z
    rm(p_value, Z)

    ###
    # sort based on the p-value
    tmp <- df
    tmp <- tmp[order(tmp$p_value, decreasing = FALSE),]

    # get equal number of genes from each side Up and Down
    if (UpDown == TRUE) {
      for (cl in groups) {
        if (verbose) {
          message(paste("Get genes for class: ",cl))
        }
        # create empty df for the up and down
        df <- data.frame(matrix(data = NA,
                                nrow = featureNo/2,
                                ncol = 2,
                                dimnames = list(c(),
                                                c("UP","DOWN"))))

        # get the classes lines
        tmp_cl <- tmp[tmp$comparison1 == cl|tmp$comparison2 == cl,]

        # NOTE: we convert the sign of Z without converting the comparison columns
        tmp_cl$Z[tmp_cl$comparison1 != cl] <- tmp_cl$Z[tmp_cl$comparison1 != cl]*-1

        # work on up genes
        tmp_cl_up <- tmp_cl[tmp_cl$Z>0,]
        tmp_cl_up <- tmp_cl_up[order(tmp_cl_up$p_value, -tmp_cl_up$Z,
                                     decreasing = FALSE),]

        for_count <- table(tmp_cl_up$gene) == (length(groups)-1) #max(table(tmp_cl$gene))
        for_count <- names(for_count[for_count == TRUE])
        tmp_cl_up <- tmp_cl_up[tmp_cl_up$gene %in% for_count,]

        if (nrow(tmp_cl_up)<((length(groups)-1)*(featureNo/2))) {
          message("Warning! Filtered genes will be less than featureNo!")
        }
        check <- !duplicated(tmp_cl_up$gene, fromLast = TRUE)
        df$UP <- tmp_cl_up$gene[check][1:(featureNo/2)]
        rm(for_count)

        # work on down genes
        tmp_cl_down <- tmp_cl[tmp_cl$Z<0,]
        tmp_cl_down <- tmp_cl_down[order(tmp_cl_down$p_value,
                                         tmp_cl_down$Z,
                                         decreasing = FALSE),]

        for_count <- table(tmp_cl_down$gene) == (length(groups)-1)#max(table(tmp_cl$gene))
        for_count <- names(for_count[for_count == TRUE])
        tmp_cl_down <- tmp_cl_down[tmp_cl_down$gene %in% for_count,]

        if (nrow(tmp_cl_down)<((length(groups)-1)*(featureNo/2))) {
          message("Warning! Filtered genes will be less than featureNo!")
        }

        check <- !duplicated(tmp_cl_down$gene, fromLast = TRUE)
        df$DOWN <- tmp_cl_down$gene[check][1:(featureNo/2)]

        final <- c(df$UP, df$DOWN)
        final <- final[!is.na(final)]

        # store the filtered genes for this class
        filtered_genes[[cl]] <- final
      }
    }

    # just get the top genes without considering which is up or down
    if (UpDown == FALSE) {
      for (cl in groups) {
        if (verbose) {
          message(paste("Get genes for class: ",cl))
        }

        # create empty df for the up and down
        df <- data.frame(matrix(data = NA,
                                nrow = featureNo,
                                ncol = 1,
                                dimnames = list(c(),
                                                c("UP_DOWN"))))

        # get the classes lines
        tmp_cl <- tmp[tmp$comparison1 == cl|tmp$comparison2 == cl,]

        # NOTE: we convert the sign of Z without converting the comparison columns
        tmp_cl$Z[tmp_cl$comparison1 != cl] <- tmp_cl$Z[tmp_cl$comparison1 != cl]*-1

        # work up and down genes based on the p-value then Z
        tmp_cl <- tmp_cl[order(tmp_cl$p_value, -abs(tmp_cl$Z),
                               decreasing = FALSE),]

        for_count <- table(tmp_cl$gene) == (length(groups)-1)#max(table(tmp_cl$gene))
        for_count <- names(for_count[for_count == TRUE])
        tmp_cl <- tmp_cl[tmp_cl$gene %in% for_count,]

        if (nrow(tmp_cl)<((length(groups)-1)*(featureNo))) {
          message("Warning! Filtered genes will be less than featureNo!")
        }

        check <- !duplicated(tmp_cl$gene, fromLast = TRUE)
        df$UP_DOWN <- tmp_cl$gene[check][1:(featureNo)]
        rm(for_count)

        final <- df$UP_DOWN
        final <- final[!is.na(final)]

        # store the filtered genes for this class
        filtered_genes[[cl]] <- final
      }
    }

    # store the filtered genes for classes
    object_tmp$OnevsrestScheme$filtered_genes <- filtered_genes
    # save the call
    object_tmp$OnevsrestScheme$calls <- param

    if (verbose) {
      message("DONE!")
    }
    return(object_tmp)
  }

  # one_vs_one filtering - Platform-wise
  if (filter == "one_vs_one" & platform_wise == TRUE) {
    # rank data
    if (verbose) {
      message("ranking data...")
    }
    D_ranked <- apply(D,2,rank, ties.method = "min")

    # make list for filtered genes for platform-wise
    plat_genes <- vector("list", length(studies))
    names(plat_genes) <- studies

    # to know check platforms has which classes
    tab_sam <- as.data.frame.matrix(table(L, plat_vector))
    involved_plat <- c()

    # sort the genes in each platform and each class alone
    for(y in studies) {
      plat_samples <- plat_vector == y
      if (verbose) {
        message(paste("Platform/study: ",y))
        message("Performing Dunn.test...")
      }
      # get data for this platform
      D_tmp <- D_ranked[,plat_samples]

      # remove genes with only Zeros because it produce error with Dunn.test
      # because the zeros will be 1s for all samples (after ranking)
      rank_same <- rowSums(D_tmp) == ncol(D_tmp)
      if (sum(rank_same)>0) {
        message(paste("These features will be removed because they rank the same in all samples from platform: ", y))
        message(paste(rownames(D_tmp)[rank_same], collapse = " "))
        D_tmp <- D_tmp[rowSums(D_tmp)>ncol(D_tmp),]
      }

      # Remove genes with NAs in all samples
      remove_na <- rowSums(is.na(D_tmp)) == ncol(D_tmp)
      if (sum(remove_na)>0) {
        message(paste("These features will be removed because they have NA values in all samples from platform: ", y))
        message(paste(rownames(D_tmp)[remove_na], collapse = " "))
        D_tmp <- D_tmp[!remove_na,]
      }

      # get which classes are there
      if (!length(unique(L[plat_samples])) > 1) {
        message("Platform ",y," was skipped due to lack of enough classes for dunn test")
        next
      }
      involved_plat <- c(involved_plat, y)

      # do dunn.test
      tmp <- do_dunn_test(
        data = D_tmp,
        data_labels = L[plat_samples])

      # get the length of the pairwise comparison
      len_com <- length(tmp[[1]]$comparisons)

      # create empty df
      df <- data.frame(matrix(NA,
                              nrow = length(tmp)*len_com,
                              ncol = 5,
                              dimnames = list(c(),
                                              c("gene", "p_value", "Z",
                                                "comparison1", "comparison2"))
      ))

      # to create the comparison1 and comparison2 column
      comparison1 <- sapply(strsplit(tmp[[1]]$comparisons, split = " - "),
                            function(x) {x[1]})
      comparison2 <- sapply(strsplit(tmp[[1]]$comparisons, split = " - "),
                            function(x) {x[2]})

      df$comparison1 <- rep(comparison1, times= length(tmp))
      df$comparison2 <- rep(comparison2, times= length(tmp))
      rm(comparison1, comparison2)

      # to create the gene column
      df$gene <- rep(names(tmp), each= len_com)

      # get the p values and the Z.statistics
      p_value <- Z <- c()

      for (g in names(tmp)) {
        p_value <- c(p_value, tmp[[g]]$P)
        Z       <- c(Z,       tmp[[g]]$Z)
      }

      # rank the p values
      df$p_value <- rank(p_value, ties.method = "min")
      df$Z       <- Z
      rm(p_value, Z)

      ###
      # sort based on the p-value
      df <- df[order(df$p_value, decreasing = FALSE),]
      plat_genes[[y]] <- df
    }

    if (verbose) {
      message("combining top filtered genes...")
    }
    # Put all the paltform-wise dfs in one df
    tmp <- do.call("rbind", plat_genes)

    # get equal number of genes from each side Up and Down
    if (UpDown == TRUE) {
      for (cl in groups) {
        if (verbose) {
          message(paste("Get genes for class: ",cl))
        }
        # create empty df for the up and down
        df <- data.frame(matrix(data = NA,
                                nrow = featureNo/2,
                                ncol = 2,
                                dimnames = list(c(),
                                                c("UP","DOWN"))))

        # get the classes lines
        tmp_cl <- tmp[tmp$comparison1 == cl | tmp$comparison2 == cl, ]

        if (nrow(tmp_cl)==0) {
          stop("class:",cl,"has no comparisons, try to switch platform-wise to FALSE!")
        }

        # NOTE: we convert the sign of Z without converting the comparison columns
        tmp_cl$Z[tmp_cl$comparison1 != cl] <- tmp_cl$Z[tmp_cl$comparison1 != cl]*-1

        # work on up genes
        tmp_cl_up <- tmp_cl[tmp_cl$Z>0,]
        tmp_cl_up <- tmp_cl_up[order(tmp_cl_up$p_value, -tmp_cl_up$Z,
                                     decreasing = FALSE),]

        # how many time the gene should be rep. in the df
        wanted_num  <- max(table(tmp_cl$gene))

        for_count <- table(tmp_cl_up$gene) == wanted_num
        for_count <- names(for_count[for_count == TRUE])
        tmp_cl_up <- tmp_cl_up[tmp_cl_up$gene %in% for_count,]

        if (nrow(tmp_cl_up)<(wanted_num*(featureNo/2))) {
          message("Warning! Filtered genes will be less than featureNo!")
        }
        check <- !duplicated(tmp_cl_up$gene, fromLast = TRUE)
        df$UP <- tmp_cl_up$gene[check][1:(featureNo/2)]
        rm(for_count)

        # work on down genes
        tmp_cl_down <- tmp_cl[tmp_cl$Z<0,]
        tmp_cl_down <- tmp_cl_down[order(tmp_cl_down$p_value,
                                         tmp_cl_down$Z,
                                         decreasing = FALSE),]

        for_count <- table(tmp_cl_down$gene) == wanted_num
        for_count <- names(for_count[for_count == TRUE])
        tmp_cl_down <- tmp_cl_down[tmp_cl_down$gene %in% for_count,]

        if (nrow(tmp_cl_down)<(wanted_num*(featureNo/2))) {
          message("Warning! Filtered genes will be less than featureNo!")
        }

        check <- !duplicated(tmp_cl_down$gene, fromLast = TRUE)
        df$DOWN <- tmp_cl_down$gene[check][1:(featureNo/2)]

        final <- c(df$UP, df$DOWN)
        final <- final[!is.na(final)]

        # store the filtered genes for this class
        filtered_genes[[cl]] <- final
      }
    }

    # just get the top genes without considering which is up or down
    if (UpDown == FALSE) {
      for (cl in groups) {
        if (verbose) {
          message(paste("Get genes for class: ",cl))
        }
        # create empty df for the up and down
        df <- data.frame(matrix(data = NA,
                                nrow = featureNo,
                                ncol = 1,
                                dimnames = list(c(),
                                                c("UP_DOWN"))))

        # get the classes lines
        tmp_cl <- tmp[tmp$comparison1 == cl|tmp$comparison2 == cl,]

        # NOTE: we convert the sign of Z without converting the comparison columns
        tmp_cl$Z[tmp_cl$comparison1 != cl] <- tmp_cl$Z[tmp_cl$comparison1 != cl]*-1

        # work up and down genes based on the p-value then Z
        tmp_cl <- tmp_cl[order(tmp_cl$p_value, -abs(tmp_cl$Z),
                               decreasing = FALSE),]

        # how many time the gene should be rep. in the df
        wanted_num    <- max(table(tmp_cl$gene))

        for_count <- table(tmp_cl$gene) == wanted_num
        for_count <- names(for_count[for_count == TRUE])
        tmp_cl <- tmp_cl[tmp_cl$gene %in% for_count,]

        if (nrow(tmp_cl)<(wanted_num*featureNo)) {
          message("Warning! Filtered genes will be less than featureNo!")
        }

        check <- !duplicated(tmp_cl$gene, fromLast = TRUE)
        df$UP_DOWN <- tmp_cl$gene[check][1:(featureNo)]
        rm(for_count)

        final <- df$UP_DOWN
        final <- final[!is.na(final)]

        # store the filtered genes for this class
        filtered_genes[[cl]] <- final
      }
    }

    # store the filtered genes for classes
    object_tmp$OnevsrestScheme$filtered_genes <- filtered_genes
    # save the call
    object_tmp$OnevsrestScheme$calls <- param

    if (verbose) {
      message("DONE!")
    }

    return(object_tmp)

  }
}


# train one vs rest scheme using SB
train_one_vs_rest_TSP <- function(data_object,
                                  filtered_genes,
                                  k_range=2:50,
                                  include_pivot=FALSE,
                                  one_vs_one_scores=FALSE,
                                  platform_wise_scores=FALSE,
                                  seed=NULL,
                                  classes,
                                  SB_arg = list(),
                                  verbose = TRUE) {

  # check the data_object class
  if (class(data_object)[1] != "multiclassPairs_object") {
    stop("This function requires multiclassPairs_object!
              Use ReadData function to generate it.")
  }

  # check the filtered_genes class
  if (class(filtered_genes)[1]  !=  "OnevsrestScheme_genes_TSP") {
    stop("filtered_genes should be OnevsrestScheme_genes_TSP object!
              Use filter_genes_TSP function to generate it!")
  }

  if (!is.logical(include_pivot) |
      !is.logical(one_vs_one_scores) |
      !is.logical(platform_wise_scores)) {
    stop("include_pivot, one_vs_one_scores, and platform_wise_scores should be logical (TRUE/FALSE) only!")
  }

  if (any(k_range<=0)) {
    stop("k_range should be positive number/range!")
  }

  # check if platform is available for platform_wise_scores
  if (platform_wise_scores == TRUE &
      is.null(data_object$data$Platform)) {
    stop("platform_wise_scores is TRUE while there is no platform information in data object!")
  }

  # get top genes
  genes <- filtered_genes[[1]]$filtered_genes

  # get the data matrix
  D <- data_object$data$Data

  # get labels
  L <- data_object$data$Labels

  # get platforms vector
  if (platform_wise_scores == TRUE) {
    plat_vector <- data_object$data$Platform
    studies     <- unique(plat_vector)
  }

  # get the classes
  if (hasArg(classes)) {
    groups <- classes

    # check if all classes are in the classifier object
    if (any(!classes %in% unique(L))) {
      message("These classes are not found in the data object:")
      message(classes[!classes %in% unique(L)])
      stop("classes argument must contain similar names to labels in data object!")
    }

    if (!all(unique(L) %in% classes)) {
      message("NOTE! Samples from these classes will be EXCLUDED:")
      message(capture.output(cat(unique(L)[!unique(L) %in% classes])))
      message("NOTE! Classifiers will be trained only for these classes:")
      message(capture.output(cat(classes[classes %in% unique(L)])))
      message("NOTE! If you want to include all classes then modify classes argument!")
    }

    # subset the data and vectors based on the classes
    D <- D[,L %in% groups]
    if (platform_wise_scores == TRUE) {
      plat_vector <- plat_vector[L %in% groups]
      studies     <- unique(plat_vector)
    }
    L <- L[L %in% groups]

  } else {
    # if no classes vector provided then take the unique of labels
    groups <- unique(L)
  }

  # create empty SB (one vs rest scheme) object
  if (verbose) {
    message("Creating new one-vs-rest scheme classifier")
  }
  object_tmp <- list(classifiers=NULL)
  class(object_tmp) <- "OnevsrestScheme_TSP"

  # train one vs rest classifier for each class
  for (cl in groups) {
    if (verbose) {
      message(paste("Class: ",cl))
    }
    # get the filtered gene
    genes_cl <- genes[[cl]]

    ### include pivot or not
    if (verbose){
      message("Pairing genes...")
    }

    if (include_pivot) {
      # restricted_pairs
      # (combine filtered genes with all genes on the dataset)
      restricted_pairs  <- as.matrix(expand.grid(genes_cl,
                                                 rownames(D),
                                                 stringsAsFactors = F))
      # remove when the same gene in both sides
      restricted_pairs <- restricted_pairs[restricted_pairs[,1]!=restricted_pairs[,2],]

      # keep only the unique rules (no need to have the same rule with swapped genes)
      restricted_pairs <- data.frame(restricted_pairs)
      restricted_pairs <- restricted_pairs[!duplicated(
        data.frame(list(do.call(pmin,restricted_pairs),
                        do.call(pmax,restricted_pairs)))),]
      restricted_pairs <- as.matrix(restricted_pairs)

      # all genes to be used as input
      genes_to_use <- rownames(D)

    } else {
      # restricted_pairs
      # (combine filtered genes with each other)
      restricted_pairs  <- as.matrix(expand.grid(genes_cl,
                                                 genes_cl,
                                                 stringsAsFactors = F))
      # remove when the same gene in both sides
      restricted_pairs <- restricted_pairs[restricted_pairs[,1]!=restricted_pairs[,2],]

      # keep only the unique rules (no need to have the same rule with swapped genes)
      restricted_pairs <- data.frame(restricted_pairs)
      restricted_pairs <- restricted_pairs[!duplicated(
        data.frame(list(do.call(pmin,restricted_pairs),
                        do.call(pmax,restricted_pairs)))),]
      restricted_pairs <- as.matrix(restricted_pairs)

      # filtered genes to be as input
      genes_to_use <- genes_cl
    }

    # create fake filtering function
    tmp_filter_fun <- function(phenoGroup, inputMat, ...) {
      return(genes_to_use)
    }


    ### score calculations
    if (verbose) {
      message("Score calculations...")
    }
    if (one_vs_one_scores == FALSE & platform_wise_scores == FALSE) {
      final_scores <- SWAP.Calculate.SignedTSPScores(
        classes = c(cl, "rest"),
        inputMat1  = as.matrix(D),
        phenoGroup = group_TSP(L, cl),
        RestrictedPairs = restricted_pairs,
        verbose = verbose)
    }


    if (one_vs_one_scores == FALSE & platform_wise_scores == TRUE) {
      # create empty list to store scores
      plat_scores <- vector("list", length(studies))
      names(plat_scores) <- studies

      for (ss in studies) {
        if (verbose) {
          message(paste("  Scores: ",cl,"in platform: ",ss, collapse = " "))
        }
        # skip calculations if there is no samples
        if (sum(plat_vector == ss & L == cl) == 0 |
            sum(plat_vector == ss & L != cl) == 0) {
          if (verbose) {
            message("Skip... due to lack of samples")
          }
          next
        }

        tmp_D <- D[,plat_vector == ss]
        tmp_L <- L[plat_vector == ss]

        tmp  <- SWAP.Calculate.SignedTSPScores(
          classes = c(cl, "rest"),
          inputMat1  = as.matrix(tmp_D),
          phenoGroup = group_TSP(tmp_L, cl),
          RestrictedPairs = restricted_pairs,
          verbose = verbose)

        plat_scores[[ss]] <- tmp$score
      }

      # take the mean of the scores
      if (verbose) {
        message("combine scores by mean")
      }
      # remove the skipped slots
      plat_scores[sapply(plat_scores, is.null)] <- NULL

      final_scores <- tmp
      final_scores$score <- rowMeans(sapply(plat_scores, unlist))
    }


    if (one_vs_one_scores == TRUE & platform_wise_scores == FALSE) {
      # create empty list to store scores
      groups_scores <- vector("list", (length(groups)-1))
      names(groups_scores) <- groups[!groups == cl]

      for (cl2 in groups[!groups == cl]) {
        if (verbose) {
          message(paste("  Scores:",cl,"vs",cl2, collapse = " "))
        }
        # skip calculations if there is no samples
        if (sum(L%in%cl) == 0 |
            sum(L%in%cl2) == 0) {
          if (verbose) {
            message("Skip... due to lack of samples")
          }
          next
        }

        tmp_D <- D[,L%in%c(cl,cl2)]
        tmp_L <- L[L%in%c(cl,cl2)]

        tmp  <- SWAP.Calculate.SignedTSPScores(
          classes = c(cl, "rest"),
          inputMat1  = as.matrix(tmp_D),
          phenoGroup = group_TSP(tmp_L, cl),
          RestrictedPairs = restricted_pairs,
          verbose = verbose)

        groups_scores[[cl2]] <- tmp$score
      }

      # take the mean of the scores
      if (verbose) {
        message("combine scores by mean")
      }
      # remove the skipped slots
      groups_scores[sapply(groups_scores, is.null)] <- NULL

      final_scores <- tmp
      final_scores$score <- rowMeans(sapply(groups_scores, unlist))
    }


    if (one_vs_one_scores == TRUE & platform_wise_scores == TRUE) {

      # create empty list to store scores
      groups_scores <- vector("list",
                              (length(groups)-1)*length(studies))
      # instead of names we will use counter
      counter <- 1

      for (ss in studies) {
        for (cl2 in groups[!groups == cl]) {
          if (verbose) {
            message(paste("  Scores:",cl,"vs",cl2,"in platform:",ss, collapse = " "))
          }
          # skip calculations if there is no samples
          if (sum(plat_vector == ss & L%in%cl) == 0 |
              sum(plat_vector == ss & L%in%cl2) == 0) {
            if (verbose) {
              message("Skip... due to lack of samples")
            }
            next
          }

          wanted_sam <- plat_vector == ss & L%in%c(cl,cl2)
          tmp_D <- D[,wanted_sam]
          tmp_L <- L[wanted_sam]

          tmp  <- SWAP.Calculate.SignedTSPScores(
            classes = c(cl, "rest"),
            inputMat1  = as.matrix(tmp_D),
            phenoGroup = group_TSP(tmp_L, cl),
            RestrictedPairs = restricted_pairs,
            verbose = verbose)

          groups_scores[[counter]] <- tmp$score
          counter <- counter + 1
        }
      }

      # take the mean of the scores
      if (verbose) {
        message("combine scores by mean")
      }
      # remove the skipped slots
      groups_scores[sapply(groups_scores, is.null)] <- NULL

      final_scores <- tmp
      final_scores$score <- rowMeans(sapply(groups_scores, unlist))
    }

    ### the fake function to input the calculated scores
    tmp_score_fun <- function(phenoGroup, inputMat1, ...) {
      return(final_scores)
    }

    # train the one-vs-rest classifier
    if (verbose) {
      message(paste("Train..."))
    }
    set.seed(seed)

    object_tmp[["classifiers"]][[cl]] <- SWAP.Train.KTSP(
      inputMat   = as.matrix(D),
      krange     = k_range,
      phenoGroup = group_TSP(L,cl),
      FilterFunc = tmp_filter_fun,
      score_fn   = tmp_score_fun,
      classes    = c(cl, "rest"),
      verbose    = verbose,
      unlist(SB_arg))
  }

  return(object_tmp)
}

# predict the classes based on one vs rest  switchBox classifier
predict_one_vs_rest_TSP <- function(classifier,
                                    Data,
                                    tolerate_missed_genes = TRUE,
                                    weighted_votes = TRUE,
                                    classes,
                                    verbose = TRUE) {

  # check the object class
  if (!class(Data)[1] %in% c("multiclassPairs_object",
                             "ExpressionSet",
                             "data.frame",
                             "matrix")) {
    stop("Data should be class:
    matrix/data.frame/ExpressionSet/multiclassPairs_object from ReadData function!")
  }

  # check classifier object
  if (class(classifier)[1] != "OnevsrestScheme_TSP") {
    stop("classifier should be OnevsrestScheme_TSP object from train_one_vs_rest_TSP function!")
  }

  # check the logical inputs
  if (!is.logical(tolerate_missed_genes) |
      !is.logical(weighted_votes)) {
    stop("tolerate_missed_genes and weighted_votes should be logical (TRUE/FALSE) only!")
  }

  # get the data matrix
  if (is.data.frame(Data)) {
    D <- Data
  }

  if (is.matrix(Data)) {
    D <- as.data.frame(Data, stringsAsFactors = FALSE)
  }

  # if the input Data is ExpressionSet object
  if (class(Data)[1] == "ExpressionSet") {
    # extract the expression matrix from the ExpressionSet
    #D <- as.matrix(exprs(Data))
    D <- as.data.frame(exprs(Data), stringsAsFactors = FALSE)
  }

  if (class(Data)[1]  ==  "multiclassPairs_object") {
    D <- Data$data$Data
    #D <- object$data$Data
  }

  if (hasArg(classes)) {
    classes <- classes

    # check if all classes are in the classifier object
    if (any(!classes %in% names(classifier[["classifiers"]]))) {
      message("These names are not found in the classifier object:")
      message(classes[!classes %in% names(classifier[["classifiers"]])])
      stop("classes names in classes argument should be similar to the names of the classifiers in classifier object!")
    }

  } else {
    # get the classes based on the names in the classifier object
    classes <- names(classifier[["classifiers"]])
  }

  # create empty score df
  scores_df <- data.frame(matrix(data = NA,
                                 nrow = ncol(D),
                                 ncol = length(classes)+2,
                                 dimnames = list(colnames(D),
                                                 c(classes,
                                                   "max_score",
                                                   "tie_flag"))),
                          stringsAsFactors = F,
                          check.names = FALSE,
                          check.rows = FALSE)

  for(cl in classes) {
    if (verbose) {
      message(paste("Get scores/votes from class: ",cl))
    }
    c <- classifier$classifiers[[cl]]

    # check if there is missed genes in the data
    # each rule has two genes, both of them should be in the sample data
    # otherwise this rule will be remove to avoid any errors
    sum  <- c$TSPs[,1] %in% rownames(D) + c$TSPs[,2] %in% rownames(D)

    if (tolerate_missed_genes == TRUE & any(sum != 2)) {
      message("These genes are missed in your data:")
      missed <- c$TSPs[,1][!c$TSPs[,1] %in% rownames(D)]
      missed <- c(missed, c$TSPs[,2][!c$TSPs[,2] %in% rownames(D)])
      message(missed)
      message("tolerate_missed_genes option is TRUE, so rules containing these genes will be Ignored")
      message(capture.output(cat(paste("Ignore",sum(sum != 2),
                                       "out of",length(sum),
                                       "rules for this class", collapse = " "))))

      # subset the classifier
      c$name  <- paste0("TSPs", sum(sum == 2))
      c$score <- c$score[sum == 2]
      c$TSPs  <- c$TSPs[sum == 2,]
    }

    if (tolerate_missed_genes == FALSE & any(sum != 2)) {
      message("Some gene are missed in your data!")
      missed <- c$TSPs[,1][!c$TSPs[,1] %in% rownames(D)]
      missed <- c(missed, c$TSPs[,2][!c$TSPs[,2] %in% rownames(D)])
      message(missed)
      stop("include these genes in your data or trun tolerate_missed_genes TRUE to exclude these rules from the classifier!")
    }

    # get scores for the predictions
    kappa <- SWAP.KTSP.Statistics(inputMat = as.matrix(D),
                                  classifier = c)

    if (weighted_votes) {
      scores_df[,cl]<-rowSums(t(t(kappa$comparisons)*c$score))/sum(c$score)
    } else {
      scores_df[,cl]<-rowMeans(kappa$comparisons)
    }
  }

  # get the prediction labels
  scores_df$max_score <- classes[max.col(scores_df[,classes],
                                         ties.method = "first")]

  # check the ties
  if (verbose) {
    message("Checking the ties")
  }
  first <- classes[max.col(scores_df[,classes], ties.method = "first")]
  last  <- classes[max.col(scores_df[,classes], ties.method = "last")]
  scores_df$tie_flag[first != last] <- "score_tie"

  # note to the use about the ties
  if (verbose) {
    if (sum(first != last)>0) {
      message(paste("Score ties found in",sum(first != last),
                    "out of",ncol(D),"samples in the data", collapse = " "))
      if (!weighted_votes) {
        message("Turning weighted_votes to TRUE may help to get less ties!")

      }
    } else {
      message("No ties found")
    }
  }

  class(scores_df) <- c(class(scores_df), "OneVsRestTSP prediction")
  return(scores_df)
}

# plot binary for classifier based on one vs rest scheme - switchBox
plot_binary_TSP <- function(Data,
                            classifier,
                            ref = NULL,
                            prediction = NULL,
                            platform = NULL,
                            classes = NULL,
                            platforms_ord = NULL,
                            top_anno = c("ref", "prediction", "platform")[1],
                            title = "",
                            binary_col = c("white", "black", "gray"),
                            ref_col = NULL,
                            pred_col = NULL,
                            platform_col = NULL,
                            show_ref = TRUE,
                            show_predictions = TRUE,
                            show_platform = TRUE,
                            show_scores = TRUE,
                            show_rule_name = TRUE,
                            legend = TRUE,
                            cluster_cols = TRUE,
                            cluster_rows = TRUE,
                            anno_height = 0.03,
                            score_height = 0.03,
                            margin = c(0, 5, 0, 5)) {


  ### get classifier ###
  # check classifier object
  if (class(classifier)[1] != "OnevsrestScheme_TSP") {
    stop("classifier should be OnevsrestScheme_TSP object from train_one_vs_rest_TSP function!")
  } else {
    C <- classifier
  }


  ### get data ###
  # check the object class
  if (!class(Data)[1] %in% c("multiclassPairs_object", "ExpressionSet",
                             "data.frame", "matrix")) {
    stop("Data should be class:
  matrix/data.frame/ExpressionSet/multiclassPairs_object from ReadData function!")
  }

  # get the data matrix
  if (is.data.frame(Data)) {
    D <- Data
  }

  if (is.matrix(Data)) {
    D <- as.data.frame(Data, stringsAsFactors = FALSE)
  }

  # if the input Data is ExpressionSet object
  if (class(Data)[1] == "ExpressionSet") {
    # extract the expression matrix from the ExpressionSet
    # D <- as.matrix(exprs(Data))
    D <- as.data.frame(exprs(Data), stringsAsFactors = FALSE)
  }

  if (class(Data)[1] == "multiclassPairs_object") {
    D <- Data$data$Data
  }

  # check if rownames is not NULL to avoid error later
  if (is.null(rownames(D))) {
    stop("Provide feature/gene names as rownames in the Data matrix!")
  }

  ### get classes ###
  if (!is.null(classes)) {
    # check if all classes are in the classifier object
    if (any(!classes %in% names(classifier[["classifiers"]]))) {
      message("These classes are not found in the classifier object:")
      message(paste0(classes[!classes %in% names(classifier[["classifiers"]])],
                     collapse = " "))
      stop("classes names in classes argument should be similar to the names of the classifiers in classifier object!")
    }
  } else {
    # get the classes based on the names in the classifier object
    classes <- names(classifier[["classifiers"]])
  }

  if (any(!names(classifier[["classifiers"]]) %in% classes)) {
    message("Because the classes argument miss these classes then these classes will be removed from the heatmap:")
    message(paste0(names(classifier[["classifiers"]])[!names(classifier[["classifiers"]])
                                                      %in% classes], collapse = " "))

  }

  ### get ref labels ###
  # if the data is object
  if (class(Data)[1] == "multiclassPairs_object") {
    # get the ref from the object
    if (is.null(ref)) {
      L <- Data$data$Labels
    }
  }

  # get the input ref Labels from the user
  if ((is.character(ref) | is.factor(ref)) & class(Data)[1] !=
      "ExpressionSet") {
    L <- as.character(ref)
  }

  # if the input Data is ExpressionSet object
  if (class(Data)[1] == "ExpressionSet") {

    # extract the Labels - in case it is stored in the
    # ExpressionSet
    if (is.character(ref) & length(ref) == 1) {
      if (ref %in% varLabels(Data)) {
        L <- as.character(pData(Data)[, ref])
      } else {
        message(capture.output(cat("Phenotype data has these variables:",
                                   varLabels(Data), fill = TRUE)))
        stop("Ref label variable is not found in the phenotype data of your ExpressionSet")
      }
    }

    # get the input Labels vector as it is
    if ((is.character(ref) | is.factor(ref)) & length(ref) !=
        1) {
      L <- as.character(ref)
    }
  }

  # check the length of the ref labels
  if (length(L) != ncol(D) & !is.null(ref)) {
    message("Number of samples: ", ncol(D))
    message("Labels length: ", length(L))
    stop("Labels vector length are not equal to
       samples in data")
  }

  # no ref labels if the user did not input ref and the input
  # is not multiclassPairs_object
  if (!is.null(ref) & class(Data)[1] != "multiclassPairs_object") {
    L <- NULL
  }


  ### get Platform labels ###
  # if the data is object
  if (class(Data)[1] == "multiclassPairs_object") {
    # get the platform from the object
    if (is.null(platform)) {
      P <- Data$data$Platform
    }
  }

  # get the input platform Labels from the user
  if ((is.character(platform) | is.factor(platform)) & class(Data)[1] !=
      "ExpressionSet") {
    P <- as.character(platform)
  }

  # if the input Data is ExpressionSet object
  if (class(Data)[1] == "ExpressionSet") {

    # extract the platform label - in case it is stored in the
    # ExpressionSet
    if (is.character(platform) & length(platform) == 1) {
      if (platform %in% varLabels(Data)) {
        P <- as.character(pData(Data)[, platform])
      } else {
        message(capture.output(cat("Phenotype data has these variables:",
                                   varLabels(Data), fill = TRUE)))
        stop("Platform/study label variable is not found in the phenotype data of your ExpressionSet")
      }
    }

    # get the input Labels vector as it is
    if ((is.character(platform) | is.factor(platform)) &
        length(platform) != 1) {
      P <- as.character(platform)
    }
  }

  # check the length of the platform labels
  if (length(P) != ncol(D) & !is.null(platform)) {
    message("Number of samples: ", ncol(D))
    message("Labels length: ", length(P))
    stop("Platform labels vector length are not equal to
       samples in data")
  }

  # no platform labels if the user did not input platform and
  # the input is not multiclassPairs_object
  if (!is.null(platform) & class(Data)[1] != "multiclassPairs_object") {
    P <- NULL
  }

  ### get platforms_ord ###
  if (!is.null(platforms_ord) & !is.null(P)) {
    # check if all platforms are in platforms_ord
    if (any(!platforms_ord %in% P)) {
      message("These platform/study in platforms_ord are not in found the platform labels:")
      message(platforms_ord[!platforms_ord %in% P])
      stop("platforms_ord argument should have similar names of the platforms/studies in the data!")
    }
  } else {
    # get the platforms_ord based on the names in the classifier
    # object
    platforms_ord <- unique(P)
  }

  ### get prediction ###
  # check if the prediction df is from the prediction function
  if (!is.null(prediction) & any(class(prediction) %in% "OneVsRestTSP prediction")) {
    pred <- prediction
  } else {
    pred <- NULL
  }

  # check the length of the platform labels
  if (!is.null(prediction)) {
    if (nrow(pred) != ncol(D)) {
      message("Number of samples in the data: ", ncol(D))
      message("Number of samples in the prediction: ", ncol(pred))
      stop("prediction should be for the same data!
     Use predict_one_vs_rest_TSP to generate it for this data!")
    }
  }

  ### checks ###
  if (is.null(pred) & is.null(L) & is.null(P)) {
    stop("No available ref, prediction, or platform labels!
     One of them atleast is needed.")
  }

  ### checks for top_anno ###
  # check if the top_anno labels are available
  if (top_anno == "ref" & is.null(L)) {
    stop("top annotation (top_anno) is ref while there is no ref labels available!")
  }
  if (top_anno == "ref" & !show_ref) {
    message("show_ref was turned to TRUE because top_anno is 'ref'!")
  }

  if (top_anno == "prediction" & is.null(pred)) {
    stop("top annotation (top_anno) is prediction while there is no prediction dataframe available! Use predict_one_vs_rest_TSP function to generate it!")
  }
  if (top_anno == "prediction" & !show_predictions) {
    message("show_predictions was turned to TRUE because top_anno is 'prediction'!")
  }

  if (top_anno == "platform" & is.null(P)) {
    stop("top annotation (top_anno) is platform while there is no platform labels available!")
  }
  if (top_anno == "platform" & !show_platform) {
    message("show_platform was turned to TRUE because top_anno is 'platform'!")
  }

  if (any(!top_anno %in% c("ref", "prediction", "platform")) |
      !is.character(top_anno) | length(top_anno) != 1) {
    stop("Top annotation argument should be character with one of these options:
     ref prediction platform")
  }

  ### title ###
  if (!is.character(title) | length(title) != 1) {
    stop("Title argument should be character input!")
  }

  ### binary heatmap colors ###
  if (!is.character(binary_col) | length(binary_col) != 3) {
    stop("binary_col should be character input with length of 3!
     Three colors are needed for rules with false, true and NAs.
     By default it is c('white','black','gray')")
  }

  ### ref anno colors ###
  if (show_ref & !is.null(L) & !is.null(ref_col)) {
    if (!is.character(ref_col) | any(!classes %in% names(ref_col))) {
      stop("ref_col should be named character vector for all classes!")
    }
  }

  ### pred anno colors ###
  if (show_predictions & !is.null(pred) & !is.null(pred_col)) {
    if (!is.character(pred_col) | any(!classes %in% names(pred_col))) {
      stop("pred_col should be named character vector for all classes!")
    }
  }

  ### platform anno colors ###
  if (show_platform & !is.null(P) & !is.null(platform_col)) {
    if (!is.character(platform_col) | any(!P %in% names(platform_col))) {
      stop("platform_col should be named character vector for all platforms/studies!")
    }
  }

  ### colors ###
  # determine the colors groups_col
  # thanks to https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
  xx_colors <- c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
                 "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabed4",
                 "#469990", "#dcbeff", "#9A6324", "#fffac8", "#800000",
                 "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9")
  xx_colors2 <- c("#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9",
                  "#469990", "#dcbeff", "#9A6324", "#fffac8", "#800000",
                  "#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
                  "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabed4")

  if (is.null(ref_col) & !is.null(L)) {
    if (length(classes)<20) {
      ref_col <- xx_colors[1:length(classes)]
    } else {
      ref_col <- sample(xx_colors, size = length(classes), replace = T)
    }
    names(ref_col) <- classes
  }
  if (is.null(pred_col) & !is.null(pred)) {
    if (length(classes)<20) {
      pred_col <- xx_colors[1:length(classes)]
    } else {
      pred_col <- sample(xx_colors, size = length(classes), replace = T)
    }
    names(pred_col) <- classes
  }
  if (is.null(platform_col) & !is.null(P)) {
    # xx_colors2 to give it a bit different colors than the classes
    if (length(platforms_ord)<20) {
      platform_col <- xx_colors2[1:length(platforms_ord)]
    } else {
      platform_col <- sample(xx_colors, size = length(platforms_ord), replace = T)
    }
    names(platform_col) <- platforms_ord
  }

  ### get info for top_anno ###
  # find the samples number to be used in the plotting
  # and get samples' names
  sam_names <- colnames(D)

  # get the labels and groups for the top anno
  if (top_anno == "ref") {
    lab        <- L
    groups     <- classes
    groups_col <- ref_col
  }
  if (top_anno == "prediction") {
    lab        <- pred$max_score
    groups     <- classes
    groups_col <- pred_col
  }
  if (top_anno == "platform") {
    lab        <- P
    groups     <- platforms_ord
    groups_col <- platform_col
  }

  ### anno ord ###
  anno_ord <- c("ref", "prediction", "platform")
  anno_ord <- anno_ord[c(!is.null(L) & show_ref,
                         !is.null(pred) & show_predictions,
                         !is.null(P) & show_platform)]
  anno_ord <- anno_ord[!anno_ord %in% top_anno]

  ### get sample order ###
  # cluster the samples in each group
  if (cluster_cols & (top_anno %in% c("ref", "prediction"))) {
    tmp <- c()
    for(i in groups){
      select_samples <- sam_names[lab==i]

      tmp_r <- C$classifiers[[i]]$TSPs

      tmp_binary <- D[tmp_r[,1],select_samples] > D[tmp_r[,2],select_samples]
      d   <- dist(t(tmp_binary[,select_samples]), method = "euclidean")
      fit <- hclust(d, method="ward.D2")
      tmp <- c(tmp, fit$labels[fit$order])
    }

    sam_ord <- order(match(sam_names, tmp))[1:length(tmp)]
    rm(tmp)
  }

  # cluster samples when platform is the top anno
  if (cluster_cols & top_anno == "platform") {

    # get the rules for all classifiers
    tmp_r <- data.frame(matrix(NA, ncol = 2, nrow = 0),
                        stringsAsFactors = FALSE)
    for (cl in classes) {
      tmp_r <- rbind(tmp_r, C$classifiers[[cl]]$TSPs)
    }

    tmp <- c()
    for(i in groups){
      select_samples <- sam_names[lab==i]

      tmp_binary <- D[tmp_r[,1],select_samples] > D[tmp_r[,2],select_samples]
      d   <- dist(t(tmp_binary[,select_samples]), method = "euclidean")
      fit <- hclust(d, method="ward.D2")
      tmp <- c(tmp, fit$labels[fit$order])
    }

    sam_ord <- order(match(sam_names, tmp))[1:length(tmp)]
    rm(tmp)
  }

  if (!cluster_cols) {
    # this will only group samples without clustering
    # based on the input data
    sam_ord <- order(match(lab, groups))
  }

  # change everything based on the new order
  # if null then will still be null
  D         <- D[,sam_ord]
  L         <- L[sam_ord]
  P         <- P[sam_ord]
  pred      <- pred[sam_ord,]
  lab       <- lab[sam_ord]
  sam_names <- sam_names[sam_ord]

  num_sam   <- ncol(D)

  # this should be after clustering find where the lines should
  # be the lines
  splits <- table(lab)[order(match(names(table(lab)), groups))]


  ### to keep the par settings from the user
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  ### plot top_anno ###
  {
    # Subtype annotation
    AreaStart <- 0.94
    SizeUnit <- anno_height
    Size <- SizeUnit * 1
    AreaEnd <- AreaStart - Size
    par(fig = c(0, 1, AreaEnd, AreaStart), mar = margin,
        mgp = c(3, 0.5, 0), new = FALSE)
    plot(c(0, 1), c(0, 1), type = "n", xaxs = "i", yaxs = "i",
         xlab = "", ylab = "", main = "", xlim = c(0, num_sam),
         ylim = c(0, 1), xaxt = "n", yaxt = "n", bty = "n")

    # headlines
    text_positions <- cumsum(splits)[1]/2
    for (i in 1:(length(cumsum(splits)) - 1)) {
      text_positions <- c(text_positions,
                          ((cumsum(splits)[i + 1] -
                              cumsum(splits)[i])/2 +
                             cumsum(splits)[i]))
    }

    # smaller headlines
    mtext(groups, side = 3, line = 0, outer = FALSE, at = text_positions,
          adj = NA, padj = NA, cex = 0.8, col = groups_col,
          font = NA)

    mtext(title, side = 3, line = -1, outer = TRUE, font = 2)


    # draw the subtypes
    axis(side = 2, at = 0.5,
         labels = c("ref"="Ref. labels",
                    "prediction"="Predictions",
                    "platform"="Platform/Study")[top_anno],
         las = 1, cex.axis = 0.7, tick = 0)
    for (f in groups) {
      for (g in which(lab == f)) {
        rect(g - 1, 0, g, 1, col = groups_col[f], border = NA)
      }
    }

    # the box and the white lines
    box(lwd = 1)
    li <- cumsum(splits)
    abline(v=li, lwd = 1.5, lty=1, col="black")
  }

  ### plot next annos ###
  for (i in anno_ord) {
    {
      # Subtype annotation
      Gap      <- 0.0
      AreaStart<- AreaStart-Size-Gap
      SizeUnit <- anno_height
      Size     <- SizeUnit*1
      AreaEnd  <- AreaStart-Size

      par(fig = c(0, 1, AreaEnd, AreaStart), mar = margin,
          mgp = c(3, 0.5, 0), new = TRUE)
      plot(c(0, 1), c(0, 1), type = "n", xaxs = "i", yaxs = "i",
           xlab = "", ylab = "", main = "", xlim = c(0, num_sam),
           ylim = c(0, 1), xaxt = "n", yaxt = "n", bty = "n")

      # draw the annotation name
      axis(side = 2, at = 0.5,
           labels = c("ref"="Ref. labels",
                      "prediction"="Predictions",
                      "platform"="Platform/Study")[i],
           las = 1, cex.axis = 0.7, tick = 0)

      if (i == "ref") {
        tmp_color <- ref_col
        tmp_lab   <- L
      }
      if (i == "prediction") {
        tmp_color <- pred_col
        tmp_lab   <- pred$max_score
      }
      if (i == "platform") {
        tmp_color <- platform_col
        tmp_lab   <- P
      }

      for (f in unique(tmp_lab)) {
        for (g in which(tmp_lab == f)) {
          rect(g - 1, 0, g, 1, col = tmp_color[f], border = NA)
        }
      }

      # the box and the white lines
      box(lwd = 1)
      li <- cumsum(splits)
      abline(v=li, lwd = 1.5, lty=1, col="black")
    }
  }

  ### from here if the top_anno is platform then groups should be classes
  if (top_anno == "platform") {
    groups <- classes
  }
  ### plot scores ###
  if (show_scores & !is.null(pred)){
    score_matrix <- pred[,groups, drop=FALSE]

    Gap      <- 0.01
    AreaStart<- AreaStart-Size-Gap
    SizeUnit <- score_height
    Size     <- SizeUnit*1
    AreaEnd  <- AreaStart-Size

    for (class in colnames(score_matrix)) {

      par(fig=c(0,1,AreaEnd,AreaStart),mar=margin,
          mgp=c(3,0.5,0),new=TRUE)

      barplot(as.numeric(score_matrix[,class]),
              col=pred_col[class],
              space=F,
              xaxs='i', yaxs='i',xlim=c(0,num_sam),border =NA,
              ylim=c(0,1),xaxt="n",yaxt="n",bty="n",ylab="",xlab="")
      box(lwd=1)
      axis(2, at =0.5,labels=paste("Scores:", class),las=1,cex.axis=0.7,tick=0)
      axis(4, at =c(0.1,0.5,0.9),labels=c("0", "0.5", "1"),
           las=1, cex.axis=0.4, tick = FALSE)

      li <- cumsum(splits)
      abline(v=li, lwd = 1.5, lty=1, col="black")

      ###
      if (class == colnames(score_matrix)[ncol(score_matrix)]) {
        next
      }
      ###
      Gap      <- 0.00
      AreaStart<- AreaStart-Size-Gap
      SizeUnit <- score_height
      Size     <- SizeUnit*1
      AreaEnd  <- AreaStart-Size
      ###
    }
  }

  ### plot binary heatmaps ###
  # to know the height of the heatmap
  Size <- (AreaEnd-(0.005*length(groups))-0.08)/length(groups)

  for (o in groups){
    tmp <- C$classifiers[[o]]$TSPs

    binary <- D[tmp[,1],] > D[tmp[,2],]
    binary <- binary + 1 # to fit with the indexes for the colors

    # cluster the rules
    if (cluster_rows) {
      d       <- dist(binary, method = "euclidean")
      fit     <- hclust(d, method="ward.D2")
      tmp_ord <- fit$labels[fit$order]
    } else {
      tmp_ord <- 1:nrow(binary)
    }

    binary <- binary[tmp_ord, ]

    ###
    Gap       <- 0.005
    AreaStart <- AreaEnd-Gap
    AreaEnd   <- AreaStart-Size
    ###

    par(fig = c(0, 1, AreaEnd, AreaStart),
        mar = margin, mgp = c(3, 0.5, 0), new=TRUE)

    myplot <- plot(c(0,1),c(0,1), type="n", xaxs='i', yaxs='i',
                   xlab = "", ylab = "", main = "",
                   xlim = c(0, num_sam), ylim = c(0, length(tmp_ord)),
                   xaxt = "n", yaxt = "n", bty = "n")

    if (show_rule_name){
      rule_names <- paste(tmp[,1],tmp[,2], sep = ">")
      axis(4, at =seq(0,(length(tmp_ord)-1),1)+0.5,
           labels=rule_names,las=1,cex.axis=0.5,tick=0)
    }

    for(f in 1:ncol(binary)){
      for(g in 1:nrow(binary)){
        rect(f-1,g,f,g-1,col=binary_col[binary[g,f]],border=F,lwd=0)
      }
    }

    # put the class name + number of rules
    axis(2, at = nrow(binary)/2, labels = paste(o, "\n", nrow(binary), "rules"),
         las = 1, cex.axis = 0.7, tick = 0, col = groups_col[o])

    box(lwd=1)

    li <- cumsum(splits)
    abline(v=li[-length(li)], lwd = 3, lty=3, col="red")
  }

  ### plot legends
  if (legend) {
    par(fig = c(0, 1, 0.02, (AreaEnd-0.01)),
        mar = margin, mgp = c(3, 0.5, 0), new=TRUE)
    plot(c(0,1),c(0,1), type="n", xaxs='i', yaxs='i',
         xlab = "", ylab = "", main = "",
         xlim = c(0, num_sam), ylim = c(0,1),
         xaxt = "n", yaxt = "n", bty = "n")

    if (!is.null(P) & show_platform) {
      legend(x = "topright", title = "Platform",
             ncol = length(platforms_ord), cex = 0.5,
             legend = platforms_ord,
             fill = platform_col)
    }

    if (!is.null(L) & show_ref) {
      title <- "Ref"
      if (!is.null(pred) & show_predictions &
          length(ref_col) == length(pred_col) &
          all(ref_col %in% pred_col)) {
        title <- "Classes"
      }
      legend(x = "topleft", title = title,
             ncol = length(classes), cex = 0.5,
             legend = names(ref_col),
             fill = ref_col)
    }

    if (!is.null(pred) & show_predictions &
        (length(ref_col) != length(pred_col) |
         any(!ref_col %in% pred_col))) {
      legend(x = "top", title = "Predictions",
             ncol = length(unique(pred$max_score)), cex = 0.5,
             legend = names(pred_col),
             fill = pred_col)
    }
  }
}

##### RF functions #####
# Gene filtering using RF
sort_genes_RF <- function (data_object,
                           featureNo_altogether,
                           featureNo_one_vs_rest,
                           rank_data = FALSE,
                           platform_wise = FALSE,
                           num.trees = 500,
                           min.node.size = 1,
                           importance = "impurity",
                           write.forest = FALSE,
                           keep.inbag = FALSE,
                           verbose = TRUE, ...) {

  ### Checks
  if (class(data_object)[1] != "multiclassPairs_object") {
    stop("This function requires multiclassPairs_object!
              Use ReadData function to generate it.")
  }
  if (!is.logical(platform_wise) |
      length(platform_wise) != 1 ){
    stop("platform_wise should be logical (TRUE/FALSE)!")
  }

  if (hasArg(featureNo_altogether)) {
    if (featureNo_altogether < 0 |
        !is.numeric(featureNo_altogether) |
        length(featureNo_altogether) != 1){
      stop("featureNo_altogether should be positive number!")
    }
  }

  if (hasArg(featureNo_one_vs_rest)) {
    if (featureNo_one_vs_rest < 0 |
        !is.numeric(featureNo_one_vs_rest) |
        length(featureNo_one_vs_rest) != 1){
      stop("featureNo_one_vs_rest should be positive number!")
    }
  }

  if (!importance[1] %in%  c('impurity', 'impurity_corrected', 'permutation')  |
      length(importance) > 1 |
      !is.character(importance)) {
    stop("importance variable should be one of these options 'impurity', 'impurity_corrected',  or 'permutation'!")
  }

  # check if platform vector in the data_object
  if (is.null(data_object$data$Platform) & platform_wise == TRUE) {
    stop("platform_wise=TRUE while there is no platform vector in the object!")
  }

  ### extract data and labels
  # get the data
  D <- data_object$data$Data

  # get the labels vector
  L <- data_object$data$Labels

  # get the classes
  groups <- unique(L)

  # get the call
  param <- match.call()

  # get platforms vector
  if (platform_wise == TRUE) {
    plat_vector <- data_object$data$Platform
    studies     <- unique(plat_vector)
  }

  # Remove genes with NAs
  if (sum(!complete.cases(D)) > 0) {
    message("These genes will be excluded from the analysis due to NAs:")
    message(capture.output(cat(rownames(D)[!complete.cases(D)])))
    D <- D[complete.cases(D),]
  }

  # rank data if user wants
  if (rank_data) {
    if (verbose) {
      message("Ranking data ...")
    }
    D <- apply(D,2,rank, ties.method = "min")
  }

  # check the distribution of classes in the platforms
  if (platform_wise == TRUE) {
    if (verbose) {
      message("Classes - Platform table:")
      print(table(Platform = data_object$data$Platform,
                  Classes  = data_object$data$Labels), right = FALSE)
      message()
    }
    for_check <- as.data.frame(table(Platform = data_object$data$Platform,
                                     Classes  = data_object$data$Labels))

    if (any(for_check$Freq == 0)) {
      message("platform-wise is TRUE!")
      message("Warning! Not all platforms/studies have all classes!")
      message("This could give unaccurate platforms-wise filtering!")
      message()
    }
  }

  # Warning if wanted filtered genes are larger than the available genes in data
  # if there is no specific number of genes return all genes sorted
  if (!hasArg(featureNo_altogether)) {
    featureNo_altogether <- nrow(D)
  }
  # if there is no specific number of genes return all genes sorted
  if (!hasArg(featureNo_one_vs_rest)) {
    featureNo_one_vs_rest <- nrow(D)
  }

  if (featureNo_altogether > nrow(D) |
      featureNo_one_vs_rest > nrow(D)) {
    message("Warning!")
    message("featureNo_altogether and/or featureNo_one_vs_rest arguments > number of genes in your data")
    message("This could mean that all genes in dataset will be included!")
    message()
  }

  ### Object
  # create empty Random Forests sorted genes object
  if (verbose) {
    message("Creating new Random Forests sorted genes object")
  }
  object_tmp <- list(
    RF_scheme = list(sorted_genes=NULL,
                     RF_classifiers=NULL,
                     calls=c()))
  class(object_tmp) <- "RandomForest_sorted_genes"
  message()


  ###
  # prepare empty list for the genes
  sorted_genes <- vector("list", length(groups)+1)
  names(sorted_genes) <- c("all", groups)

  ###
  # filtering
  if (verbose) {
    message("Building RF for filtering...")
  }
  # gene filtering without platform consideration
  if (platform_wise == FALSE) {

    # make list to store classifiers
    RF_classifiers <- sorted_genes

    if (featureNo_altogether > 0) {
      rf_all <- ranger(x=t(D),
                       y=factor(L),
                       num.trees = num.trees,
                       importance = importance,
                       write.forest = write.forest,
                       keep.inbag = keep.inbag,
                       min.node.size = min.node.size, ...)

      RF_classifiers[["all"]] <- rf_all

      if (verbose) {
        message(paste("RF: all classes",
                      "| num trees:", num.trees,
                      "| min node size:", min.node.size,
                      "| error:", round(rf_all$prediction.error,3), collapse = " "))
      }
      tmp1 <- sort(rf_all$variable.importance, decreasing = TRUE)
      tmp1 <- names(tmp1)
      tmp1 <- tmp1[1:featureNo_altogether]
      tmp1 <- tmp1[!is.na(tmp1)]
      sorted_genes[["all"]] <- tmp1
      rm(tmp1)
    }

    #Run the basic gene RF for subtypes separately (One vs Rest)
    if (featureNo_one_vs_rest > 0){

      for (cl in groups) {

        rfa_cl <- ranger(y = group_TSP(label = L, my_group = cl),
                         x = t(D),
                         num.trees = num.trees,
                         importance = importance,
                         write.forest = write.forest,
                         keep.inbag = keep.inbag,
                         min.node.size = min.node.size, ...)

        RF_classifiers[[cl]] <- rfa_cl

        if (verbose) {
          message(paste("RF: class", cl,
                        "| num trees:", num.trees,
                        "| min node size:", min.node.size,
                        "| error:", round(rfa_cl$prediction.error,3), collapse = " "))
        }
        tmp1 <- sort(rfa_cl$variable.importance, decreasing = TRUE)
        tmp1 <- names(tmp1)
        tmp1 <- tmp1[1:featureNo_one_vs_rest]
        tmp1 <- tmp1[!is.na(tmp1)]
        sorted_genes[[cl]] <- tmp1
      }
    }

    object_tmp$RF_scheme$sorted_genes <- sorted_genes
    object_tmp$RF_scheme$RF_classifiers <- RF_classifiers
    object_tmp$RF_scheme$calls          <- param
    return(object_tmp)
  }

  # platform-wise gene filtering
  if (platform_wise == TRUE) {

    # make list for sorted genes for platform-wise
    plat_genes <- vector("list", length(studies))
    names(plat_genes) <- studies

    # fill each study with list of classes for genes
    for (y in studies) {
      plat_genes[[y]] <- sorted_genes
    }

    # make list to store classifiers
    RF_classifiers <- plat_genes

    # this counter to follow how many platforms for each class
    counter <- sorted_genes
    #counter[sapply(counter, is.null)] <- c()

    # sort the genes in each platform and each class alone
    for(y in studies) {
      plat_samples <- plat_vector == y
      if (verbose) {
        message()
        message(paste("Platform/study: ",y))
      }
      if (featureNo_altogether > 0) {

        rf_all <- ranger(x=t(D[, plat_samples]),
                         y=factor(L[plat_samples]),
                         importance = importance,
                         write.forest = write.forest,
                         keep.inbag = keep.inbag,
                         min.node.size = min.node.size, ...)

        RF_classifiers[[y]][["all"]] <- rf_all

        if (verbose) {
          message(paste("RF: all classes",
                        "| num trees:", num.trees,
                        "| min node size:", min.node.size,
                        "| error:", round(rf_all$prediction.error,3), collapse = " "))
        }
        tmp1 <- sort(rf_all$variable.importance, decreasing = TRUE)
        tmp1 <- names(tmp1)
        plat_genes[[y]][["all"]] <- tmp1
      }

      #Run the basic gene RF for subtypes separately (One vs Rest)
      if (featureNo_one_vs_rest > 0){
        for (cl in groups){

          tmp_check <- as.character(group_TSP(label = L[plat_samples],
                                              my_group = cl))

          if (length(unique(tmp_check)) == 1) {
            if (verbose) {
              message("RF: skip class",cl,"in platform",y)
            }
            plat_genes[[y]] <- plat_genes[[y]][names(plat_genes[[y]]) != cl]
            next
          }

          rfa_cl <- ranger(y = group_TSP(label = L[plat_samples], my_group = cl),
                           x = t(D[, plat_samples]),
                           num.trees = num.trees,
                           importance = importance,
                           write.forest = write.forest,
                           keep.inbag = keep.inbag,
                           min.node.size=min.node.size, ...)

          RF_classifiers[[y]][[cl]] <- rfa_cl
          counter[[cl]] <- c(counter[[cl]], y)
          if (verbose) {
            message(paste("RF: class", cl,
                          "| num trees:", num.trees,
                          "| min node size:", min.node.size,
                          "| error:", round(rfa_cl$prediction.error,3), collapse = " "))
          }
          tmp1 <- sort(rfa_cl$variable.importance, decreasing = TRUE)
          tmp1 <- names(tmp1)
          plat_genes[[y]][[cl]] <- tmp1
        }
      }
    }

    # get the top genes in all platforms
    # fill altogether slot
    if (featureNo_altogether > 0){
      tmp <- as.vector(t(rbind(sapply(plat_genes, '[[', "all"))))

      # be sure it mentioned same to platforms number
      tmp2 <- table(tmp) == length(studies)
      tmp2 <- names(tmp2[tmp2 == TRUE])
      tmp  <- tmp[tmp %in% tmp2]

      # last mentioned in the list is last ranked
      tmp2 <- !duplicated(tmp, fromLast = TRUE)
      tmp <- tmp[tmp2][1:featureNo_altogether]

      # store it in the class genes
      sorted_genes[["all"]] <- tmp[!is.na(tmp)]
    } # otherwise it is already NULL


    # fill one vs rest for each class
    if (featureNo_one_vs_rest > 0) {

      for (cl in groups) {
        # get the sorted genes for that class from
        # platforms where this class is available
        tmp <- sapply(plat_genes[counter[[cl]]], '[[', cl)

        # get them in a vector (in the right order first first 2nd 2nd ...)
        tmp <- as.vector(t(rbind(tmp)))

        # be sure it mentioned same to platforms number
        tmp2 <- table(tmp) == length(counter[[cl]])
        tmp2 <- names(tmp2[tmp2 == TRUE])
        tmp  <- tmp[tmp %in% tmp2]

        # last mentioned in the list is last ranked
        tmp2 <- !duplicated(tmp, fromLast = TRUE)
        tmp <- tmp[tmp2][1:featureNo_one_vs_rest]

        # store it in the class genes
        sorted_genes[[cl]] <- tmp[!is.na(tmp)]
      }
    } # otherwise it is already NULL for each class

    # fill the object
    object_tmp$RF_scheme$sorted_genes <- sorted_genes
    object_tmp$RF_scheme$RF_classifiers <- RF_classifiers
    object_tmp$RF_scheme$calls          <- param
    return(object_tmp)
  }

  # end
}

# after sorting genes RF this gives an idea of how many genes you need to use to generate specific number of rules (NOTE without consideration of gene replication in rules, this need to be done after sorting the rules)
summary_genes_RF <- function(sorted_genes_RF,
                             genes_altogether,
                             genes_one_vs_rest) {

  if (class(sorted_genes_RF)[1] != "RandomForest_sorted_genes") {
    stop("This function requires RandomForest_sorted_genes object!
              Use sort_genes_RF function to generate it.")
  }

  if (any(genes_altogether < 0) |
      !is.numeric(genes_altogether) |
      length(genes_altogether) != length(genes_one_vs_rest)){
    stop("genes_altogether should a vector with zero or positive numbers and with the same length of genes_one_vs_rest vector!")
  }

  if (any(genes_one_vs_rest < 0) |
      !is.numeric(genes_one_vs_rest) |
      length(genes_altogether) != length(genes_one_vs_rest)){
    stop("genes_one_vs_rest should a vector with zero or positive numbers and with the same length of genes_altogether vector!")
  }

  n_all <- length(sorted_genes_RF[[1]]$sorted_genes$all)
  n_1_r <- length(sorted_genes_RF[[1]]$sorted_genes[[2]])

  # create empty df to store the results

  results_df <- data.frame(
    matrix(data = NA,
           nrow = length(genes_altogether),
           ncol = 7,
           dimnames = list(c(1:length(genes_altogether)),
                           c("genes_altogether_in_object",
                             "genes_1_vs_r_in_object",
                             "n_classes",
                             "from_altogether",
                             "from_one_vs_rest",
                             "n_unique_genes",
                             "n_rules"))),
    stringsAsFactors = FALSE)

  results_df[ , "genes_altogether_in_object"] <- n_all
  results_df[ , "genes_1_vs_r_in_object"]     <- n_1_r

  # prepare empty list for the genes
  empty_list_copy <- vector("list", length(sorted_genes_RF$RF_scheme$sorted_genes))
  names(empty_list_copy) <- names(sorted_genes_RF$RF_scheme$sorted_genes)

  # get the classes
  groups <- names(sorted_genes_RF$RF_scheme$sorted_genes)[-1]

  results_df[ , "n_classes"] <- length(groups)

  for (i in 1:length(genes_altogether)) {

    empty_list <- empty_list_copy

    results_df[i,"from_altogether"]  <- genes_alt    <- genes_altogether[i]
    results_df[i,"from_one_vs_rest"] <- genes_1_vs_R <- genes_one_vs_rest[i]

    if (!is.null(sorted_genes_RF[[1]]$sorted_genes$all)){
      if (genes_alt > 0) {
        tmp <- sorted_genes_RF[[1]]$sorted_genes$all
        tmp <- tmp[1:genes_alt]
        empty_list$all <- tmp[!is.na(tmp)]
      }
    }

    if (genes_1_vs_R > 0) {
      for (cl in groups) {
        if (is.null(sorted_genes_RF[[1]]$sorted_genes[[cl]])){
          next
        }
        tmp <- sorted_genes_RF[[1]]$sorted_genes[[cl]]
        tmp <- tmp[1:genes_1_vs_R]
        empty_list[[cl]] <- tmp[!is.na(tmp)]
      }
    }

    all_genes <- unique(as.vector(unlist(empty_list)))

    results_df[i,"n_unique_genes"] <- length(all_genes)
    results_df[i,"n_rules"]        <- choose(length(all_genes),2)
  }
  return(results_df)
}

# rules filtering using RF
sort_rules_RF <- function (data_object,
                           sorted_genes_RF,
                           genes_altogether = 200,
                           genes_one_vs_rest = 200,
                           run_altogether = TRUE,
                           run_one_vs_rest = TRUE,
                           platform_wise = FALSE,
                           num.trees = 500,
                           min.node.size = 1,
                           importance = "impurity",
                           write.forest = FALSE,
                           keep.inbag = FALSE,
                           verbose = TRUE, ...) {

  ### Checks
  if (class(data_object)[1] != "multiclassPairs_object") {
    stop("This function requires multiclassPairs_object!
              Use ReadData function to generate it.")
  }

  if (class(sorted_genes_RF)[1] != "RandomForest_sorted_genes") {
    stop("This function requires RandomForest_sorted_genes object!
              Use sort_genes_RF function to generate it.")
  }

  if (!is.logical(platform_wise) |
      length(platform_wise) != 1){
    stop("platform_wise should be logical (TRUE/FALSE)!")
  }

  if (genes_altogether < 0 |
      !is.numeric(genes_altogether) |
      length(genes_altogether) != 1){
    stop("genes_altogether should be zero or positive number!")
  }

  if (genes_one_vs_rest < 0 |
      !is.numeric(genes_one_vs_rest) |
      length(genes_one_vs_rest) != 1){
    stop("genes_one_vs_rest should be zero or positive number!")
  }

  if (!importance[1] %in%  c('impurity', 'impurity_corrected', 'permutation') |
      length(importance) > 1 |
      !is.character(importance)) {
    stop("importance variable should be one of these options 'impurity', 'impurity_corrected',  or 'permutation'!")
  }

  # check if platform vector in the object
  if (is.null(data_object$data$Platform) & platform_wise == TRUE) {
    stop("platform_wise=TRUE while there is no platform vector in the object!")
  }

  # checks for run_altogether and run_one_vs_rest
  if (!is.logical(run_altogether) |
      length(run_altogether) != 1){
    stop("run_altogether argument should be logical (TRUE/FALSE)!")
  }

  if (!is.logical(run_one_vs_rest) |
      length(run_one_vs_rest) != 1){
    stop("run_one_vs_rest argument should be logical (TRUE/FALSE)!")
  }

  if (!run_altogether & !run_one_vs_rest) {
    stop("run_altogether and run_one_vs_rest can not be both FALSE!
    One or both of them can be TRUE!")
  }

  ### extract data and labels
  # get the data
  D <- data_object$data$Data

  # get the labels vector
  L <- data_object$data$Labels

  # get the classes
  groups <- unique(L)

  # get the call
  param <- match.call()

  # get platforms vector
  if (platform_wise == TRUE) {
    plat_vector <- data_object$data$Platform
    studies     <- unique(plat_vector)
  }

  ### additional checks
  # Warning if wanted filtered genes > than the available genes
  if (genes_altogether > length(sorted_genes_RF[[1]]$sorted_genes$all) &
      run_altogether == TRUE) {
    message("NOTE!")
    message("genes_altogether > number of genes in sorted genes object")
    message("This means all genes in 'all' slot will be used!")
    message()
  }

  # Warning if wanted sorted genes > than the available genes
  tmp <- sapply(sorted_genes_RF[[1]]$sorted_genes[groups], length)
  if (any(genes_one_vs_rest > tmp) &
      run_one_vs_rest == TRUE) {
    message("NOTE!")
    message("genes_one_vs_rest > number of genes in sorted genes object")
    message("This means all genes will be used for these classes:")
    message(capture.output(cat(names(tmp)[genes_one_vs_rest > tmp])))
    message()
  }
  rm(tmp)

  # Remove genes with NAs
  if (sum(!complete.cases(D)) > 0) {
    message("These genes will be excluded from the analysis due to NAs:")
    message(capture.output(cat(rownames(D)[!complete.cases(D)])))
    D <- D[complete.cases(D),]
  }

  # check the distribution of classes in the platforms
  if (platform_wise == TRUE) {
    if (verbose) {
      message("Classes - Platform table:")
      print(table(Platform = data_object$data$Platform,
                  Classes  = data_object$data$Labels), right = FALSE)
      message()
    }
    for_check <- as.data.frame(table(Platform = data_object$data$Platform,
                                     Classes  = data_object$data$Labels))

    if (any(for_check$Freq == 0)) {
      message("platform-wise is TRUE!")
      message("Warning! Not all platforms/studies have all classes!")
      message("This could give unaccurate platforms-wise filtering!")
      message()
    }
  }

  ### Object
  # create empty Random Forests sorted genes object
  if (verbose) {
    message("Creating new Random Forests sorted rules object")
    message()
  }
  object_tmp <- list(
    RF_scheme = list(sorted_genes=NULL,
                     sorted_rules=NULL,
                     gene_repetition=NULL,
                     RF_classifiers=NULL,
                     calls=c()))
  class(object_tmp) <- "RandomForest_sorted_rules"

  ###
  # prepare empty list for the genes
  sorted_genes <- vector("list", length(groups) + 1)
  names(sorted_genes) <- c("all", groups)

  empty_list <- sorted_genes

  ###
  # get the wanted genes
  if (verbose) {
    message("Extracting the needed genes from the sorted genes object...")
  }
  # subset the sorted genes based on the needed number of genes
  if (!is.null(sorted_genes_RF[[1]]$sorted_genes$all)){
    if (genes_altogether > 0) {
      tmp <- sorted_genes_RF[[1]]$sorted_genes$all
      tmp <- tmp[1:genes_altogether]
      sorted_genes$all <- tmp[!is.na(tmp)]
    }
  }

  if (genes_one_vs_rest > 0) {
    for (cl in groups) {
      if (is.null(sorted_genes_RF[[1]]$sorted_genes[[cl]])){
        next
      }
      tmp <- sorted_genes_RF[[1]]$sorted_genes[[cl]]
      tmp <- tmp[1:genes_one_vs_rest]
      sorted_genes[[cl]] <- tmp[!is.na(tmp)]
    }
  }

  # merge sorted genes
  if (verbose) {
    message("Merge sorted genes (all + classes)...")
  }
  all_genes <- unique(as.vector(unlist(sorted_genes)))
  if (verbose) {
    message(paste(length(all_genes),"unique gene will be used for rules production",
                  collapse = " "))
  }
  # check if there is enough genes
  if (length(all_genes) < 2) {
    stop("No enough genes to make any rules!
         Check the arguments and the sorted genes object!")
  }

  if (length(all_genes) < 2*length(groups)) {
    message("Maybe there is no enough genes to make rules to seperate all classes!\nCheck the arguments and the sorted genes object!\nTry to increase the number of genes to get enough rules!")
  }

  # combine all sorted genes
  if (verbose) {
    message("Combine genes to produce all possible pairs...")
  }
  pairs     <- combn(all_genes, 2)
  if (verbose) {
    message(paste(length(pairs)/2,"pairs", collapse = " "))
    message()
  }
  # get binary matrix for training RF
  binary    <- D[pairs[1,],] < D[pairs[2,],]
  rownames(binary) <- paste0(pairs[1,],"__",pairs[2,])

  # filtering
  if (verbose) {
    message("Building RF for filtering...")
  }
  # gene filtering without platform consideration
  if (platform_wise == FALSE) {

    # make list to store classifiers and rules
    RF_classifiers   <- empty_list
    sorted_rules     <- empty_list
    gene_repetition  <- empty_list

    # run RF for all classes together
    if (run_altogether) {

      rf_all <- ranger(x=t(binary),
                       y=factor(L),
                       num.trees = num.trees,
                       importance = importance,
                       write.forest = write.forest,
                       keep.inbag = keep.inbag,
                       min.node.size = min.node.size, ...)

      RF_classifiers[["all"]] <- rf_all
      if (verbose) {
        message(paste("RF: all classes",
                      "| num trees:", num.trees,
                      "| min node size:", min.node.size,
                      "| error:", round(rf_all$prediction.error,3), collapse = " "))
      }
      # order based on the importance
      tmp <- names(sort(rf_all$variable.importance, decreasing = TRUE))

      # get the gene repetition in the rules list for altogether
      times <- matrix(unlist(strsplit(tmp, "__")), nrow = 2, byrow = FALSE)

      # store the rules as two columns
      sorted_rules[["all"]] <- data.frame(Gene1=times[1,],
                                          Gene2=times[2,],
                                          row.names = tmp,
                                          stringsAsFactors = FALSE)

      # get the repetition of the genes in the rules
      times <- ave(times, times, FUN = seq_along)
      gene_repetition[["all"]] <- data.frame(Gene1=times[1,],
                                             Gene2=times[2,],
                                             row.names = tmp,
                                             stringsAsFactors = FALSE)

      rm(tmp, times)
    }

    # run RF for classes separately (One vs Rest)
    if (run_one_vs_rest) {
      for (cl in groups) {

        rfa_cl <- ranger(y = group_TSP(label = L, my_group = cl),
                         x = t(binary),
                         num.trees = num.trees,
                         importance = importance,
                         write.forest = write.forest,
                         keep.inbag = keep.inbag,
                         min.node.size = min.node.size, ...)

        RF_classifiers[[cl]] <- rfa_cl
        if (verbose) {
          message(paste("RF: class", cl,
                        "| num trees:", num.trees,
                        "| min node size:", min.node.size,
                        "| error:", round(rfa_cl$prediction.error,3), collapse = " "))
        }
        # order based on the importance
        tmp <- names(sort(rfa_cl$variable.importance, decreasing = TRUE))

        # get the gene repetition in the rules list for altogether
        times <- matrix(unlist(strsplit(tmp, "__")), nrow = 2, byrow = FALSE)

        # store the rules as two columns
        sorted_rules[[cl]] <- data.frame(Gene1=times[1,],
                                         Gene2=times[2,],
                                         row.names = tmp,
                                         stringsAsFactors = FALSE)

        # get the repetition of the genes in the rules
        times <- ave(times, times, FUN = seq_along)
        gene_repetition[[cl]] <- data.frame(Gene1=times[1,],
                                            Gene2=times[2,],
                                            row.names = tmp,
                                            stringsAsFactors = FALSE)

        rm(tmp, times)
      }
    }

    object_tmp$RF_scheme$sorted_genes    <- sorted_genes
    object_tmp$RF_scheme$sorted_rules    <- sorted_rules
    object_tmp$RF_scheme$gene_repetition <- gene_repetition
    object_tmp$RF_scheme$RF_classifiers  <- RF_classifiers
    object_tmp$RF_scheme$calls           <- param
    return(object_tmp)
  }

  # platform-wise gene filtering
  if (platform_wise == TRUE) {

    # list to store the sorted rules
    sorted_rules     <- empty_list
    gene_repetition  <- empty_list

    # make list for sorted rules for platform-wise
    plat_rules <- vector("list", length(studies))
    names(plat_rules) <- studies

    # fill each study with list of classes for genes
    for (y in studies) {
      plat_rules[[y]] <- empty_list
    }

    # copy the list to store classifiers
    RF_classifiers <- plat_rules

    # this counter to follow how many platforms for each class
    counter <- empty_list

    # sort the genes in each platform and each class alone
    for(y in studies) {
      plat_samples <- plat_vector == y
      if (verbose) {
        message()
        message(paste("Platform/study: ",y))
      }
      # run RF for all classes together
      if (run_altogether) {

        rf_all <- ranger(x=t(binary[, plat_samples]),
                         y=factor(L[plat_samples]),
                         importance = importance,
                         write.forest = write.forest,
                         keep.inbag = keep.inbag,
                         min.node.size = min.node.size, ...)

        RF_classifiers[[y]][["all"]] <- rf_all
        if (verbose) {
          message(paste("RF: all classes",
                        "| num trees:", num.trees,
                        "| min node size:", min.node.size,
                        "| error:", round(rf_all$prediction.error,3), collapse = " "))
        }
        tmp1 <- sort(rf_all$variable.importance, decreasing = TRUE)
        tmp1 <- names(tmp1)
        plat_rules[[y]][["all"]] <- tmp1
        rm(tmp1)
      }

      #Run RF for classes separately (One vs Rest)
      if (run_one_vs_rest){
        for (cl in groups){

          tmp_check <- as.character(group_TSP(label = L[plat_samples],
                                              my_group = cl))

          if (length(unique(tmp_check)) == 1) {
            if (verbose) {
              message("RF: skip class",cl,"in platform",y)
            }
            plat_rules[[y]] <- plat_rules[[y]][names(plat_rules[[y]]) != cl]
            next
          }

          rfa_cl <- ranger(y = group_TSP(label = L[plat_samples], my_group = cl),
                           x = t(binary[, plat_samples]),
                           num.trees = num.trees,
                           importance = importance,
                           write.forest = write.forest,
                           keep.inbag = keep.inbag,
                           min.node.size=min.node.size, ...)

          RF_classifiers[[y]][[cl]] <- rfa_cl
          counter[[cl]] <- c(counter[[cl]], y)
          if (verbose) {
            message(paste("RF: class", cl,
                          "| num trees:", num.trees,
                          "| min node size:", min.node.size,
                          "| error:", round(rfa_cl$prediction.error,3), collapse = " "))
          }
          tmp1 <- sort(rfa_cl$variable.importance, decreasing = TRUE)
          tmp1 <- names(tmp1)
          plat_rules[[y]][[cl]] <- tmp1
          rm(tmp1)
        }
      }
    }

    # get the top genes in all platforms
    # fill altogether slot
    if (run_altogether){
      tmp <- as.vector(t(rbind(sapply(plat_rules, '[[', "all"))))

      # be sure it mentioned same to platforms number
      tmp2 <- table(tmp) == length(studies)
      tmp2 <- names(tmp2[tmp2 == TRUE])
      tmp  <- tmp[tmp %in% tmp2]

      # last mentioned in the list is last ranked
      tmp2 <- !duplicated(tmp, fromLast = TRUE)
      tmp  <- tmp[tmp2]#[1:rules_altogether]

      # get the gene repetition in the rules list for altogether
      times <- matrix(unlist(strsplit(tmp, "__")), nrow = 2, byrow = FALSE)

      # store the rules as two columns
      sorted_rules[["all"]] <- data.frame(Gene1=times[1,],
                                          Gene2=times[2,],
                                          row.names = tmp,
                                          stringsAsFactors = FALSE)

      # get the repetition of the genes in the rules
      times <- ave(times, times, FUN = seq_along)
      gene_repetition[["all"]] <- data.frame(Gene1=times[1,],
                                             Gene2=times[2,],
                                             row.names = tmp,
                                             stringsAsFactors = FALSE)

      rm(tmp, times)
    } # otherwise it is already NULL

    # fill one vs rest for each class
    if (run_one_vs_rest) {
      for (cl in groups) {
        # get the sorted genes for that class from
        # platforms where this class is available
        tmp <- sapply(plat_rules[counter[[cl]]], '[[', cl)

        # get them in a vector (in the right order first first 2nd 2nd ...)
        tmp <- as.vector(t(rbind(tmp)))

        # be sure it mentioned same to platforms number
        tmp2 <- table(tmp) == length(counter[[cl]])
        tmp2 <- names(tmp2[tmp2 == TRUE])
        tmp  <- tmp[tmp %in% tmp2]

        # last mentioned in the list is last ranked
        tmp2 <- !duplicated(tmp, fromLast = TRUE)

        tmp <- tmp[tmp2]#[1:rules_one_vs_rest]

        # get the gene repetition in the rules list for altogether
        times <- matrix(unlist(strsplit(tmp, "__")), nrow = 2, byrow = FALSE)

        # store the rules as two columns
        sorted_rules[[cl]] <- data.frame(Gene1=times[1,],
                                         Gene2=times[2,],
                                         row.names = tmp,
                                         stringsAsFactors = FALSE)

        # get the repetition of the genes in the rules
        times <- ave(times, times, FUN = seq_along)
        gene_repetition[[cl]] <- data.frame(Gene1=times[1,],
                                            Gene2=times[2,],
                                            row.names = tmp,
                                            stringsAsFactors = FALSE)

        rm(tmp, times)

      }
    } # otherwise it is already NULL for each class

    # fill the object
    object_tmp$RF_scheme$sorted_genes       <- sorted_genes
    object_tmp$RF_scheme$sorted_rules       <- sorted_rules
    object_tmp$RF_scheme$gene_repetition    <- gene_repetition
    object_tmp$RF_scheme$RF_classifiers     <- RF_classifiers
    object_tmp$RF_scheme$calls              <- param
    return(object_tmp)
  }

  # end
}

# function to optimize RF parameters
optimize_RF <- function(data_object,
                        sorted_rules_RF,
                        parameters,
                        overall=c("Accuracy","Kappa",
                                  "AccuracyLower","AccuracyUpper",
                                  "AccuracyNull","AccuracyPValue"
                                  ,"McnemarPValue")[1:2],
                        byclass=c("Sensitivity","Specificity",
                                  "Pos Pred Value","Neg Pred Value",
                                  "Precision","Recall",
                                  "F1","Prevalence",
                                  "Detection Rate","Detection Prevalence",
                                  "Balanced Accuracy")[c(11)],
                        seed = 123456,
                        test_object = NULL,
                        impute = TRUE,
                        impute_reject = 0.67,
                        verbose = FALSE) {

  ### Checks
  if (class(data_object)[1] != "multiclassPairs_object") {
    stop("data_object argument requires multiclassPairs_object!
  Use ReadData function to generate it.")
  }

  if (class(sorted_rules_RF)[1] != "RandomForest_sorted_rules") {
    stop("This function requires RandomForest_sorted_rules object!
  Use filter_rules_RF function to generate it.")
  }

  if (class(test_object)[1] != "multiclassPairs_object" & !is.null(test_object)) {
    stop("test_object requires multiclassPairs_object!
  Use ReadData function to generate it.
  If test_object is NULL then the same training data will be used.
")
  }

  if (!is.data.frame(parameters)) {
    stop("parameters should be a dataframe with column names match argument names in train_RF function!")
  }

  # overall
  if (any(!overall%in%c("Accuracy","Kappa",
                        "AccuracyLower","AccuracyUpper",
                        "AccuracyNull","AccuracyPValue"
                        ,"McnemarPValue"))) {
    stop("overall argument should be a vector with one or more of these variables:
         c('Accuracy','Kappa',
         'AccuracyLower','AccuracyUpper','AccuracyNull',
         'AccuracyPValue','McnemarPValue')")
  }

  # by class
  if (any(!byclass%in%c("Sensitivity","Specificity",
                        "Pos Pred Value","Neg Pred Value",
                        "Precision","Recall",
                        "F1","Prevalence",
                        "Detection Rate","Detection Prevalence",
                        "Balanced Accuracy"))) {
    stop("byclass argument should be a vector with one or more of these variables:
         c('Sensitivity','Specificity',
           'Pos Pred Value','Neg Pred Value',
           'Precision','Recall',
           'F1','Prevalence',
           'Detection Rate','Detection Prevalence',
           'Balanced Accuracy')")
  }

  # get the classes
  groups <- unique(data_object$data$Labels)

  #### prepare output summary df
  # calculate how many columns and row we need for the results df
  # nrow is same  of parameters
  nr <- nrow(parameters)

  # row names
  r_nam <- paste0("Trial_",1:nr)

  # columns should be parameters then n_gene, n_rules then overall things then by class
  # get the col names
  c_nam <- colnames(parameters)
  to_add <- c("n_genes","n_rules", overall,
              as.vector(outer(groups, byclass, paste, sep=".")))
  c_nam <- c(c_nam, to_add)
  nc <- length(c_nam)

  # create the summary df
  out_df <- data.frame(matrix(NA,
                              nrow = nr,
                              ncol = nc,
                              dimnames = list(c(r_nam),
                                              c(c_nam))),
                       check.names = FALSE)

  out_df[,colnames(parameters)] <- parameters


  ##### prepare the output object
  # prepare list for confusion matrices and df
  res_list <- list(summary=NULL,
                   confusionMatrix = vector("list", nr),
                   errors = vector("list", nr),
                   calls = match.call())
  names(res_list$confusionMatrix) <- r_nam
  names(res_list$errors) <- r_nam

  class(res_list) <- "optimize_RF_output"

  # print message about the seed
  if(verbose){
    if (is.null(seed)) {
      message("seed is NULL! It is recommended to use a seed to have reproducible results!")
    } else {
      message("Used seed is ",seed)
    }
  }

  # for loop for the trials
  for (i in 1:nrow(parameters)) {
    print(paste("Trial:",i))

    # get arguments
    args1 <- as.list(parameters[i,])
    args2 <- list(data_object = data_object,
                  sorted_rules_RF = sorted_rules_RF)
    args  <- c(args2,args1)
    args["verbose"] <- verbose
    args["seed"]    <- seed

    # train
    RF_classifier <- try(do.call(train_RF, args))

    if (any(class(RF_classifier)=="try-error") |
        any(class(RF_classifier)!="rule_based_RandomForest")) {
      # some code to store that there is an error
      out_df[i,to_add] <- "error"
      res_list$errors[[i]] <- RF_classifier
      message("Error produced by this trial... skip trial number ",i)
      next()
    }

    # get the number of rules and genes in the model
    out_df[i,"n_genes"] <- length(RF_classifier$RF_scheme$genes)
    out_df[i,"n_rules"] <- nrow(RF_classifier$RF_scheme$rules)

    # predict on test data
    # if the user input a different test data
    if (!is.null(test_object)) {

      pred <- predict_RF(classifier = RF_classifier,
                         Data = test_object,
                         impute = TRUE,
                         verbose = verbose,
                         impute_reject = impute_reject)

      pred    <- pred$predictions

      # in case some samples were rejected due to the imputation
      # get sample names
      if (is.factor(pred)) {
        wanted_sam <- names(pred)
      }
      if (is.matrix(pred)){
        wanted_sam <- rownames(pred)
      }

      # get the ref labels for these samples
      wanted_sam <- order(match(colnames(test_object$data$Data),
                                wanted_sam))[1:length(wanted_sam)]
      ref_lab <- test_object$data$Labels[wanted_sam]

    } else {

      # get the training predictions
      pred    <- RF_classifier$RF_scheme$RF_classifier$predictions
      ref_lab <- data_object$data$Labels
    }

    # get the prediction labels
    # if the classifier trained using probability	= FALSE
    if (is.factor(pred)) {
      pred <- as.character(pred)
    }

    # if the classifier trained using probability	= TRUE
    if (is.matrix(pred)) {
      pred <- colnames(pred)[max.col(pred)]
    }

    # produce the confusion matrix by Caret package
    con <- confusionMatrix(data =factor(pred,
                                        levels = groups),
                           reference = factor(ref_lab,
                                              levels = groups),
                           mode = "everything")

    res_list$confusionMatrix[[i]] <- con

    for (o in overall) {
      out_df[i,o] <- con$overall[[o]]
    }

    for (b in byclass) {
      for (cl in groups) {
        out_df[i,paste(cl,b, sep = ".")] <- con$byClass[paste("Class:",cl),b]
      }
    }

    rm(RF_classifier)
  }

  res_list$summary <- out_df
  res_list$summary$seed <- seed
  return(res_list)
}


# train RF classifier
train_RF <- function (data_object,
                      sorted_rules_RF,

                      gene_repetition = 1,
                      rules_altogether = 200,
                      rules_one_vs_rest = 200,

                      run_boruta = FALSE,
                      plot_boruta = FALSE,
                      boruta_args = list(doTrace = 1),

                      num.trees = 500,
                      min.node.size = 1,
                      importance = "impurity",
                      write.forest = TRUE,
                      keep.inbag = TRUE,
                      probability = TRUE,
                      verbose = TRUE, ...) {

  ### Checks
  if (class(data_object)[1] != "multiclassPairs_object") {
    stop("This function requires multiclassPairs_object!
              Use ReadData function to generate it.")
  }

  if (class(sorted_rules_RF)[1] != "RandomForest_sorted_rules") {
    stop("This function requires RandomForest_sorted_rules object!
              Use filter_rules_RF function to generate it.")
  }

  if (write.forest != TRUE) {
    stop("No model will be saved! Please set write.forest to TRUE")
  }

  if (hasArg(rules_altogether)) {
    if (rules_altogether < 0 |
        !is.numeric(rules_altogether) |
        length(rules_altogether) != 1){
      stop("rules_altogether should be positive number!")
    }
  }

  if (hasArg(rules_one_vs_rest)) {
    if (rules_one_vs_rest < 0 |
        !is.numeric(rules_one_vs_rest) |
        length(rules_one_vs_rest) != 1){
      stop("rules_one_vs_rest should be positive number!")
    }
  }

  if (!importance[1] %in%  c('impurity', 'impurity_corrected', 'permutation') |
      length(importance) > 1 |
      !is.character(importance)) {
    stop("importance variable should be one of these options 'impurity', 'impurity_corrected',  or 'permutation'!")
  }

  ### extract data and labels
  # get the data
  D <- data_object$data$Data

  # get the labels vector
  L <- data_object$data$Labels

  # get the classes
  groups <- unique(L)

  # get the call
  param <- match.call()

  ### additional checks
  # Remove genes with NAs
  if (sum(!complete.cases(D)) > 0) {
    message("These genes will be excluded from the analysis due to NAs:")
    message(capture.output(cat(rownames(D)[!complete.cases(D)])))
    D <- D[complete.cases(D),]
  }

  # Warning if wanted sorted rules > than the available rules - for altogether
  if (hasArg(rules_altogether)) {
    if (rules_altogether > length(sorted_rules_RF[[1]]$sorted_rules$all)) {
      message("NOTE!")
      message("rules_altogether > number of rules in sorted rules object")
      message("This means all rules in 'all' slot will be used!")
      message()
    }
  }
  # else {
  #   rules_altogether <- length(sorted_rules_RF[[1]]$sorted_rules$all)
  # }

  # Warning if wanted sorted rules > than the available rules - for classes
  if (hasArg(rules_one_vs_rest)) {
    tmp <- sapply(sorted_rules_RF[[1]]$sorted_rules[groups],length)
    if (any(rules_one_vs_rest > tmp)) {
      message("NOTE!")
      message("rules_one_vs_rest > number of rules in sorted rules object")
      message("This means all rules will be used for these classes:")
      message(capture.output(cat(names(tmp)[rules_one_vs_rest > tmp])))
      message()
    }
    rm(tmp)
  }
  # else {
  #   rules_one_vs_rest <- max(sapply(sorted_rules_RF[[1]]$sorted_rules[groups],length))
  # }

  # check the distribution of classes in the data
  if (verbose) {
    message("Classes - table:")
    print(table(Classes  = L), right = FALSE)
    message()
  }

  ### Object
  # create empty Random Forests sorted genes object
  if (verbose) {
    message("Creating new Random Forest object")
    message()
  }
  object_tmp <- list(
    RF_scheme = list(genes=NULL,
                     rules=NULL,
                     mode=NULL,
                     boruta=NULL,
                     RF_classifier=NULL,
                     calls=c()))
  class(object_tmp) <- "rule_based_RandomForest"


  ###
  # get the wanted genes
  if (verbose) {
    message("Extracting the needed rules from the sorted rules object...")
  }
  # subset the sorted genes based on the needed number of genes
  rules <- data.frame(V1=character(),
                      V2=character(),
                      stringsAsFactors=FALSE)

  # from altogether
  if (rules_altogether > 0) {

    if (is.null(sorted_rules_RF[[1]]$sorted_rules$all)) {
      if (verbose) {
        message("No rules are available form all classes slot")
      }
    } else {

      # from all
      tmp <- sorted_rules_RF[[1]]$sorted_rules$all
      #tmp <- matrix(unlist(strsplit(tmp, "__")), nrow = 2, byrow = FALSE)

      # get the repetition of the genes in the rules
      #times <- ave(tmp, tmp, FUN = seq_along)
      times <- sorted_rules_RF[[1]]$gene_repetition$all

      # keep only genes repeated specific number of times
      keep <- as.integer(times[,1]) <= gene_repetition &
        as.integer(times[,2]) <= gene_repetition

      tmp <- tmp[keep, , drop=FALSE]

      # let the user know how many rules we collected from all
      if (verbose) {
        message("Altogether rules:")
        message(paste("available rules:",
                      nrow(sorted_rules_RF[[1]]$sorted_rules$all), "rules",
                      collapse = " "))
        message(paste("removing the rules with repeated genes (",
                      gene_repetition,"times allowed ):",
                      nrow(tmp),"rules left",
                      collapse = " "))
        message(paste("get", min(rules_altogether, nrow(tmp)), "rules",
                      collapse = " "))
        message()
      }
      # get the wanted number of rules
      tmp <- tmp[1:min(rules_altogether, nrow(tmp)),]

      # get
      rules <- rbind(rules, tmp)
      rm(tmp)
    }
  }

  # from one_vs_rest
  if (rules_one_vs_rest > 0) {
    for (cl in groups) {

      if (is.null(sorted_rules_RF[[1]]$sorted_rules[[cl]])) {
        if (verbose) {
          message("No rules are available form class: ",cl)
          message()
        }
        next
      }

      tmp <- sorted_rules_RF[[1]]$sorted_rules[[cl]]
      #tmp <- matrix(unlist(strsplit(tmp, "__")), nrow = 2, byrow = FALSE)

      # get the repetition of the genes in the rules
      #times <- ave(tmp, tmp, FUN = seq_along)
      times <- sorted_rules_RF[[1]]$gene_repetition[[cl]]

      # keep only genes repeated specific number of times
      keep <- as.integer(times[,1]) <= gene_repetition &
        as.integer(times[,2]) <= gene_repetition

      tmp <- tmp[keep, , drop=FALSE]

      # let the user know how many rules we collected from all
      if (verbose) {
        message(paste(cl,"class rules:",
                      collapse = " "))
        message(paste("available rules:",
                      nrow(sorted_rules_RF[[1]]$sorted_rules[[cl]]), "rules",
                      collapse = " "))
        message(paste("removing the rules with repeated genes (",
                      gene_repetition,"times allowed )...",
                      nrow(tmp)," rules left",
                      collapse = " "))
        message(paste("get", min(rules_one_vs_rest, nrow(tmp)), "rules",
                      collapse = " "))
        message()
      }

      # get the wanted number of rules
      tmp <- tmp[1:min(rules_one_vs_rest, nrow(tmp)),]

      # get
      rules <- rbind(rules, tmp)
    }
  }

  # get the unique rules from the pooled rules
  # thanks for https://stackoverflow.com/a/25298863
  rules <- rules[!duplicated(data.frame(list(do.call(pmin,rules),
                                             do.call(pmax,rules)))),]

  gene1 <- rules[,1]
  gene2 <- rules[,2]

  for_flip <- rowSums(D[gene1,] < D[gene2,])
  for_flip <- for_flip > (ncol(D)/2)

  from1 <- gene1[for_flip]
  from2 <- gene2[for_flip]

  gene1[for_flip] <- from2
  gene2[for_flip] <- from1

  # store the genes and rules
  genes <- unique(c(gene1,gene2))
  rules <- cbind(gene1,
                 gene2,
                 rule=paste0(gene1,"<",gene2))

  # give how many unique rules we have after merging
  if (verbose) {
    message("Pooling rules ...")
    message(paste("There are",
                  nrow(rules),
                  "unique rules after pooling",
                  collapse = " "))
  }

  # get binary matrix for training RF
  binary    <- D[gene1,] < D[gene2,]
  rownames(binary) <- paste0(gene1,"__",gene2)

  ###
  # run burota
  if (run_boruta) {
    if (verbose) {
      message("Run Boruta...")
    }

    # doTrace argument should be 1
    # if (verbose) {
    #   message("doTrace argument for Boruta turned to 1 by the tool!")
    # }
    # boruta_args[["doTrace"]] <- 1

    # seed
    if (!is.null(list(...)[["seed"]])) {
      set.seed(list(...)[["seed"]])
    }

    args_bor <- list(x = t(binary),
                     y =  factor(L, levels = unique(L)))
    args_bor  <- c(args_bor, boruta_args)
    out <- do.call(Boruta, args_bor)

    # let the user know how many rules left
    if (verbose) {
      message(paste(" Reject",sum(out$finalDecision == "Rejected"),"rules",
                    collapse = " "))

      message(paste(" Use ",sum(out$finalDecision != "Rejected"),
                    "rules for the final RF classifier...",
                    collapse = " "))
    }

    if (plot_boruta) {
      if(is.null(out$ImpHistory)){
        if (verbose) {
          message("holdHistory is FALSE, so no plots will be generated!")
        }
      } else {
        # plots
        if (verbose) {
          message("Plot feature importance history by Boruta...")
        }
        print(plot(out, sort=TRUE))
        print(plotImpHistory(out))
      }
    }

    if (sum(out$finalDecision != "Rejected") == 0) {
      stop("No rules left after burota")
    }

    # get good rules
    good_rules <- names(out$finalDecision[out$finalDecision != "Rejected"])

    # subset for good rules
    binary     <- binary[good_rules,]

    # save
    object_tmp$RF_scheme$boruta <- out

    # store the genes and rules
    tmp <- matrix(unlist(strsplit(good_rules, "__")),
                  ncol = 2, byrow = TRUE)

    genes <- unique(c(tmp[,1], tmp[,2]))

    rules <- cbind(gene1=tmp[,1], gene2=tmp[,2],
                   rule=paste0(tmp[,1],"<",tmp[,2]))
  }

  # train the model
  # Note: any seed will be used here also
  rf_all <- ranger(x = t(binary),
                   y = factor(L),
                   num.trees = num.trees,
                   importance = importance,
                   write.forest = write.forest,
                   keep.inbag = keep.inbag,
                   min.node.size = min.node.size,
                   probability = probability, ...)
  if (verbose) {
    message(paste("RF: Done!",
                  "| total genes:", length(genes),
                  "| total rules:", nrow(rules),
                  "| num trees:", num.trees,
                  "| min node size:", min.node.size,
                  "| error:", round(rf_all$prediction.error,3), collapse = " "))
  }

  # store the "centroids" to be used to impute missing values in testing data
  #Get the most common value function
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }

  #Get the most common value for each class (groups) Get it from the TRAINING_LABELS vector
  ModeVector <- lapply(groups,function(x){
    apply(binary[, which(L==x)], 1, getmode)
  })
  names(ModeVector) <- groups
  ModeVector <- do.call("cbind", ModeVector)


  # store in object
  object_tmp$RF_scheme$genes <- genes
  object_tmp$RF_scheme$rules <- data.frame(rules,
                                           stringsAsFactors = FALSE)
  object_tmp$RF_scheme$mode <- ModeVector
  object_tmp$RF_scheme$RF_classifier <- rf_all
  object_tmp$RF_scheme$calls <- param

  return(object_tmp)
}

# predict function for RF
predict_RF <- function(classifier,
                       Data,
                       impute = FALSE,
                       impute_reject=0.67,
                       verbose = TRUE) {

  # check the object class
  if (!class(Data)[1] %in% c("multiclassPairs_object",
                             "ExpressionSet",
                             "data.frame",
                             "matrix")) {
    stop("Data should be class:
    matrix/data.frame/ExpressionSet/multiclassPairs_object from ReadData function!")
  }

  # check classifier object
  if (class(classifier)[1] != "rule_based_RandomForest") {
    stop("classifier should be rule_based_RandomForest object from train_RF function!")
  }

  if (!is.numeric(impute_reject) |
      !length(impute_reject) == 1 |
      any(impute_reject >= 1) |
      any(impute_reject <= 0)) {
    stop("impute_reject argument should be a number between 0 and 1!")
  }

  # get the data matrix
  if (is.data.frame(Data)) {
    D <- Data
  }

  if (is.matrix(Data)) {
    D <- as.data.frame(Data, stringsAsFactors = FALSE)
  }

  # if the input Data is ExpressionSet object
  if (class(Data)[1] == "ExpressionSet") {
    # extract the expression matrix from the ExpressionSet
    D <- as.data.frame(exprs(Data), stringsAsFactors = FALSE)
  }

  if (class(Data)[1]  ==  "multiclassPairs_object") {
    D <- Data$data$Data
  }

  # extract the genes and the rules
  genes <- classifier$RF_scheme$genes
  rules <- classifier$RF_scheme$rules

  # check if all genes are in the data
  if (any(!genes %in% rownames(D))) {

    if (verbose){
      message("These genes are not found in the data:")
      message(capture.output(cat(genes[!genes %in% rownames(D)])))
      message("Gene names should as rownames and sample names as columns!")
      message("Check the genes in classifier object to see all the needed genes.")
    }

    if (impute == FALSE) {
      stop("All genes should be in the data with no NA values! Or you can turn impute argument to TRUE to impute missed genes to the closest class for each sample!")
    }

    if (impute == TRUE & verbose) {
      message("Missed genes will be imputed to the closest class for each sample!")
    }
  }

  # create empty matrix for the data
  complete <- data.frame(matrix(data = NA,
                                nrow = length(genes),
                                ncol = ncol(D),
                                dimnames = list(genes, colnames(D))),
                         check.names = FALSE,
                         stringsAsFactors = FALSE)

  # fill it with data
  found <- genes[genes %in% rownames(D)]
  complete[found,colnames(D)] <- D[found,]

  # Remove genes with NAs
  if (sum(!complete.cases(complete)) > 0) {
    if (verbose){
      message("These genes have NAs:")
      message(paste(rownames(complete)[!complete.cases(complete)],
                    collapse = " "))
    }

    if (impute == FALSE) {
      message("Turn impute to TRUE to impute NAs to the closest class for each sample with NAs!")
      stop("Gene which is used in the classifier should not have NAs!")
    }

    if (impute == TRUE & verbose) {
      message("These genes will be imputed to the closest class for each sample with NAs")
    }
  }

  # produce the binary matrix
  binary <- complete[rules$gene1,] < complete[rules$gene2,]
  rownames(binary) <- paste0(rules$gene1,
                             "__",
                             rules$gene2)


  #Impute if needed
  if (impute) {
    # get the mode values from the training data - stored in the classifier object
    mode_df <- classifier$RF_scheme$mode

    # to store the index for the samples to be removed
    # due to lack of a lot of rules
    to_remove_sam <- c()

    for(i in 1:ncol(binary)){

      # get which rules are missed in this sample
      is_na <- is.na(binary[,i])

      # if everything is OK then go to the next sample
      if (sum(is_na) == 0) {
        next
      }

      # skip the sample if it misses >0.67 of rules
      if (sum(is_na) > (nrow(binary)*impute_reject)) {
        to_remove_sam <- c(to_remove_sam,i)
        next()
      } else {
        # give warning if the sample misses >0.5 of the rules
        if (sum(is_na) > (nrow(binary)*0.5) & verbose) {
          message("More than the half of the rules need imputation for this sample:")
          message(colnames(binary)[i])
          message("This could affect the prediction accuracy for this sample!")
        }
      }

      # remove the NAs before find the dist
      ok_rules <- names(which(is_na==FALSE))
      sam      <- binary[ok_rules,i, drop=FALSE]

      dist_mat <- as.matrix(dist(t(cbind(sam, mode_df[ok_rules,])),
                                 method="binary"))

      # remove the first because it is the sample itself
      closest  <- which.min(dist_mat[-1,1])
      closest  <- names(closest)[1]

      # get the rules those need imputation for this sample
      impute_rules <- names(which(is_na==TRUE))

      # get the mode values as imputations
      binary[impute_rules, i] <- mode_df[impute_rules, closest]
    }

    # tell the user that we skipped these samples
    if (length(to_remove_sam)>0) {
      message("#####")
      message("More than two thirds of the rules are missed in ",
              length(to_remove_sam),
              " sample(s), because of that these sample(s) were removed from the prediction:")
      message(paste0(colnames(binary)[to_remove_sam],
                     collapse = " "))
      message("#####")

      binary  <- binary[,-to_remove_sam, drop=FALSE]
    }

    if (any(dim(binary)== 0)) {
      stop("No samples left!")
    }
  }

  # predict by original ranger function
  results <- predict(classifier[[1]]$RF_classifier,
                     data = t(binary))

  # give the prediction the sample names
  if (is.matrix(results$predictions)) {
    rownames(results$predictions) <- colnames(binary)

    # get the highest score
    pred <- as.data.frame(results$predictions, stringsAsFactors = FALSE)

    # get the prediction labels
    results$predictions_classes <- colnames(pred)[max.col(pred,
                                                          ties.method = "first")]

    names(results$predictions_classes) <- rownames(results$predictions)

    # to generate warnings if there is ties
    first <- colnames(pred)[max.col(pred,
                                    ties.method = "first")]
    last  <- colnames(pred)[max.col(pred,
                                    ties.method = "last")]
    if (sum(first != last)>0) {
      message(paste("Score ties were found in", sum(first != last),
                    "out of",nrow(pred),"samples in the data",
                    collapse = " "))

    }
  }

  if (is.factor(results$predictions)) {
    names(results$predictions) <- colnames(binary)
  }
  #
  return(results)
}


# plot binary for classifier based on binary RF scheme
plot_binary_RF <- function(Data,
                           classifier,
                           ref = NULL,
                           prediction = NULL,
                           as_training = FALSE,
                           platform = NULL,
                           classes = NULL,
                           platforms_ord = NULL,
                           top_anno = c("ref", "prediction", "platform")[1],
                           title = "",
                           binary_col = c("white", "black", "gray"),
                           ref_col = NULL,
                           pred_col = NULL,
                           platform_col = NULL,
                           show_ref = TRUE,
                           show_predictions = TRUE,
                           show_platform = TRUE,
                           show_scores = TRUE,
                           show_rule_name = TRUE,
                           legend = TRUE,
                           cluster_cols = TRUE,
                           cluster_rows = TRUE,
                           anno_height = 0.03,
                           score_height = 0.03,
                           margin = c(0, 5, 0, 5)) {


  ### get classifier ###
  # check classifier object
  if (class(classifier)[1] != "rule_based_RandomForest") {
    stop("classifier should be rule_based_RandomForest object from train_RF function!")
  } else {
    C <- classifier
  }


  ### get data ###
  # check the object class
  if (!class(Data)[1] %in% c("multiclassPairs_object", "ExpressionSet",
                             "data.frame", "matrix")) {
    stop("Data should be class:
  matrix/data.frame/ExpressionSet/multiclassPairs_object from ReadData function!")
  }

  # get the data matrix
  if (is.data.frame(Data)) {
    D <- Data
  }

  if (is.matrix(Data)) {
    D <- as.data.frame(Data, stringsAsFactors = FALSE)
  }

  # if the input Data is ExpressionSet object
  if (class(Data)[1] == "ExpressionSet") {
    # extract the expression matrix from the ExpressionSet
    # D <- as.matrix(exprs(Data))
    D <- as.data.frame(exprs(Data), stringsAsFactors = FALSE)
  }

  if (class(Data)[1] == "multiclassPairs_object") {
    D <- Data$data$Data
  }

  # check if rownames is not NULL to avoid error later
  if (is.null(rownames(D))) {
    stop("Provide feature/gene names as rownames in the Data matrix!")
  }

  ### get classes ###
  if (is.matrix(classifier$RF_scheme$RF_classifier$predictions)) {
    tmp_n <- colnames(classifier$RF_scheme$RF_classifier$predictions)
  }
  if (is.factor(classifier$RF_scheme$RF_classifier$predictions)) {
    tmp_n <- levels(classifier$RF_scheme$RF_classifier$predictions)
  }

  if (!is.null(classes)) {
    # check if all classes are in the classifier object
    if (any(!classes %in% tmp_n)) {
      message("These classes are not found in the classifier object:")
      message(paste0(classes[!classes %in% tmp_n],
                     collapse = " "))
      stop("classes names in classes argument should be similar to the names of the classifiers in classifier object!")
    }
  } else {
    # get the classes based on the names in the classifier object
    classes <- tmp_n
  }

  if (any(!tmp_n %in% classes)) {
    message("Because the classes argument miss these classes then these classes will be removed from the heatmap:")
    message(paste0(tmp_n[!tmp_n %in% classes], collapse = " "))
  }

  ### get ref labels ###
  # if the data is object
  if (class(Data)[1] == "multiclassPairs_object") {
    # get the ref from the object
    if (is.null(ref)) {
      L <- Data$data$Labels
    }
  }

  # get the input ref Labels from the user
  if ((is.character(ref) | is.factor(ref)) & class(Data)[1] !=
      "ExpressionSet") {
    L <- as.character(ref)
  }

  # if the input Data is ExpressionSet object
  if (class(Data)[1] == "ExpressionSet") {

    # extract the Labels - in case it is stored in the
    # ExpressionSet
    if (is.character(ref) & length(ref) == 1) {
      if (ref %in% varLabels(Data)) {
        L <- as.character(pData(Data)[, ref])
      } else {
        message(capture.output(cat("Phenotype data has these variables:",
                                   varLabels(Data), fill = TRUE)))
        stop("Ref label variable is not found in the phenotype data of your ExpressionSet")
      }
    }

    # get the input Labels vector as it is
    if ((is.character(ref) | is.factor(ref)) & length(ref) !=
        1) {
      L <- as.character(ref)
    }
  }

  # no ref labels if the user did not input ref and the input
  # is not multiclassPairs_object
  if (is.null(ref) & class(Data)[1] != "multiclassPairs_object") {
    L <- NULL
  }

  # check the length of the ref labels
  if (length(L) != ncol(D) & !is.null(ref)) {
    message("Number of samples: ", ncol(D))
    message("Labels length: ", length(L))
    stop("Labels vector length are not equal to
       samples in data")
  }

  ### get Platform labels ###
  # if the data is object
  if (class(Data)[1] == "multiclassPairs_object") {
    # get the platform from the object
    if (is.null(platform)) {
      P <- Data$data$Platform
    }
  }

  # get the input platform Labels from the user
  if ((is.character(platform) | is.factor(platform)) & class(Data)[1] !=
      "ExpressionSet") {
    P <- as.character(platform)
  }

  # if the input Data is ExpressionSet object
  if (class(Data)[1] == "ExpressionSet") {

    # extract the platform label - in case it is stored in the
    # ExpressionSet
    if (is.character(platform) & length(platform) == 1) {
      if (platform %in% varLabels(Data)) {
        P <- as.character(pData(Data)[, platform])
      } else {
        message(capture.output(cat("Phenotype data has these variables:",
                                   varLabels(Data), fill = TRUE)))
        stop("Platform/study label variable is not found in the phenotype data of your ExpressionSet")
      }
    }

    # get the input Labels vector as it is
    if ((is.character(platform) | is.factor(platform)) &
        length(platform) != 1) {
      P <- as.character(platform)
    }
  }

  # no platform labels if the user did not input platform and
  # the input is not multiclassPairs_object
  if (is.null(platform) & class(Data)[1] != "multiclassPairs_object") {
    P <- NULL
  }

  # check the length of the platform labels
  if (length(P) != ncol(D) & !is.null(platform)) {
    message("Number of samples: ", ncol(D))
    message("Labels length: ", length(P))
    stop("Platform labels vector length are not equal to
       samples in data")
  }

  ### get platforms_ord ###
  if (!is.null(platforms_ord) & !is.null(P)) {
    # check if all platforms are in platforms_ord
    if (any(!platforms_ord %in% P)) {
      message("These platform/study in platforms_ord are not in found the platform labels:")
      message(platforms_ord[!platforms_ord %in% P])
      stop("platforms_ord argument should have similar names of the platforms/studies in the data!")
    }
  } else {
    # get the platforms_ord based on the names in the classifier
    # object
    platforms_ord <- unique(P)
  }


  ### get prediction ###
  pred <- NULL

  # if as_training is true then extract the prediction labels from the classifier
  if (show_predictions & as_training) {

    tmp_here <- classifier$RF_scheme$RF_classifier$predictions

    if (!is.null(prediction)) {
      message("prediction object will be ignored because as_training is TRUE!")
      message("prediction will be extracted from the classifier object for the training data!")
    }

    if (is.matrix(tmp_here)) {
      pred <- as.data.frame(tmp_here,
                            stringsAsFactors = FALSE)

      # get the prediction labels
      pred$max_score <- colnames(pred)[max.col(pred, ties.method = "first")]
      prediction <- pred
    }

    if (is.factor(tmp_here)) {
      pred <- data.frame(matrix(data = 0,
                                nrow = length(tmp_here),
                                ncol = length(levels(tmp_here))+1,
                                dimnames = list(names(tmp_here),
                                                c(levels(tmp_here),
                                                  "max_score"))
      ))
      pred$max_score = as.character(tmp_here)
      prediction <- pred

      # just to make sure not plot scores when there is no scores
      # in input prediction object
      show_scores <- FALSE
      message("show_scores turned FALSE because no scores are provided in the prediction object!")
      message("You need to train your RF model with probability = TRUE to get scores when you predict classes!")
    }
    rm(tmp_here)

  } else {
    # check if the prediction df is from the prediction function
    if (!is.null(prediction) & any(class(prediction) %in% "ranger.prediction")) {

      if (is.matrix(prediction$predictions)) {
        pred <- as.data.frame(prediction$predictions, stringsAsFactors = FALSE)

        # get the prediction labels
        pred$max_score <- colnames(pred)[max.col(pred, ties.method = "first")]
        prediction <- pred
      }

      if (is.factor(prediction$predictions)) {
        pred <- data.frame(matrix(data = 0,
                                  nrow = length(prediction$predictions),
                                  ncol = length(levels(prediction$predictions))+1,
                                  dimnames = list(names(prediction$predictions),
                                                  c(levels(prediction$predictions),
                                                    "max_score"))
        ))
        pred$max_score = as.character(prediction$predictions)
        prediction <- pred

        # just to make sure not plot scores when there is no scores
        # in input prediction object
        show_scores <- FALSE
        message("show_scores turned FALSE because no scores are provided in the prediction object!")
        message("You need to train your RF model with probability = TRUE to get scores when you predict classes!")
      }
    }
  }

  # check the length of the platform labels
  if (!is.null(pred)) {
    if (nrow(pred) != ncol(D)) {
      message("Number of samples in the data: ", ncol(D))
      message("Number of samples in the prediction: ", nrow(pred))
      message("This could be due to skipped samples during the imputation step in predict_RF function!")
      stop("Predictions should be for the same data!
     Use predict_RF to generate it for this data!")
    }
  }
  ### checks ###
  if (is.null(pred) & is.null(L) & is.null(P)) {
    stop("No available ref, prediction, or platform labels!
     One of them atleast is needed!")
  }

  ### checks for top_anno ###
  # check if the top_anno labels are available
  if (top_anno == "ref" & is.null(L)) {
    message("top_anno can be one of these three:  'ref', 'prediction', 'platform'")
    stop("top annotation (top_anno) is ref while there is no ref labels available!")
  }
  if (top_anno == "ref" & !show_ref) {
    message("show_ref was turned to TRUE because top_anno is 'ref'!")
  }

  if (top_anno == "prediction" & is.null(pred)) {
    stop("top annotation (top_anno) is prediction while there is no prediction dataframe available!
         Use predict_RF function to generate it or use as_training to extract predictions from the classifier object if the plot is for training data!")
  }

  if (top_anno == "prediction" & !show_predictions) {
    message("show_predictions was turned to TRUE because top_anno is 'prediction'!")
  }

  if (top_anno == "platform" & is.null(P)) {
    stop("top annotation (top_anno) is platform while there is no platform labels available!")
  }
  if (top_anno == "platform" & !show_platform) {
    message("show_platform was turned to TRUE because top_anno is 'platform'!")
  }

  if (any(!top_anno %in% c("ref", "prediction", "platform")) |
      !is.character(top_anno) | length(top_anno) != 1) {
    stop("Top annotation argument should be character with one of these options:
     ref prediction platform")
  }

  ### title ###
  if (!is.character(title) | length(title) != 1) {
    stop("Title argument should be character input!")
  }

  ### binary heatmap colors ###
  if (!is.character(binary_col) | length(binary_col) != 3) {
    stop("binary_col should be character input with length of 3!
     Three colors are needed for rules with false, true and NAs.
     By default it is c('white','black','gray')")
  }

  ### ref anno colors ###
  if (show_ref & !is.null(L) & !is.null(ref_col)) {
    if (!is.character(ref_col) | any(!classes %in% names(ref_col))) {
      stop("ref_col should be named character vector for all classes!")
    }
  }

  ### pred anno colors ###
  if (show_predictions & !is.null(pred) & !is.null(pred_col)) {
    if (!is.character(pred_col) | any(!classes %in% names(pred_col))) {
      stop("pred_col should be named character vector for all classes!")
    }
  }

  ### platform anno colors ###
  if (show_platform & !is.null(P) & !is.null(platform_col)) {
    if (!is.character(platform_col) | any(!P %in% names(platform_col))) {
      stop("platform_col should be named character vector for all platforms/studies!")
    }
  }

  ### colors ###
  # determine the colors groups_col
  # thanks to https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
  xx_colors <- c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
                 "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabed4",
                 "#469990", "#dcbeff", "#9A6324", "#fffac8", "#800000",
                 "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9")
  xx_colors2 <- c("#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9",
                  "#469990", "#dcbeff", "#9A6324", "#fffac8", "#800000",
                  "#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
                  "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabed4")

  if (is.null(ref_col) & !is.null(L)) {
    if (length(classes)<20) {
      ref_col <- xx_colors[1:length(classes)]
    } else {
      ref_col <- sample(xx_colors, size = length(classes), replace = T)
    }
    names(ref_col) <- classes
  }
  if (is.null(pred_col) & !is.null(pred)) {
    if (length(classes)<20) {
      pred_col <- xx_colors[1:length(classes)]
    } else {
      pred_col <- sample(xx_colors, size = length(classes), replace = T)
    }
    names(pred_col) <- classes
  }
  if (is.null(platform_col) & !is.null(P)) {
    # xx_colors2 to give it a bit different colors than the classes
    if (length(platforms_ord)<20) {
      platform_col <- xx_colors2[1:length(platforms_ord)]
    } else {
      platform_col <- sample(xx_colors, size = length(platforms_ord), replace = T)
    }
    names(platform_col) <- platforms_ord
  }

  ### get info for top_anno ###
  # find the samples number to be used in the plotting
  # and get samples' names
  sam_names <- colnames(D)

  # get the labels and groups for the top anno
  if (top_anno == "ref") {
    lab        <- L
    groups     <- classes
    groups_col <- ref_col
  }
  if (top_anno == "prediction") {
    lab        <- pred$max_score
    groups     <- classes
    groups_col <- pred_col
  }
  if (top_anno == "platform") {
    lab        <- P
    groups     <- platforms_ord
    groups_col <- platform_col
  }

  ### anno ord ###
  anno_ord <- c("ref", "prediction", "platform")
  anno_ord <- anno_ord[c(!is.null(L) & show_ref,
                         !is.null(pred) & show_predictions,
                         !is.null(P) & show_platform)]
  anno_ord <- anno_ord[!anno_ord %in% top_anno]

  ### get sample order ###
  # cluster the samples in each group
  if (cluster_cols & (top_anno %in% c("ref", "prediction"))) {
    tmp <- c()
    for(i in groups){
      select_samples <- sam_names[lab==i]

      tmp_r <- C$RF_scheme$rules

      tmp_binary <- D[tmp_r[,1],select_samples] < D[tmp_r[,2],select_samples]
      d   <- dist(t(tmp_binary[,select_samples]), method = "euclidean")
      fit <- hclust(d, method="ward.D2")
      tmp <- c(tmp, fit$labels[fit$order])
    }

    sam_ord <- order(match(sam_names, tmp))[1:length(tmp)]
    rm(tmp)
  }

  # cluster samples when platform is the top anno
  if (cluster_cols & top_anno == "platform") {

    # get the rules for all classifiers
    # tmp_r <- data.frame(matrix(NA, ncol = 2, nrow = 0),
    #                     stringsAsFactors = FALSE)

    tmp_r <- C$RF_scheme$rules

    tmp <- c()
    for(i in groups){
      select_samples <- sam_names[lab==i]

      tmp_binary <- D[tmp_r[,1],select_samples] < D[tmp_r[,2],select_samples]
      d   <- dist(t(tmp_binary[,select_samples]), method = "euclidean")
      fit <- hclust(d, method="ward.D2")
      tmp <- c(tmp, fit$labels[fit$order])
    }

    sam_ord <- order(match(sam_names, tmp))[1:length(tmp)]
    rm(tmp)
  }

  if (!cluster_cols) {
    # this will only group samples without clustering
    # based on the input data
    sam_ord <- order(match(lab, groups))
  }

  # change everything based on the new order
  # if null then will still be null
  D         <- D[,sam_ord]
  L         <- L[sam_ord]
  P         <- P[sam_ord]
  pred      <- pred[sam_ord,]
  lab       <- lab[sam_ord]
  sam_names <- sam_names[sam_ord]

  num_sam   <- ncol(D)

  # this should be after clustering find where the lines should
  # be the lines
  splits <- table(lab)[order(match(names(table(lab)), groups))]

  ### to keep the par settings from the user
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  ### plot top_anno ###
  {
    # Subtype annotation
    AreaStart <- 0.94
    SizeUnit <- anno_height
    Size <- SizeUnit * 1
    AreaEnd <- AreaStart - Size
    par(fig = c(0, 1, AreaEnd, AreaStart), mar = margin,
        mgp = c(3, 0.5, 0), new = FALSE)
    plot(c(0, 1), c(0, 1), type = "n", xaxs = "i", yaxs = "i",
         xlab = "", ylab = "", main = "", xlim = c(0, num_sam),
         ylim = c(0, 1), xaxt = "n", yaxt = "n", bty = "n")

    # headlines
    text_positions <- cumsum(splits)[1]/2
    for (i in 1:(length(cumsum(splits)) - 1)) {
      text_positions <- c(text_positions,
                          ((cumsum(splits)[i + 1] -
                              cumsum(splits)[i])/2 +
                             cumsum(splits)[i]))
    }

    # smaller headlines
    mtext(groups, side = 3, line = 0, outer = FALSE, at = text_positions,
          adj = NA, padj = NA, cex = 0.8, col = groups_col,
          font = NA)

    mtext(title, side = 3, line = -1, outer = TRUE, font = 2)


    # draw the subtypes
    axis(side = 2, at = 0.5,
         labels = c("ref"="Ref. labels",
                    "prediction"="Predictions",
                    "platform"="Platform/Study")[top_anno],
         las = 1, cex.axis = 0.7, tick = 0)
    for (f in groups) {
      for (g in which(lab == f)) {
        rect(g - 1, 0, g, 1, col = groups_col[f], border = NA)
      }
    }

    # the box and the white lines
    box(lwd = 1)
    li <- cumsum(splits)
    abline(v=li, lwd = 1.5, lty=1, col="black")
  }

  ### plot next annos ###
  for (i in anno_ord) {
    {
      # Subtype annotation
      Gap      <- 0.0
      AreaStart<- AreaStart-Size-Gap
      SizeUnit <- anno_height
      Size     <- SizeUnit*1
      AreaEnd  <- AreaStart-Size

      par(fig = c(0, 1, AreaEnd, AreaStart), mar = margin,
          mgp = c(3, 0.5, 0), new = TRUE)
      plot(c(0, 1), c(0, 1), type = "n", xaxs = "i", yaxs = "i",
           xlab = "", ylab = "", main = "", xlim = c(0, num_sam),
           ylim = c(0, 1), xaxt = "n", yaxt = "n", bty = "n")

      # draw the annotation name
      axis(side = 2, at = 0.5,
           labels = c("ref"="Ref. labels",
                      "prediction"="Predictions",
                      "platform"="Platform/Study")[i],
           las = 1, cex.axis = 0.7, tick = 0)

      if (i == "ref") {
        tmp_color <- ref_col
        tmp_lab   <- L
      }
      if (i == "prediction") {
        tmp_color <- pred_col
        tmp_lab   <- pred$max_score
      }
      if (i == "platform") {
        tmp_color <- platform_col
        tmp_lab   <- P
      }

      for (f in unique(tmp_lab)) {
        for (g in which(tmp_lab == f)) {
          rect(g - 1, 0, g, 1, col = tmp_color[f], border = NA)
        }
      }

      # the box and the white lines
      box(lwd = 1)
      li <- cumsum(splits)
      abline(v=li, lwd = 1.5, lty=1, col="black")
    }
  }

  ### from here if the top_anno is platform then groups should be classes
  if (top_anno == "platform") {
    groups <- classes
  }
  ### plot scores ###
  if (show_scores & !is.null(pred)){
    score_matrix <- pred[,groups, drop=FALSE]

    Gap      <- 0.01
    AreaStart<- AreaStart-Size-Gap
    SizeUnit <- score_height
    Size     <- SizeUnit*1
    AreaEnd  <- AreaStart-Size

    for (class in colnames(score_matrix)) {

      par(fig=c(0,1,AreaEnd,AreaStart),mar=margin,
          mgp=c(3,0.5,0),new=TRUE)

      barplot(as.numeric(score_matrix[,class]),
              col=pred_col[class],
              space=F,
              xaxs='i', yaxs='i',xlim=c(0,num_sam),border =NA,
              ylim=c(0,1),xaxt="n",yaxt="n",bty="n",ylab="",xlab="")
      box(lwd=1)
      axis(2, at =0.5,labels=paste("Scores:", class),las=1,cex.axis=0.7,tick=0)
      axis(4, at =c(0.1,0.5,0.9),labels=c("0", "0.5", "1"),
           las=1, cex.axis=0.4, tick = FALSE)

      li <- cumsum(splits)
      abline(v=li, lwd = 1.5, lty=1, col="black")

      ###
      if (class == colnames(score_matrix)[ncol(score_matrix)]) {
        next
      }
      ###
      Gap      <- 0.00
      AreaStart<- AreaStart-Size-Gap
      SizeUnit <- score_height
      Size     <- SizeUnit*1
      AreaEnd  <- AreaStart-Size
      ###
    }
  }

  ### plot binary heatmaps ###
  # to know the height of the heatmap
  Size <- AreaEnd-0.08-0.08

  tmp <- C$RF_scheme$rules

  binary <- D[tmp[,1],] < D[tmp[,2],]
  binary <- binary + 1 # to fit with the indexes for the colors

  # cluster the rules
  if (cluster_rows) {
    d       <- dist(binary, method = "euclidean")
    fit     <- hclust(d, method="ward.D2")
    tmp_ord <- fit$labels[fit$order]
  } else {
    tmp_ord <- 1:nrow(binary)
  }

  binary <- binary[tmp_ord, ]

  ###
  Gap       <- 0.005
  AreaStart <- AreaEnd-Gap
  AreaEnd   <- AreaStart-Size
  ###

  par(fig = c(0, 1, AreaEnd, AreaStart),
      mar = margin, mgp = c(3, 0.5, 0), new=TRUE)

  myplot <- plot(c(0,1),c(0,1), type="n", xaxs='i', yaxs='i',
                 xlab = "", ylab = "", main = "",
                 xlim = c(0, num_sam), ylim = c(0,length(tmp_ord)),
                 xaxt = "n", yaxt = "n", bty = "n")

  if (show_rule_name){
    rule_names <- paste(tmp[,1],tmp[,2], sep = "<")
    axis(4, at =seq(0,(length(tmp_ord)-1),1)+0.5,
         labels=rule_names,las=1,cex.axis=0.5,tick=0)
  }

  for(f in 1:ncol(binary)){
    for(g in 1:nrow(binary)){
      rect(f-1,g,f,g-1,col=binary_col[binary[g,f]],border=F,lwd=0)
    }
  }

  # put the class name + number of rules
  axis(2, at = nrow(binary)/2, labels = paste0("All classes\n",
                                               nrow(binary), "rules"),
       las = 1, cex.axis = 0.7, tick = 0, adj=1)
  # axis(2, at = (nrow(binary)/2)+1, labels = paste(nrow(binary), "rules"),
  #      las = 1, cex.axis = 0.7, tick = 0, adj=1)

  box(lwd=1)

  li <- cumsum(splits)
  abline(v=li[-length(li)], lwd = 3, lty=3, col="red")

  ### plot legends
  if (legend) {
    par(fig = c(0, 1, 0.02, (AreaEnd-0.01)),
        mar = margin, mgp = c(3, 0.5, 0), new=TRUE)
    plot(c(0,1),c(0,1), type="n", xaxs='i', yaxs='i',
         xlab = "", ylab = "", main = "",
         xlim = c(0, num_sam), ylim = c(0, 1),
         xaxt = "n", yaxt = "n", bty = "n")

    if (!is.null(P) & show_platform) {
      legend(x = "topright", title = "Platform",
             ncol = length(platforms_ord), cex = 0.5,
             legend = platforms_ord,
             fill = platform_col)
    }

    if (!is.null(L) & show_ref) {
      title <- "Ref"
      if (!is.null(pred) & show_predictions &
          length(ref_col) == length(pred_col) &
          all(ref_col %in% pred_col)) {
        title <- "Classes"
      }
      legend(x = "topleft", title = title,
             ncol = length(classes), cex = 0.5,
             legend = names(ref_col),
             fill = ref_col)
    }

    if (!is.null(pred) & show_predictions &
        (length(ref_col) != length(pred_col) |
         any(!ref_col %in% pred_col))) {
      legend(x = "top", title = "Predictions",
             ncol = length(unique(pred$max_score)), cex = 0.5,
             legend = names(pred_col),
             fill = pred_col)
    }
  }
}

# coclustering plot for RF training data
cocluster_RF <- function(object,
                         classifier,
                         title = "",
                         top_anno = c("ref","platform")[1],
                         classes = NULL,
                         sam_order = NULL,
                         ref_col = NULL,
                         platform_col = NULL,
                         platforms_ord = NULL,
                         show_platform = TRUE,
                         cluster_cols = FALSE,
                         legend = TRUE,
                         anno_height = 0.03,
                         margin = c(0, 5, 0, 5)){
  ### get classifier ###
  # check classifier object
  if (class(classifier)[1] != "rule_based_RandomForest") {
    stop("classifier should be rule_based_RandomForest object from train_RF function!")
  }

  ### get data ###
  # check the object class
  if (!class(object)[1] == "multiclassPairs_object") {
    stop("Object should be a multiclassPairs_object from ReadData function!")
  }

  C       <- classifier$RF_scheme$RF_classifier
  rules   <- classifier$RF_scheme$rules
  D       <- object$data$Data

  L       <- object$data$Labels
  P       <- object$data$Platform

  ### get classes ###
  tmp_n <- colnames(classifier$RF_scheme$RF_classifier$predictions)
  if (!is.null(classes)) {
    # check if all classes are in the classifier object
    if (any(!classes %in% tmp_n)) {
      message("These classes are not found in the classifier object:")
      message(paste0(classes[!classes %in% tmp_n],
                     collapse = " "))
      stop("classes names in classes argument should be similar to the names of the classifiers in classifier object!")
    }
  } else {
    # get the classes based on the names in the classifier object
    classes <- tmp_n
  }

  #
  if (is.null(C$inbag.counts)) {
    stop("call ranger with keep.inbag = TRUE")
  }

  ### title ###
  if (!is.character(title) | length(title) != 1) {
    stop("Title argument should be character input!")
  }

  ### ref anno colors ###
  if (!is.null(ref_col)) {
    if (!is.character(ref_col) | any(!classes %in% names(ref_col))) {
      stop("ref_col should be named character vector for all classes!")
    }
  }

  ### platform anno colors ###
  if (show_platform & !is.null(P) & !is.null(platform_col)) {
    if (!is.character(platform_col) | any(!P %in% names(platform_col))) {
      stop("platform_col should be named character vector for all platforms/studies!")
    }
  }

  ### get platforms_ord ###
  if (!is.null(platforms_ord) & !is.null(P)) {
    # check if all platforms are in platforms_ord
    if (any(!platforms_ord %in% P)) {
      message("These platform/study in platforms_ord are not in found the platform labels:")
      message(platforms_ord[!platforms_ord %in% P])
      stop("platforms_ord argument should have similar names of the platforms/studies in the data!")
    }
  } else {
    # get the platforms_ord based on the names in the classifier
    # object
    platforms_ord <- unique(P)
  }

  ### get the same order from the user ###
  if (!is.null(sam_order)) {
    if (all(!sam_order %in% colnames(D))) {
      stop("sam_order argument should have similar names of the samples in the data!")
    }
  }

  ### get binary matrix ###
  binary <- D[rules$gene1, ] < D[rules$gene2, ]
  rownames(binary) <- paste0(rules$gene1,"__",rules$gene2)
  binary <- as.data.frame(t(binary))

  #
  pred  <- predict(C, binary, type = "terminalNodes")$predictions
  # ntree <- ncol(pred)

  prox  <- data.frame(matrix(NA, nrow(pred), nrow(pred),
                             dimnames = list(colnames(D), colnames(D))),
                      stringsAsFactors = FALSE,
                      check.rows = FALSE,
                      check.names = FALSE)
  n     <- nrow(prox)

  # Get inbag counts
  inbag = simplify2array(C$inbag.counts)

  for (i in 1:n) {
    for (j in 1:n) {
      # Use only trees where both obs are OOB
      tree_idx = inbag[i, ] == 0 & inbag[j, ] == 0
      prox[i, j] = sum(pred[i, tree_idx] == pred[j, tree_idx]) / sum(tree_idx)
    }
  }

  ### colors ###
  # determine the colors groups_col
  # thanks to https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
  xx_colors <- c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
                 "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabed4",
                 "#469990", "#dcbeff", "#9A6324", "#fffac8", "#800000",
                 "#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9")
  xx_colors2 <- c("#aaffc3", "#808000", "#ffd8b1", "#000075", "#a9a9a9",
                  "#469990", "#dcbeff", "#9A6324", "#fffac8", "#800000",
                  "#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
                  "#911eb4", "#42d4f4", "#f032e6", "#bfef45", "#fabed4")

  if (is.null(ref_col)) {
    if (length(classes)<20) {
      ref_col <- xx_colors[1:length(classes)]
    } else {
      ref_col <- sample(xx_colors, size = length(classes), replace = T)
    }
    names(ref_col) <- classes
  }

  if (is.null(platform_col) & !is.null(P)) {
    # xx_colors2 to give it a bit different colors than the classes
    if (length(unique(P))<20) {
      platform_col <- xx_colors2[1:length(unique(P))]
    } else {
      platform_col <- sample(xx_colors, size = length(unique(P)), replace = T)
    }
    names(platform_col) <- unique(P)
  }

  ### get info for top_anno ###
  # find the samples number to be used in the plotting
  # and get samples' names
  sam_names <- colnames(prox)

  # get the labels and groups for the top anno
  if (top_anno == "ref") {
    lab        <- L
    groups     <- classes
    groups_col <- ref_col
  }

  if (top_anno == "platform") {
    lab        <- P
    groups     <- platforms_ord
    groups_col <- platform_col
  }

  ### anno ord ###
  anno_ord <- c("ref", "platform")
  anno_ord <- anno_ord[c(TRUE, !is.null(P) & show_platform)]
  anno_ord <- anno_ord[!anno_ord %in% top_anno]

  ### get sample order ###
  # cluster the samples in each group
  if (cluster_cols & (top_anno %in% c("ref"))) {
    tmp <- c()
    for(i in groups){
      select_samples <- sam_names[lab==i]

      d   <- dist(t(prox[,select_samples]), method = "euclidean")
      fit <- hclust(d, method="ward.D2")
      tmp <- c(tmp, fit$labels[fit$order])
    }

    sam_ord <- order(match(sam_names, tmp))[1:length(tmp)]
    rm(tmp)
  }

  # cluster samples when platform is the top anno
  if (cluster_cols & top_anno == "platform") {

    tmp <- c()
    for(i in groups){
      select_samples <- sam_names[lab==i]

      d   <- dist(t(prox[,select_samples]), method = "euclidean")
      fit <- hclust(d, method="ward.D2")
      tmp <- c(tmp, fit$labels[fit$order])
    }

    sam_ord <- order(match(sam_names, tmp))[1:length(tmp)]
    rm(tmp)
  }

  if (!cluster_cols & is.null(sam_order)) {
    # this will only group samples without clustering
    # based on the input data
    length_sel_groups <- sum(lab %in% groups)
    sam_ord <- order(match(lab, groups))[1:length_sel_groups]
  }

  if (!is.null(sam_order)) {
    # based on the user order
    sam_ord <- order(match(sam_names, sam_order))[1:length(sam_order)]
  }

  # change everything based on the new order
  # if null then will still be null
  D         <- D[,sam_ord]
  L         <- L[sam_ord]
  P         <- P[sam_ord]
  lab       <- lab[sam_ord]
  prox      <- prox[rev(sam_ord),sam_ord]

  sam_names <- sam_names[sam_ord]

  num_sam   <- ncol(D)


  # this should be after clustering find where the lines should
  # be the lines
  splits <- table(lab)[order(match(names(table(lab)), groups))]

  ### to keep the par settings from the user
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  ### plot top_anno ###
  {
    # Subtype annotation
    AreaStart <- 0.94
    SizeUnit  <- anno_height
    Size <- SizeUnit * 1
    AreaEnd <- AreaStart - Size
    par(fig = c(0, 1, AreaEnd, AreaStart), mar = margin,
        mgp = c(3, 0.5, 0), new = FALSE)
    plot(c(0, 1), c(0, 1), type = "n", xaxs = "i", yaxs = "i",
         xlab = "", ylab = "", main = "", xlim = c(0, num_sam),
         ylim = c(0, 1), xaxt = "n", yaxt = "n", bty = "n")

    if (is.null(sam_order)){
      # headlines
      text_positions <- cumsum(splits)[1]/2
      for (i in 1:(length(cumsum(splits)) - 1)) {
        text_positions <- c(text_positions,
                            ((cumsum(splits)[i + 1] -
                                cumsum(splits)[i])/2 +
                               cumsum(splits)[i]))
      }

      # smaller headlines
      mtext(groups, side = 3, line = 0, outer = FALSE, at = text_positions,
            adj = NA, padj = NA, cex = 0.8, col = groups_col,
            font = NA)
    }
    mtext(title, side = 3, line = -1, outer = TRUE, font = 2)


    # draw the subtypes
    axis(side = 2, at = 0.5,
         labels = c("ref"="Ref. labels",
                    "prediction"="Predictions",
                    "platform"="Platform/Study")[top_anno],
         las = 1, cex.axis = 0.7, tick = 0)
    for (f in groups) {
      for (g in which(lab == f)) {
        rect(g - 1, 0, g, 1, col = groups_col[f], border = NA)
      }
    }

    # the box and the white lines
    box(lwd = 1)
    li <- cumsum(splits)
    abline(v=li, lwd = 1.5, lty=1, col="black")
  }

  ### plot next annos ###
  for (i in anno_ord) {
    {
      # Subtype annotation
      Gap      <- 0.0
      AreaStart<- AreaStart-Size-Gap
      SizeUnit <- anno_height
      Size     <- SizeUnit*1
      AreaEnd  <- AreaStart-Size

      par(fig = c(0, 1, AreaEnd, AreaStart), mar = margin,
          mgp = c(3, 0.5, 0), new = TRUE)
      plot(c(0, 1), c(0, 1), type = "n", xaxs = "i", yaxs = "i",
           xlab = "", ylab = "", main = "", xlim = c(0, num_sam),
           ylim = c(0, 1), xaxt = "n", yaxt = "n", bty = "n")

      # draw the annotation name
      axis(side = 2, at = 0.5,
           labels = c("ref"="Ref. labels",
                      "platform"="Platform/Study")[i],
           las = 1, cex.axis = 0.7, tick = 0)

      if (i == "ref") {
        tmp_color <- ref_col
        tmp_lab   <- L
      }
      if (i == "platform") {
        tmp_color <- platform_col
        tmp_lab   <- P
      }

      for (f in unique(tmp_lab)) {
        for (g in which(tmp_lab == f)) {
          rect(g - 1, 0, g, 1, col = tmp_color[f], border = NA)
        }
      }

      # the box and the white lines
      box(lwd = 1)
      li <- cumsum(splits)
      abline(v=li, lwd = 1.5, lty=1, col="black")
    }
  }

  ### plot binary heatmaps ###
  # to know the height of the heatmap
  Size <- AreaEnd-0.08-0.08

  ###
  Gap       <- 0.005
  AreaStart <- AreaEnd-Gap
  AreaEnd   <- AreaStart-Size
  ###

  par(fig = c(0, 1, AreaEnd, AreaStart),
      mar = margin, mgp = c(3, 0.5, 0), new=TRUE)

  myplot <- plot(c(0,1),c(0,1), type="n", xaxs='i', yaxs='i',
                 xlab = "", ylab = "", main = "",
                 xlim = c(0, num_sam), ylim = c(0, num_sam),
                 xaxt = "n", yaxt = "n", bty = "n")

  HM_colors <- colorRampPalette(c("white","darkblue"))(100)

  prox <- round(prox, 2)
  prox <- prox*100

  prox <- as.matrix(prox)
  prox[which(prox<1)]   <- 1
  prox[which(prox>100)] <- 100

  for(f in 1:ncol(prox)){
    for(g in 1:nrow(prox)){
      rect(f-1,g,f,g-1,col=HM_colors[prox[g,f]],border=FALSE,lwd=0)
    }
  }

  box(lwd=1)

  ### plot legends
  if (legend) {
    par(fig = c(0, 1, 0.02, (AreaEnd-0.01)),
        mar = margin, mgp = c(3, 0.5, 0), new=TRUE)
    plot(c(0,1),c(0,1), type="n", xaxs='i', yaxs='i',
         xlab = "", ylab = "", main = "",
         xlim = c(0, num_sam), ylim = c(0, 1),
         xaxt = "n", yaxt = "n", bty = "n")

    if (!is.null(P) & show_platform) {
      legend(x = "topright", title = "Platform/study",
             ncol = length(platforms_ord), cex = 0.5,
             legend = platforms_ord,
             fill = platform_col)
    }

    legend(x = "topleft", title = "Ref labels",
           ncol = length(classes), cex = 0.5,
           legend = names(ref_col),
           fill = ref_col)
  }
}

##### print functions #####
# print object function
print.multiclassPairs_object <- function(x, ...) {


  # print info about the input data and labels
  cat("multiclassPairs object\n")

  cat("*Data:\n")
  cat("   Number of samples:", ncol(x$data$Data),"\n")
  cat("   Number of genes/features:", nrow(x$data$Data),"\n")
  message()

  cat("*Labels and Platforms:\n")
  cat("   Classes:", unique(x$data$Labels),"\n")
  if (!is.null(x$data$Platform)) {
    cat("   Platforms/studies:", unique(x$data$Platform), "\n")
  } else {
    cat("   Platforms/studies: NULL\n")
  }
  message()

  # samples distribution
  cat("*Samples count table:")
  if (is.null(x$data$Platform)) {
    print(table(x$data$Labels))
  } else {
    print(table(x$data$Platform, x$data$Labels))
  }
}

# print function for SB genes objects
print.OnevsrestScheme_genes_TSP <- function(x, ...) {
  # print info about the input data and labels
  cat("sorted genes for One-vs-rest Scheme:\n")

  # print what is not empty in the other slots
  for (y in names(x)) {
    if (all(sapply(x[[y]], is.null))) {
      cat("  Object contains: NULL\n")
    } else {
      cat("  Object contains:\n")
      for (i in names(x[[y]])) {
        if (!is.null(x[[y]][[i]])) {
          cat("     ","-",i,"\n")

          if (i == "filtered_genes") {
            for (z in names(x[[y]][[i]])) {
              cat("        ","- class:",z,":",
                  length(x[[y]][[i]][[z]]),"genes\n")
            }
          }

          if (i == "calls") {
            cat("            ")
            cat(gsub(capture.output(cat(capture.output(x[[y]][[i]]))),
                     pattern = paste0(c(",", ",     "), collapse = "|"),
                     replacement = ",\n            "))
          }

        } else { # if NULL
          cat("     ","-",i,": NULL\n")
        }
        message()
      }
    }
  }
}

# print classifier object - SB
print.OnevsrestScheme_TSP <- function(x, ...) {
  # print info about the input data and labels
  cat("multiclassPairs - One vs rest Scheme\n")

  # print what is not empty in the other slots
  cat("*Classifier:\n")
  for (y in names(x)) {
    if (all(sapply(x[[y]], is.null))) {
      cat("  classifier: NULL\n")
    } else {
      cat("  contains binary classifiers:\n")
      for (i in names(x[[y]])) {
        if (!is.null(x[[y]][[i]])) {
          cat("     - Class:",i,"...",nrow(x[[y]][[i]]$TSPs),"rules\n")
        }
      }
    }
  }

}

# print function objects - RF
print.RandomForest_sorted_genes <- function(x, ...) {
  # print info about the input data and labels
  cat("multiclassPairs - Random Forest scheme\n")

  # print what is not empty in the other slots
  for (y in names(x)) {
    if (all(sapply(x[[y]], is.null))) {
      cat("  Object contains: NULL\n")
    } else {
      cat("  Object contains:\n")
      for (i in names(x[[y]])) {
        if (!is.null(x[[y]][[i]])) {
          cat("     ","-",i,"\n")

          if (i == "sorted_genes") {
            for (z in names(x[[y]][[i]])) {
              cat("        ","- class:",z,":",
                  length(x[[y]][[i]][[z]]),"genes\n")
            }
          }

          if (i == "calls") {
            cat("            ")
            cat(gsub(capture.output(cat(capture.output(x[[y]][[i]]))),
                     pattern = paste0(c(",", ",     "), collapse = "|"),
                     replacement = ",\n            "))
          }

        } else { # if NULL
          cat("     ","-",i,": NULL\n")
        }
        message()
      }
    }
  }
}

# print function for the sorted rules of RF
print.RandomForest_sorted_rules <- function(x, ...) {
  # print info about the input data and labels
  cat("multiclassPairs - Random Forest scheme\n")

  # print what is not empty in the other slots
  for (y in names(x)) {
    if (all(sapply(x[[y]], is.null))) {
      cat("  Object contains: NULL\n")
    } else {
      cat("  Object contains:\n")
      for (i in names(x[[y]])) {
        if (!is.null(x[[y]][[i]])) {
          cat("     ","-",i,"\n")

          if (i == "sorted_genes") {
            for (z in names(x[[y]][[i]])) {
              cat("        ","- class:",z,":",
                  length(x[[y]][[i]][[z]]),"genes\n")
            }
          }

          if (i == "sorted_rules") {
            for (z in names(x[[y]][[i]])) {
              cat("        ","- class:",z,":",
                  length(x[[y]][[i]][[z]])," sorted rules\n")
            }
          }

          if (i == "calls") {
            cat("            ")
            cat(gsub(capture.output(cat(capture.output(x[[y]][[i]]))),
                     pattern = paste0(c(",", ",     "), collapse = "|"),
                     replacement = ",\n            "))
          }

        } else { # if NULL
          cat("     ","-",i,": NULL\n")
        }
        message()
      }
    }
  }
}

# print function for optimize_RF results object
print.optimize_RF_output <- function(x, ...) {
  # print info about the input data and labels
  cat("optimize_RF output\n")
  cat("    Number of trials:", nrow(x$summary),"\n")
  cat("    Number of errors:", sum(!sapply(x$errors, is.null)),"\n")
  cat("
  This object contains:
    - Summary table (access it by using this_object$summary)
    - List of confusion matrices for each successful trial
    - List of errors produced by failed trials\n")
  cat("\n")
  cat("  Generate by:\n")
  cat("    ",gsub(capture.output(cat(capture.output(x$calls))),
                  pattern = paste0(c(",", ",     "), collapse = "|"),
                  replacement = ",\n            "))
  # deparse(substitute(x))
}

# print function for the RF classifier
print.rule_based_RandomForest <- function(x, ...) {
  # print info about the input data and labels
  cat("multiclassPairs - Rule based Random Forest classifier\n")

  # print what is not empty in the other slots
  for (y in names(x)) {
    if (all(sapply(x[[y]], is.null))) {
      cat("  Object contains: NULL\n")
    } else {
      cat("  Object contains:\n")
      for (i in names(x[[y]])) {
        if (!is.null(x[[y]][[i]])) {
          cat("     ","-",i)

          if (i == "genes") {
            cat(":", length(x[[y]][[i]]),"genes\n")
          }

          if (i == "rules") {
            cat(":", nrow(x[[y]][[i]]),"rules\n")
          }

          if (i == "calls") {
            cat(": ")
            cat(gsub(capture.output(cat(capture.output(x[[y]][[i]]))),
                     pattern = paste0(c(",", ",     "), collapse = "|"),
                     replacement = ",\n            "))
          }

          if (i == "RF_classifier"| i== "boruta" | i== "mode") {
            cat("\n")
          }

        } else { # if NULL
          cat("     ","-",i,": NULL\n")
        }
        message()
      }
    }
  }
}
