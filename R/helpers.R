#' Stratify hypotheses based on increasing value of the covariate
#'
#'  Hypotheses are stratified into nbins different strata of (approximately) equal size based on
#' increasing value of the covariate
#'
#' @param covariate Numeric vector of ordinal covariates based on which the stratification will be done.
#' @param nbins Integer, number of groups/strata into which p-values will be split based on covariate.
#' @param ties.method Character specifying how ties are treated, see \code{\link{rank}} function.
#' @param seed Integer, specifies random seed to be used when ties.method=="random".
#' @return A factor with nbins different levels, each entry corresponds to the stratum the i-th hypothesis
#'  was assigned to.
#' @examples
#' covariates <- runif(100)
#' groups <- groups_by_filter(covariates, 10)
#' table(groups)
#' @export
groups_by_filter <- function(covariate, nbins, ties.method = "random", seed = NULL) {
  if (!is.null(seed) && ties.method == "random") {
    # http://stackoverflow.com/questions/14324096/setting-seed-locally-not-globally-in-r?rq=1
    tmp <- runif(1)
    old <- .Random.seed
    on.exit({
      .Random.seed <<- old
    })
    set.seed(as.integer(seed))
  }
  rfs <- rank(covariate, ties.method = ties.method) / length(covariate)
  as.factor(ceiling(rfs * nbins))

  # breaks <- quantile(covariate, probs = seq(0, 1, length.out = nbins+1))
  # cut(covariate, breaks=breaks, include.lowest=TRUE)
}

groups_by_cut <- function(covariate, nbins) {
  nvar <- ncol(covariate)
  nbins_dim <- rep(nbins, nvar)
  #nbins_dim <- rep(1, nvar)
  #i <- 1
  #while (prod(nbins_dim) < nbins) {
  #  pos <- i %% nvar + 1
  #  nbins_dim[pos] <- nbins_dim[pos] + 1
  #  i <- i + 1
  #} # when nvar is too large, it is not possible to consider every variable

  # nbins_dim <- floor(log(nbins)/log(nvar)) #TODO
  # nbins_dim <- max(nbins,2)
  # max1

  groups <- sapply(seq_len(nvar), function(i) {
    nbins_dim_i <- nbins_dim[i]
    covariate_i <- covariate[, i]
    breaks <- quantile(covariate_i, probs = seq(0, 1, length.out = nbins_dim_i + 1)) #interact?
    breaks <- unique(breaks)
    as.integer(cut(covariate_i, breaks = breaks, include.lowest = TRUE))
    ## BurStMisc::ntile
  })

  groups <- groups %>%
    as.data.frame() %>%
    tidyr::unite("group", seq_len(ncol(groups)), sep = "-") %>%
    mutate(group = as.factor(group))

  return(groups)
}


#' @importFrom ranger ranger
group_by_forest_BocaLeek2 <- function(pvalues, covariates, folds, nbins, ntrees, taus = NULL, ntaus =5, maxdepth = NULL, min.node.size = NULL, seed = NULL) {
  m <- length(pvalues)
  # calculate parameters for forest
  if(is.null(maxdepth)){
    maxdepth <- log2(nbins) # TODO should depend on nvar
    maxdepth <- ceiling(maxdepth)
  }
  #maxdepth <- 2 #TODO
  # min.node.size <- 2^(maxdepth+1)-1 #TODO
  if(is.null(min.node.size)) min.node.size <- floor(0.9 * (m / 8))
  
  if(is.null(ntaus)) ntaus <- length(taus) 
  taus_specified <- !is.null(taus)
  mtry <- ceiling(0.9 * ncol(covariates))  # a lot of noise data => high mtry
  
  groups <- matrix(0, nrow = m, ncol = as.numeric(ntrees * ntaus))
  
  for (i in unique(folds)) {
    pvalues_other_folds <- pvalues[folds != i]
    if(!taus_specified){
      quantile_seq <- seq(0, 1, length.out = ntaus+2)[2:(ntaus+1)]
      taus <- quantile(pvalues_other_folds, quantile_seq)
    } 
    for (j in seq_along(taus)) {
      tau <- taus[j]
      # binary indicator from Boca and leek/storey
      data <- data.frame(
        indic = (pvalues >= tau) / (1 - tau),
        covariates,
        folds
      )
      
      data_holdout_fold <- data[folds == i, ]
      data_other_folds <- data[folds != i, ]

      # build forest based on other folds
      forest_other_fold <- ranger::ranger(
        formula = indic ~ . - folds - indic,
        data = data_other_folds,
        num.trees = ntrees,
        mtry = mtry,
        min.node.size = min.node.size, 
        max.depth = maxdepth,
        splitrule = "variance",
        importance = "none",
        #oob.error	= FALSE,  #TODO
        seed = seed
      )

      # predict terminal nodes of holdout_fold
      predict_groups <- predict(forest_other_fold,
        predict.all = TRUE,
        data = data_holdout_fold,
        type = "terminalNodes"
      )
      groups[folds == i, (j - 1) * ntrees + seq_len(ntrees)] <- predict_groups$predictions
    }
  }
  groups <- as.data.frame(groups)
  groups[] <- lapply(groups, as.factor)
  colnames(groups) <- paste0("group", seq_along(groups))
  return(groups)
}


#' @importFrom randomForestSRC rfsrc
group_by_forest_BocaLeek <- function(pvalues, covariates, folds, nbins, ntrees,nsplit =NULL, taus = NULL, ntaus =5, maxdepth = NULL, min.node.size = NULL, seed = NULL) {
  m <- length(pvalues)
  # calculate parameters for forest
  if(is.null(maxdepth)){
    maxdepth <- log2(nbins) # TODO should depend on nvar
    maxdepth <- ceiling(maxdepth)
  }
  
  #maxdepth <- 2 #TODO
  # min.node.size <- 2^(maxdepth+1)-1 #TODO
  if(is.null(min.node.size)) min.node.size <- floor(0.9 * (m / 8))
  if(is.null(nsplit)) nsplit <- 4
  if(is.null(ntaus)) ntaus <- length(taus) 
  taus_specified <- !is.null(taus)
  mtry <- ceiling(0.9 * ncol(covariates))  # a lot of noise data => high mtry
  
  groups <- matrix(0, nrow = m, ncol = as.numeric(ntrees * ntaus))
  
  for (i in unique(folds)) {
    pvalues_other_folds <- pvalues[folds != i]
    if(!taus_specified){
      quantile_seq <- seq(0, 1, length.out = ntaus+2)[2:(ntaus+1)]
      taus <- quantile(pvalues_other_folds, quantile_seq)
    } 
    for (j in seq_along(taus)) {
      tau <- taus[j]
      # binary indicator from Boca and leek/storey
      data <- data.frame(
        indic = (pvalues >= tau) / (1 - tau),
        covariates,
        folds
      )
      
      data_holdout_fold <- data[folds == i, !names(data) %in% c("indic", "folds"), drop = FALSE]
      data_other_folds <- data[folds != i, !names(data) %in% c("folds"), drop = FALSE]
      
      # build forest based on other folds
      forest_other_fold <- randomForestSRC::rfsrc(
        indic ~ . - folds - indic,
        data = data_other_folds,
        ntree = ntrees,
        mtry = mtry,
        nodesize = min.node.size,
        nodedepth = maxdepth,
        splitrule = "mse", #TODO
        nsplit = nsplit, 
        block.size = F,
        forest.wt = F,
        seed = seed)
      
      #useful diagnostics
      if(FALSE){
        var.used <- predict(forest_other_fold, var.used = "by.tree")$var.used 
        var_used_colsum <- colSums(var.used)
        correct_splits <- var_used_colsum[1]/sum(var_used_colsum)
        print(paste("percentage: ",correct_splits))
        
        print(paste("number correct splits:", mean(var.used[,"cov1"])))
      }

      # predict terminal nodes of holdout_fold
      predict_groups <- stats::predict(forest_other_fold, data_holdout_fold, membership = T)
      groups[folds == i, (j - 1) * ntrees + seq_len(ntrees)] <- predict_groups$membership
    }
  }
  #browser()
  groups <- as.data.frame(groups)
  groups[] <- lapply(groups, as.factor)
  colnames(groups) <- paste0("group", seq_along(groups))
  return(groups)
}

group_by_forest_RFCDE <- function(pvalues, covariates, folds, nbins, nbasis, ntrees, basis_system = basis_system, mtry = 3L) {
  m <- length(pvalues)
  covariates <- as.matrix(covariates)
  node_size <- floor(m / (2 * nbins))
  groups <- matrix(0, nrow = nrow(covariates), ncol = ntrees)

  for (i in unique(folds)) {
    holdout_fold <- which(folds == i) # TODO does that make sense? hold out fold?
    pvalues_other_folds <- pvalues[folds != i]
    covariates_other_folds <- covariates[folds != i, ]

    forest_other_folds <- RFCDE(
      x_train = pvalues_other_folds, z_train = covariates_other_folds, n_trees = ntrees, mtry = mtry,
      node_size = node_size, n_basis = nbasis, basis_system = basis_system
    )

    forest_other_folds$rcpp$set_leaves_id()

    for (j in holdout_fold) {
      groups[j, ] <- forest_other_folds$rcpp$traverse_forest(covariates[j, ]) # TODO as matrix
    }
  }

  groups <- as.data.frame(groups)
  groups[] <- lapply(groups, as.factor)
  colnames(groups) <- paste0("group", seq_along(groups))
  return(groups)
}

group_by_forest_Julia <- function(pvalues, covariates, folds, nbins, ntrees, nbasis) {
  # JuliaCall::julia_setup(JULIA_HOME = "/Applications/Julia-1.6.app/Contents/Resources/julia/bin/")
  # JuliaCall::julia_command("cd(\"/Users/default/Google Drive/currentDocumants/Studium/Master/3.Semester/Masterarbeit/Code/IndependentHypothesisWeightingTrees.jl\")")
  # JuliaCall::julia_source("/Users/default/Google Drive/currentDocumants/Studium/Master/3.Semester/Masterarbeit/Code/IndependentHypothesisWeightingTrees.jl/example/wrapper.jl")

  nbins <- as.integer(nbins)
  ntrees <- as.integer(ntrees)
  nbasis <- as.integer(nbasis) # must be > 2!

  m <- length(pvalues)
  covariates <- as.matrix(covariates)
  groups <- matrix(0, nrow = nrow(covariates), ncol = ntrees)

  for (i in unique(folds)) {
    pvalues_other_folds <- pvalues[folds != i]
    covariates_other_folds <- covariates[folds != i, , drop = F]
    covariates_holdout_fold <- covariates[folds == i, , drop = F]

    groups[folds == i, ] <- JuliaCall::julia_call("learn_discretizations_forest_train", pvalues_other_folds, covariates_other_folds, covariates_holdout_fold, ntrees, nbins, nbasis)
  }

  groups <- as.data.frame(groups)
  groups[] <- lapply(groups, as.factor)
  colnames(groups) <- paste0("group", seq_along(groups))
  return(groups)
}

# given list of adjusted p-values and original p-values
# calculate actual threshold!
padj_to_threshold <- function(padj, pvals, alpha) {
  filtered_pvals <- pvals[padj <= alpha]
  ifelse(length(filtered_pvals) == 0, 0, max(filtered_pvals))
}


mydiv <- function(x, y) {
  ifelse(x == 0, 0,
    ifelse(y == 0, 1, pmin(x / y, 1))
  )
}

# returns a threshold!
filter_pvals_for_optim <- function(pvals, alpha, ntests = length(pvals)) {
  nrjs <- sum(p.adjust(pvals, method = "BH", n = ntests) < alpha)
  min(1, 10 * nrjs / length(pvals))
}

#'  Data-driven threshold of Benjamini Hochberg Procedure
#'
#'       Given pvalues and a nominal significance level alpha, this function returns the
#'   rejection threshold of the Benjamini-Hochberg procedure, i.e. a value t_BH such that p-values with
#'   P_i <= t_BH get rejected by the procedure.
#'
#' @param pvals Numeric, vector of p-values
#' @param alpha Numeric in [0,1], significance level of the multiple testing procedure
#' @param mtests Integer, total number of hypothesis tests; only set this (to non-default) when you know what you are doing!
#'
#' @return A numeric in [0,1], threshold of the BH procedure
#'
#' @examples
#' pvalues <- c(runif(1000), rbeta(1000, 0.5, 7)) # generate some p-values
#' adj_pvalues <- p.adjust(pvalues, method = "BH") # calculate adjusted p-values
#' t_BH <- get_bh_threshold(pvalues, 0.1) # get rejection threshold at alpha=0.1
#' all((pvalues <= t_BH) == (adj_pvalues <= 0.1)) # equivalence of two formulations
#' @export
get_bh_threshold <- function(pvals, alpha, mtests = length(pvals)) {
  m <- length(pvals)
  pvals <- sort(pvals)
  prejected <- which(pvals <= (1:m) / mtests * alpha)
  ifelse(length(prejected) == 0, 0, pvals[prejected[which.max(prejected)]])
}


get_bh_thresholds <- function(unadj_p, filterstat, nbins, alpha) {
  t <- get_bh_threshold(unadj_p, alpha)
  grps <- groups_by_filter(filterstat, nbins)
  pv_list <- split(unadj_p, grps)
  sapply(pv_list, function(ps) max(0, ps[ps <= t]))
}

get_wbh_weights <- function(obj) {
  weighted_pv <- pvalues(obj) / weights(obj, levels_only = FALSE)
  t <- get_bh_threshold(weighted_pv, alpha(obj))
  m_groups <- table(groups_factor(obj))
  grps <- groups_factor(obj)
  pv_list <- split(pvalues(obj), grps)
  ts <- sapply(pv_list, function(ps) max(0, ps[ps <= t]))
  ts * sum(m_groups) / sum(ts * m_groups)
}


lsl_pi0_est <- function(pvalue) {
  n <- length(pvalue)
  ls <- (n:1) / (1 - sort(pvalue))
  ls_diff <- ls[-c(1, 2)] - ls[-c(1, n)]
  index <- min(which(ls_diff > 0)) + 2
  if (index == Inf) {
    pi0 <- 1
  } else {
    pi0 <- min(1, (1 + floor(ls[index])) / n)
  }
  pi0
}

fill_nas_reorder <- function(reduced_vector, nna, order) {
  if (length(nna) == 1 && nna) {
    full_vector <- reduced_vector[order]
  } else {
    full_vector <- rep(NA, length(nna))
    full_vector[nna] <- reduced_vector[order]
  }
  full_vector
}

fill_nas_reorder_dataframe <- function(reduced_df, nna, order) {
  if (length(nna) == 1 && nna) {
    full_df <- reduced_df[order, , drop = F]
  } else {
    full_df <- data.frame(matrix(NA, nrow = length(nna), ncol = ncol(reduced_df)))
    colnames(full_df) <- colnames(reduced_df)
    full_df[nna, ] <- reduced_df[order, ]
  }
  full_df
}

#' regularize weights
#'
#'
#' @param ws weights to be regularized
#' @param penalty Integer, number of groups/strata into which p-values will be split based on covariates.
#' @param gamma regularization parameter, similar function as lambda. gamma = 1 corresponds to no regularization.
#' @return regularized weights
#' @examples
#' TODO
#' covariates <- runif(100)
#' groups <- groups_by_filter(covariates, 10)
#' table(groups)
#' @export
posteori_regularization <- function(ws, m_groups, penalty = "total variation", gamma = 1) {
  if (gamma == 0) {
    uniform_ws <- rep(1, nbins)
    return(uniform_ws)
  } else if (gamma == 1) {
    return(ws)
  }

  if (penalty == "total variation") {
    browser()
    reg_ws <- loess(order(ws) ~ ws, data = ws, span = 0.10) # 10% smoothing span
    reg_ws <- thresholds_to_weights(reg_ws, m_groups)
    return(reg_ws)
  } else if (penalty == "uniform deviation") {
    uniform_ws <- rep(1, length(ws))
    reg_ws <- (1 - gamma) * uniform_ws + gamma * ws
    return(reg_ws)
  }
}
