# functions to assist in working with weights


total_variation <- function(ws) {
  sum(abs(diff(ws)))
}

uniform_deviation <- function(ws) {
  sum(abs(ws - 1))
}

thresholds_to_weights <- function(ts, m_groups) {
  if (length(ts) != length(m_groups)) stop("incosistent number of elements for ts and m_groups")

  nbins <- length(ts)
  m <- sum(m_groups)

  if (all(ts == .0)) {
    rep(1, nbins)
  } else {
    ts * m / sum(m_groups * ts)
  }
}

thresholds_to_weights_full <- function(ts) {
  m <- length(ts)

  if (all(ts == .0)) {
    rep(1, m)
  } else {
    ts * m / sum(ts)
  }
}


# normalize weights so that their sum will be equal to length(ws)
# also make sure no negative weights appear
normalize_weights <- function(ws) {
  abs(ws) * length(ws) / sum(abs(ws))
}

#' regularize weights
#'
#'
#' @param ws weights to be regularized
#' @param penalty Integer, number of groups/strata into which p-values will be split based on covariate.
#' @param gamma regularization parameter, similar function as lambda. gamma = 1 corresponds to no regularization.
#' @return regularized weights
#' @examples
#' TODO
#' covariates <- runif(100)
#' groups <- groups_by_filter(covariates, 10)
#' table(groups)
#' @export
regularize_weights <- function(ws, m_groups, penalty = "total variation", gamma = 1) {
  #if (gamma == 0) {
  #  uniform_ws <- rep(1, length(ws))
  #  return(uniform_ws)
  #} else if (gamma == 1) {
  #  return(ws)
  #}

  if (penalty == "total variation") {
    # try more advanced later
    #browser()  
    #ws_df <- data.frame(ws = ws, order_ws = order(ws))
    #reg_ws_lo <- loess(order_ws ~ ws, data=ws_df, degree =1, span = 1) # 10% smoothing span
    #reg_ws <- predict(reg_ws_lo, ws_df)

    G <- length(ws)
    identity_matrix <- diag(G)
    moving_window_matrix <- diag(G)
    for (i in 1:G) {
      for (j in 1:G) {
        if ((j == 1 & i %in% 1:2) |
          (j == G & i %in% (G - 1):G)) {
          moving_window_matrix[i, j] <- 0.5
        } else if (abs(i - j) <= 1) {
          moving_window_matrix[i, j] <- 1 / 3
        }
      }
    }

    mixed_matrix<- (1 - gamma) * moving_window_matrix + gamma * identity_matrix
    reg_ws <-   mixed_matrix %*% ws
    
    reg_ws <- thresholds_to_weights(reg_ws, m_groups)
    return(reg_ws)
  } else if (penalty == "uniform deviation") {
    uniform_ws <- rep(1, length(ws))
    reg_ws <- (1 - gamma) * uniform_ws + gamma * ws
    return(reg_ws)
  }
}

#toy_ws <- c(1, 1.5, 0.5,1,1,1,1,1,1,1,1)
#m_groups <- rep(200, length(toy_ws))
#toy_ws <- thresholds_to_weights(toy_ws, m_groups)
#regularize_weights(toy_ws, m_groups, gamma = 0.5)
