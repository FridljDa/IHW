library(tictoc)
library(stringr)
devtools::load_all()

m <- 1e5
r <- 1
ntrees <- 100

# regular simulation
wasserman_normal_sim <- function(m, pi0, xi_min, xi_max, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  X <- runif(m, min = xi_min, max = xi_max)
  H <- rbinom(m, 1, 1 - pi0)
  Z <- rnorm(m, H * X)
  pvalue <- 1 - pnorm(Z)
  simDf <- data.frame(pvalue = pvalue, filterstat = X, H = H, Z = Z, noise = 0)
}

sim <- wasserman_normal_sim(m, 0.85, 0, 3, seed = 1)

cov_num_vec <- sim$filterstat



#set forest parameter
set.seed(1)

# Capture the console output in a string variable
output <- capture.output({
  ihw_forest_num_vec <-  replicate(r, ihw(sim$pvalue, cov_num_vec, .1, stratification_method = "forest", ntrees = ntrees))
})



output <- paste0(output, collapse = "")
# Vector of patterns
beginning <- c("group_by_forest","lpsymphony", "constr_matrix","fdrtool::gcmlcm","filtered_sorted_pvalues","sorted_weights","sorted_weighted_pvalues","sorted_adj_p")

# Initialize an empty list to store the matching times
matching_times <- list()

# Loop through each pattern
for (pattern in beginning) {
  # Create the pattern string dynamically
  pattern_string <- paste0(pattern, ":\\s+(\\d+\\.\\d+) sec")
  
  # Find matches in the string using the pattern
  matches <- str_match_all(output, pattern_string)
  
  # Extract the matching times and convert them to numeric
  times <- sapply(matches, function(x) as.numeric(x[, 2]))
  
  # Sum the matching times and store the result in the list
  matching_times[[pattern]] <- sum(times)
}




cat("Wassermann simulation set-up:\n")
cat("Number of hypothesis:", m, "\n")
cat("Number of replicates:", r, "\n")
cat("Number of trees:", ntrees, "\n")
cat("\n")
cat("\n")



cat("Measured times: \n")
# Print the matching times
for (pattern in beginning) {
  cat(pattern, matching_times[[pattern]], "\n")
}

row <- data.frame(matching_times) 
row$m <- m
row$r <- r
row$ntress <- ntrees

csv_fname <- "measure_time_scripts/data/measured_time.csv"
write.table(row, 
            csv_fname, 
            sep = ",", 
            #col.names = FALSE,
            row.names = FALSE,
            col.names = !file.exists(csv_fname), 
            append = T)
