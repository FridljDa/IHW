library(RFCDE)
set.seed(42)
generate_data <- function(n) {
  x_relevant <- matrix(runif(n * 10), n, 10)
  x_irrelevant <- matrix(runif(n * 10), n, 10)
  z <- rnorm(n, rowSums(x_relevant), 1)
  return(list(x = cbind(x_relevant, x_irrelevant), z = z))
}
n_train <- 10000
n_test <- 4
train_data <- generate_data(n_train)
x_train <- train_data$x
z_train <- train_data$z

