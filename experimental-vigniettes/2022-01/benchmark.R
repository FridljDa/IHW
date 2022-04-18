#benchmark
# Reference manual
# http://bioconductor.org/packages/release/bioc/manuals/IHW/man/IHW.pdf
#library("IHW")
# clear memory
rm(list = ls(all = TRUE))

dirname(rstudioapi::getSourceEditorContext()$path)

#script.dir <- "/Users/default/Google Drive/currentDocumants/Studium/Master/3.Semester/Masterarbeit/R/bioconductor/R"
script.dir <- "/g/huber/users/fridljand/R/IHW/bioconductor/R"

scripts <- list.files(script.dir)
lapply(scripts, function(script) {source(file.path(script.dir, script))})
## ---TODO delete---
#library(gurobi)
library(rbenchmark)
library(magrittr)
library(ggplot2)

#Rcpp::sourceCpp("/Users/default/Google Drive/currentDocumants/Studium/Master/3.Semester/Masterarbeit/R/bioconductor/C++/grenador.cpp")
Rcpp::sourceCpp("/g/huber/users/fridljand/R/IHW/bioconductor/C++/grenador.cpp")

set.seed(1)

ms <- 200000
nfoldsI <- 5

benchmarks <- lapply(ms, function(m){
  
  X <- runif(m, min = 0, max = 2.5) # covariate
  H <- rbinom(m, 1, 0.1) # hypothesis true or false
  Z <- rnorm(m, H * X) # Z-score
  # .Random.seed <- save.seed
  pvalue <- 1 - pnorm(Z) # pvalue
  
  sorted_folds <- sample(1:nfoldsI, m, replace = TRUE) 
  sorted_folds <- as.factor(sorted_folds)
  adjustment_type<- "BH" #"bonferroni" "BH" 
  
  benchmark <- benchmark(
    "lpsymphony" = {
      ihw_fdr2 <- ihw(pvalue, X, alpha =.1, lambdas = Inf, covariate_type = "ordinal", adjustment_type = adjustment_type, lp_solver = "lpsymphony", folds = sorted_folds) 
    },
    "algorithm5" = {
      ihw_fdr1 <- ihw(pvalue, X, alpha = .1, lambdas = Inf, covariate_type = "ordinal", adjustment_type = adjustment_type, lp_solver = "algorithm5", folds = sorted_folds) 
    },
    replications = 1,
    columns = c("test",  "elapsed"
                #,"replications","relative", "user.self", "sys.self"
                ))

  benchmark$m <- m
  return(benchmark)
})
benchmarks <- data.table::rbindlist(benchmarks)

g1 <- ggplot(benchmarks, aes(x = m, y = elapsed, color = test)) +
  geom_line(size = 1.5) +
  scale_x_continuous(trans='log10')
ggsave("/g/huber/users/fridljand/R/IHW/benchmark.png", g1#, height = 4, width = 8
       )
