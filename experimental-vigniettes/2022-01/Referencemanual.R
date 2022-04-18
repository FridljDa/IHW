# Reference manual
# http://bioconductor.org/packages/release/bioc/manuals/IHW/man/IHW.pdf
#library("IHW")
# clear memory
rm(list = ls(all = TRUE))

dirname(rstudioapi::getSourceEditorContext()$path)

script.dir <- "/g/huber/users/fridljand/R/IHW/bioconductor/R"
#script.dir <- "smb://huber/huber/users/fridljand/R/IHW/bioconductor/R"

scripts <- list.files(script.dir)
lapply(scripts, function(script) {source(file.path(script.dir, script))})

## ---TODO delete---
#library(gurobi)
library(rbenchmark)

#Rcpp::sourceCpp("/Users/default/Google Drive/currentDocumants/Studium/Master/3.Semester/Masterarbeit/R/bioconductor/C++/grenador.cpp")
Rcpp::sourceCpp("/g/huber/users/fridljand/R/IHW/bioconductor/C++/grenador.cpp")

set.seed(1)
m <- 20000
nfoldsI <- 2
X <- runif(m, min = 0, max = 20) # covariate
H <- rbinom(m, 1, 0.1) # hypothesis true or false
Z <- rnorm(m, H * X) # Z-score
# .Random.seed <- save.seed
pvalue <- 1 - pnorm(Z) # pvalue

sorted_folds <- sample(1:nfoldsI, m, replace = TRUE)
adjustment_type<- "bonferroni" #"bonferroni" "BH" 

benchmark <- benchmark(
  "algorithm5" = {
    ihw_fdr1 <- ihw(pvalue, X, alpha = .1, lambdas = Inf, adjustment_type = adjustment_type, lp_solver = "algorithm5", folds = sorted_folds) 
  },
  "lpsymphony" = {
    ihw_fdr2 <- ihw(pvalue, X, alpha =.1, lambdas = Inf, adjustment_type = adjustment_type, lp_solver = "lpsymphony", folds = sorted_folds) 
  },
  #"gurobi" = {
  #  ihw_fdr3 <- ihw(pvalue, X, alpha =.1, lambdas = Inf, adjustment_type = adjustment_type, lp_solver = "gurobi", folds = sorted_folds ) 
  # d},
  replications = 1,
  columns = c("test", "replications", "elapsed",
              "relative", "user.self", "sys.self"))
  
  #create plots
  g1 <- IHW::plot(ihw_fdr1) + ggplot2::ggtitle("algorithm5")
  g2 <- IHW::plot(ihw_fdr2)+ ggplot2::ggtitle("lpsymphony")
  #g3 <- IHW::plot(ihw_fdr3)+ ggplot2::ggtitle("gurobi")
  
  #set mutual y range
  y_min <- min(c(weights(ihw_fdr1), weights(ihw_fdr2)))#, weights(ihw_fdr3)))
  y_max <- max(c(weights(ihw_fdr1), weights(ihw_fdr2)))#, weights(ihw_fdr3)))
  
  g1 <- g1 + ggplot2::ylim(y_min, y_max)
  g2 <- g2 + ggplot2::ylim(y_min, y_max)
  #g3 <- g3+ ggplot2::ylim(y_min, y_max)
  
  g_combined <- ggpubr::ggarrange(g1, g2,#g3,
                                  ncol = 1)
  g_combined <- g_combined + ggplot2::ggtitle(adjustment_type)
  g_combined
  
  rejections1 <- IHW::rejections(ihw_fdr1) #algorithm 5
  rejections2 <- IHW::rejections(ihw_fdr2) #lpsymphony
  print(paste("algorithm 5 rejects", rejections1, "hypothesis"))
  print(paste("lpsymphony rejects", rejections2, "hypothesis"))

  if(sum(abs(weights(ihw_fdr1)- weights(ihw_fdr2))) == 0){
    print("lpsymphony and algorithm5 yield the same weights")
  }else{
    print("lpsymphony and algorithm5 do NOT yield the same weights")
  }

  #testthat::test_that("supdif & lu_finder", {
  #  testthat::expect_equal(sum(abs(weights(ihw_fdr1)- weights(ihw_fdr2))),
  #                         0,
  #                         tolerance = 1e-10)
  #})
  

