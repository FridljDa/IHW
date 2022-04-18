rm(list = ls(all = TRUE))

#install profvis
if (!require("profvis")) install.packages("profvis")

#load code locally, so that profvis can step into
script.dir <-   "/g/huber/users/fridljand/R/IHW/bioconductor/R"
scripts <- list.files(script.dir)
lapply(scripts, function(script) {source(file.path(script.dir, script))})

#generate simple data with 2 Million p values
set.seed(1)
m <- 2000000
nfoldsI <- 2
X <- runif(m, min = 0, max = 2.5) # covariate
H <- rbinom(m, 1, 0.1) # hypothesis true or false
Z <- rnorm(m, H * X)
pvalue <- 1 - pnorm(Z) # pvalue

sorted_folds <- sample(1:nfoldsI, m, replace = TRUE)
adjustment_type<- "BH" 

#run profvis
profvis::profvis({
   ihw_fdr1 <- ihw(pvalue, X, alpha = .1, lambdas = Inf, adjustment_type = adjustment_type, lp_solver = "lpsymphony", folds = sorted_folds) #lpsymphony implementation
})

#7100ms overall runtime using "as.factor"
#4700ms overall runtime using cut
