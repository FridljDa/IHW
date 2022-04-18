library(ggplot2)
library(dplyr)       # data wrangling
library(rpart)       # performing regression trees
library(rpart.plot)  # plotting regression trees
library(ipred)       # bagging
library(caret)       # bagging

#Data generation as in 5.3
sim <- IHWStatsPaper::beta_unif_sim(m=10000, mus_slope=1,one_sided_tests=FALSE, prob_one_sided = 0.25)
Xs <- sim$Xs
Ps <- sim$Ps
Hs <- sim$Hs

#Data generation as in Reference manual
#m <- 20000
#X1 <- runif(m, min = 0, max = 2.5) # covariate
#X2 <- runif(m, min = 0, max = 2.5) # covariate
#Xs <- as.data.frame(cbind(X1, X2) )
#Hs <- rbinom(m, 1, 0.1) # hypothesis true or false
#Z <- rnorm(m, Hs* (X1+X2)) # Z-score
#Ps <- 1 - pnorm(Z) # pvalue

#extremely informative, seperated data
#m <- 20000
#nbinsX <- 5
#X1 <- rep(seq_len(nbinsX), each = m/nbinsX) + rnorm(m)# covariate 1
#X2 <- rep(seq_len(nbinsX), each = m/nbinsX) + rnorm(m)# covariate 2
#Xs <- as.data.frame(cbind(X1, X2) )
#Hs <- rbinom(m, 1, 0.1) # hypothesis true or false
#Z <- rnorm(m, Hs* X1*X2) # Z-score
#Ps <- 1 - pnorm(Z) # pvalue

#plot input
ggplot(data = cbind(Xs, Ps),
       aes(X1, X2, color = Ps)) + geom_point()

#blocks
#Xs_disc <- interaction(cut(Xs[,1],5), cut(Xs[,2],5))
#ihw_nmeth_fit <- IHW::ihw(Ps, Xs_disc, alpha=0.1, lambdas=Inf)
#ihw_nmeth_res <- IHW::rejected_hypotheses(ihw_nmeth_fit)
#ihw_nmeth_fdp <- IHWStatsPaper::fdp_eval(Hs, ihw_nmeth_res)

#ws <- ihw_nmeth_fit@weights
#groups <- levels(ihw_nmeth_fit@df$group)
#ws_fold1 <- ws[,1]

#data1 = cbind(Xs, weights = as.numeric(IHW::weights(ihw_nmeth_fit)))

#g1 <- ggplot(data = data1,
#       aes(X1, X2, color = weights)) + geom_point()

#regression tree
tau <-  0.9 #https://www.biorxiv.org/content/10.1101/035675v4.full.pdf fixed thresholds between 0.8 and 0.95 are often used. 
Ys <- (Ps>= tau)/(1-tau)
 
data2 <- cbind(Ys, Xs, Ps)

ggplot(data = data2, aes(X1, X2, color = Ps)) + geom_point()

m1 <- rpart(
  formula = Ps ~ X1 + X2,
  data    = data2,
  method  = "anova",
  control = rpart.control(minbucket = 20, cp = 0.02, xval = 10)
)#control = list(minsplit = 10, maxdepth = 12, xval = 10)
#rpart.control(minbucket = 20, cp = 0.0002, maxsurrogate = 0, usesurrogate = 0, xval = 10)
rpart.plot(m1)
