# clear memory
rm(list = ls(all = TRUE))
library(ggplot2)
library(cowplot)
library(magrittr)
library(dplyr)

dirname(rstudioapi::getSourceEditorContext()$path)

#script.dir <- "/Users/default/Google Drive/currentDocumants/Studium/Master/3.Semester/Masterarbeit/R/bioconductor/R"
script.dir <- "/g/huber/users/fridljand/R/IHW/bioconductor/R"

scripts <- list.files(script.dir)
lapply(scripts, function(script) {source(file.path(script.dir, script))})



#' Evaluate multiple testing procedure
#
#' @param Hs Vector with indicators of alternatives (1) and true nulls (0)
#' @param rjs Vector with indicator of rejected hypotheses
#' @return Data frame with columns `rjs` (total rejections), `pow` (Power), `FDP` (False discovery proportion)
#' @export
fdp_eval <- function(Hs, rjs){
  rjs_total <- sum(rjs)
  pow <- sum(rjs*Hs)/max(1,sum(Hs))
  FDP <- sum(rjs*(1-Hs))/max(1,rjs_total)
  FWER <- sum( (1-Hs)*rjs) > 0
  data.frame(rjs=rjs_total, pow=pow, FDP=FDP, FWER=FWER)
}

#' Simulation: Misspecified conditional Beta-uniform mixture model
#
#' @param m Number of hypotheses (default: m=10000)
#' @param mus_slope Numeric (default:1.5) parameter bar(beta) in equation (12)
#' @param one_sided_tests Bool (default:FALSE), if true adds some nulls that are strictly superuniform
#' @param prob_one_sided Numeric (default:0.25) proportion of nulls that are strictly superuniform
#'
#' @return Data frame with columns `Hs` (null or alternative), `Ps` (p-value), `Xs` (side-information),
#'         `alphas` (parameters of alternative distribution), `pi1s` (probabilities of being from the alternative distribution),
#'         `oracle_lfdr` (oracle local fdr)
#' @export
beta_unif_sim <- function(m=10000, mus_slope=1.5, one_sided_tests=FALSE, prob_one_sided=0.25){
  Xs <- matrix(runif(m*2, 0,1), ncol=2)
  colnames(Xs) <- c("X1", "X2")
  
  pi1s <- ifelse( Xs[,1]^2 + Xs[,2]^2 <= 1, 0.02, 0.4)
  mus <- pmax(1.3, sqrt(Xs) %*% c(1,1)*mus_slope)
  
  mu_alphas <- 1/mus
  
  Hs <- stats::rbinom(m, size=1, prob=pi1s)
  Ps <- stats::runif(m)*(1-Hs) + stats::rbeta(m, mu_alphas, 1)*Hs
  Xs <- data.frame(Xs)
  if (one_sided_tests){
    Hs_alt <-  1- (1-Hs)*stats::rbinom(m, size=1, prob=prob_one_sided)
    Ps[Hs_alt == 0] <- stats::rbeta(sum(Hs_alt == 0), 1, 0.5)
    oracle_lfdr_null <- (1-pi1s)*( 1-prob_one_sided + prob_one_sided*stats::dbeta(Ps, 1, 0.5) )
  } else{
    oracle_lfdr_null <- 1-pi1s
  }
  oracle_lfdr_alternative <- pi1s*dbeta(Ps, mu_alphas, 1)
  oracle_lfdr <- oracle_lfdr_null/(oracle_lfdr_null+oracle_lfdr_alternative)
  list(Xs=Xs, Ps=Ps, Hs=Hs, alphas=mu_alphas, pi1s=pi1s, oracle_lfdrs=oracle_lfdr)
}

#Rcpp::sourceCpp("/Users/default/Google Drive/currentDocumants/Studium/Master/3.Semester/Masterarbeit/R/bioconductor/C++/grenador.cpp")
Rcpp::sourceCpp("/g/huber/users/fridljand/R/IHW/bioconductor/C++/grenador.cpp")

#generate data
m <- 20000
nfoldsI <- 5
sorted_folds <- sample(1:nfoldsI, m, replace = TRUE)
adjustment_type<- "bonferroni" #"bonferroni" "BH" 

betamix_sim_combs <- expand.grid(mu_slope = seq(1, 3, length=4), seed = 1:10)

fdp_eval <- apply(betamix_sim_combs,1, function(experiment){
  mu_slope <- experiment[[1]]
  seed <- experiment[[2]]
  sim <- beta_unif_sim(m=m, mus_slope=mu_slope,
                       one_sided_tests=FALSE, prob_one_sided = 0.25)
  Xs <- sim$Xs
  Xs <- unlist(Xs[,1])
  Ps <- sim$Ps
  Hs <- sim$Hs
  mu_alphas <- sim$alphas
  pi1s <- sim$pi1s
  oracle_lfdrs <- sim$oracle_lfdrs
  
  ihw_fdr1 <- ihw(Ps, Xs, alpha = .1, lambdas = Inf, adjustment_type = adjustment_type, lp_solver = "algorithm5", folds = sorted_folds) 
  ihw_fdr2 <- ihw(Ps, Xs, alpha =.1, lambdas = Inf, adjustment_type = adjustment_type, lp_solver = "lpsymphony", folds = sorted_folds) 
 
  fdp_eval_alg <- fdp_eval(Hs,  IHW::rejected_hypotheses(ihw_fdr1))
  fdp_eval_alg$method <- "algorithm 5"
  fdp_eval_lp <- fdp_eval(Hs,  IHW::rejected_hypotheses(ihw_fdr2))
  fdp_eval_lp$method <- "lpsymphony"
  
  fdp_eval <- rbind(fdp_eval_alg, fdp_eval_lp)
  fdp_eval$mu_slope <- mu_slope

  return(fdp_eval)
}) 

fdp_eval <- data.table::rbindlist(fdp_eval)
#TODO mean over seed
fdp_eval <-  fdp_eval %>% 
  group_by(method, mu_slope) %>%
  summarize(FDR = mean(FDP, na.rm=TRUE),
            Power=mean(pow, na.rm=TRUE),
            FWER = mean(FWER, na.rm = TRUE),
            n=sum(!is.na(FDP)),
            se_power = sd(pow, na.rm=TRUE)/sqrt(n)) %>%
  arrange(mu_slope, desc(Power)) %>% 
  ungroup()

fdr_plot <- ggplot(fdp_eval, aes(x=mu_slope, y=FWER,shape=method, col=method)) + 
  geom_line() + geom_point() + 
  ylim(0, 0.15) + 
  geom_hline(yintercept=0.1, linetype=2)+
  xlab("mu_slope")+
  theme_cowplot() +
  theme(legend.position = "none")

fdr_plot

power_plot <- ggplot(fdp_eval, aes(x=mu_slope, y=Power, shape=method, col=method)) + 
  geom_line() + geom_point() + 
  xlab("mu_slope")+
  scale_y_log10() + 
  ylab("Power (log scale)") + 
  theme_cowplot()  +
  theme(legend.title=element_blank())  #+  theme(legend.position = "none")

power_plot

ggsave("/g/huber/users/fridljand/R/IHW/fdr_plot.png", fdr_plot#, height = 4, width = 8
)
