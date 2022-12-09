# source("./R/JUMP.R")
library(JUMP)
library(splines)

## Data simulation
n = 100 # number of replications
m = 10000
xi00 = 0.6
xi01 = 0.175
xi10 = 0.175
xi11 = 0.05
mu1 = 2
mu2 = 2
sigma1 = 1
sigma2 = 1
alphas <- 0.05

fdp_JUMP  <- matrix(NA, n, length(alphas))
fdp_BH    <- matrix(NA, n, length(alphas))
fdp_MaxP  <- matrix(NA, n, length(alphas))

pd_JUMP  <- matrix(NA, n, length(alphas))
pd_BH    <- matrix(NA, n, length(alphas))
pd_MaxP  <- matrix(NA, n, length(alphas))

for(i in 1:n){
  h = sample(0:3, m, replace = TRUE, prob = c(xi00, xi01, xi10, xi11))
  states1 = rep(0, m)
  states1[which(h==2|h==3)] = 1
  states2 = rep(0, m)
  states2[which(h==1|h==3)] = 1

  stat1 = rnorm(m, states1*mu1, sigma1)
  stat2 = rnorm(m, states2*mu2, sigma2)

  p1 = 1 - pnorm(stat1, mean = 0, sd = sigma1)
  p2 = 1 - pnorm(stat2, mean = 0, sd = sigma2)

  truth <- states1 * states2

  # BH
  pvals1.bh <- p.adjust(p1, method = "BH")
  pvals2.bh <- p.adjust(p2, method = "BH")
  
  # MaxP
  maxp <- pmax(p1, p2)
  pvals.maxp <- p.adjust(maxp, method = "BH")
  
  for(j in 1:length(alphas)){
    alpha = alphas[j]

    # JUMP
    jump.obj <- JUMP(p1, p2, alpha)
    jump.thr <- jump.obj$jump.thr
    p.max <- jump.obj$p.max
    fdp_JUMP[i,j] <- sum(p.max <= jump.thr & !truth)/max(sum(p.max <= jump.thr), 1)
    pd_JUMP[i,j]  <- sum(p.max <= jump.thr & truth) / sum(truth)

    # MaxP
    fdp_MaxP[i,j] <- sum(pvals.maxp <= alpha & !truth)/ max(sum(pvals.maxp <= alpha), 1)
    pd_MaxP[i,j]  <- sum(pvals.maxp <= alpha & truth) / sum(truth)

    # BH
    fdp_BH[i,j] <- sum(pvals1.bh <= alpha & pvals2.bh <= alpha & !truth)/max(sum(pvals1.bh <= alpha & pvals2.bh <= alpha), 1)
    pd_BH[i,j] <- sum(pvals1.bh <= alpha & pvals2.bh <= alpha & truth) / sum(truth)
  }
  print(i)
}

fdr_MaxP = colMeans(fdp_MaxP, na.rm = TRUE)
fdr_BH = colMeans(fdp_BH, na.rm = TRUE)
fdr_JUMP = colMeans(fdp_JUMP, na.rm = TRUE)

power_ MaxP = colMeans(pd_MaxP, na.rm = TRUE)
power_BH = colMeans(pd_BH, na.rm = TRUE)
power_JUMP = colMeans(pd_JUMP, na.rm = TRUE)
