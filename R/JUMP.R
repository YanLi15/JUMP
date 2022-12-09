#' Estimate the proportions of the composite null hypotheses from the data
#'
#' The function provides conservative estimates of the proportions defined in the four-group model
#'
#' @param pvals1 a numeric vector of the p-values from study 1
#' @param pvals2 a numeric vector of the p-values from study 2
#' @param lambda the value of the the tuning parameter to estimate
#'
#' @return a list of estimates of the proportions of four groups corresponding to the joint hidden states
#'
#' @importFrom splines bs
#' @export

xi.est <- function(pvals1, pvals2, lambda = seq(0.01, 0.8, 0.01)){
  require(splines)
  m = length(pvals1)

  xi00.hat = c(); pi0.hat1 = c(); pi0.hat2 = c()
  for (i in 1:length(lambda)) {
    xi00.hat[i] <- sum(pvals1[1:m]>=lambda[i] & pvals2>=lambda[i]) / (m*(1-lambda[i])^2)
    pi0.hat1[i] <- sum(pvals1>=lambda[i]) / (m * (1-lambda[i]))
    pi0.hat2[i] <- sum(pvals2>=lambda[i]) / (m * (1-lambda[i]))
  }
  # fitting a cubic spline following Storey and Tibshirani (2003)
  fit1 <- lm(xi00.hat~., data.frame(cbind(xi00.hat, bs(lambda))))
  fit2 <- lm(pi0.hat1~., data.frame(cbind(pi0.hat1, bs(lambda))))
  fit3 <- lm(pi0.hat2~., data.frame(cbind(pi0.hat2, bs(lambda))))

  pred1 = predict(fit1)
  diff1 = abs(diff(pred1))

  pred2 = predict(fit2)
  diff2 = abs(diff(pred2))

  pred3 = predict(fit3)
  diff3 = abs(diff(pred3))

  xi00.hat = as.numeric(pred1[which.min(diff1)])
  xi01.hat = max(0, as.numeric(pred2[which.min(diff2)] - xi00.hat))
  xi10.hat = max(0, as.numeric(pred3[which.min(diff3)] - xi00.hat))
  xi11.hat = as.numeric(1 - xi00.hat -xi01.hat - xi10.hat)

  return(list(xi00.hat = xi00.hat,
              xi01.hat = xi01.hat,
              xi10.hat = xi10.hat,
              xi11.hat = xi11.hat))
}


#' False discovery rate control for replicability analysis
#'
#' The function implements an FDR control procedure for replicability analysis of spatially
#' variable genes in two spatially resolved transcriptomic studies. The maximum p-value
#' statistic is calculated from the paired p-values, and the threshold for rejection of
#' replicability null hypotheses is estimated using a step-up procedure.
#'
#' @param pvals1 a numeric vector of the p-values from study 1
#' @param pvals2 a numeric vector of the p-values from study 2
#' @param alpha the significance level
#' @param lambda the value of the the tuning parameter to estimate
#'
#' @return a list with the elements
#' \item{p.max} the maximum p-value statistics for replicability analysis
#' \item{jump.thr} the threshold to control FDR of replicability null hypotheses estimated for the maximum p-value
#'
#' @export

JUMP <- function(pvals1, pvals2, alpha = 0.05, lambda = seq(0.01, 0.8, 0.01)){
  xi.hat <- xi.est(pvals1, pvals2, lambda)
  xi00.hat <- xi.hat$xi00.hat
  xi01.hat <- xi.hat$xi01.hat
  xi10.hat <- xi.hat$xi10.hat

  m = length(pvals1)
  p.max <- pmax(pvals1, pvals2)

  thresholds <- sort(p.max, decreasing = TRUE)
  fdp.hat <- c()
  for (i in 1:length(thresholds)){
    thr = thresholds[i]
    fdp.hat[i] <- m*(xi00.hat*thr^2 + (xi01.hat+xi10.hat)*thr)/max(sum(p.max <= thr), 1)
    if(fdp.hat[i] <= alpha) break
  }

  return(list(p.max = p.max, jump.thr = thr))
}
