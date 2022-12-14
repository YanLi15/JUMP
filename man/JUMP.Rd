% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/JUMP.R
\name{JUMP}
\alias{JUMP}
\title{False discovery rate control for replicability analysis}
\usage{
JUMP(pvals1, pvals2, alpha = 0.05, lambda = seq(0.01, 0.8, 0.01))
}
\arguments{
\item{pvals1}{a numeric vector of the p-values from study 1}

\item{pvals2}{a numeric vector of the p-values from study 2}

\item{alpha}{the significance level}

\item{lambda}{the value of the the tuning parameter to estimate}
}
\value{
a list with the elements
\item{p.max} the maximum p-value statistics for replicability analysis
\item{jump.thr} the threshold to control FDR of replicability null hypotheses estimated for the maximum p-value
}
\description{
The function implements an FDR control procedure for replicability analysis of spatially variable genes in two spatially resolved transcriptomic studies. The maximum p-value statistic is calculated from the paired p-values, and the threshold for rejection of replicability null hypotheses is estimated using a step-up procedure.
}
\author{
Yan Li
}
