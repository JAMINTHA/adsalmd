\name{estimate_t0}
\alias{estimate_t0}
\title{Estimate initial temperature T0   [Equation 3]}
\usage{
estimate_t0(init_sol, W, adj, row_W, col_W, SC_min, Pop_min, Pop_max, Pop_tar, r, P0 = 0.8, n_samples = 100)
}
\arguments{
\item{init_sol}{Integer vector, starting solution.}
\item{W}{Numeric matrix.}
\item{adj}{Integer/logical matrix.}
\item{row_W}{Numeric vector.}
\item{col_W}{Numeric vector.}
\item{SC_min}{Numeric.}
\item{Pop_min}{Numeric.}
\item{Pop_max}{See Description.}
\item{Pop_tar}{See Description.}
\item{r}{Integer, penalty exponent.}
\item{P0}{Numeric, target acceptance probability (default 0.8).}
\item{n_samples}{Integer, perturbations used (default 100).}
}
\value{
Numeric scalar T0 > 0.
}
\description{
T0 = -E[|Delta F|] / ln(P0) E[|Delta F|] is the mean absolute fitness change across n_samples random perturbations of the initial solution.
}
