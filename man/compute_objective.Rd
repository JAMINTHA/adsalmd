\name{compute_objective}
\alias{compute_objective}
\title{AdSA-ALMD objective function F(x)   [Equation 5]}
\usage{
compute_objective(sol, W, row_W, col_W, SC_min, Pop_min, Pop_max = Inf, Pop_tar = NULL, r, adj, check_contiguity = TRUE)
}
\arguments{
\item{sol}{Integer vector of assignments (length N).}
\item{W}{Numeric matrix (N x N).}
\item{row_W}{Numeric vector (length N).}
\item{col_W}{Numeric vector (length N).}
\item{SC_min}{Numeric.}
\item{Pop_min}{Numeric.}
\item{Pop_max}{See Description.}
\item{Pop_tar}{See Description.}
\item{r}{Integer, penalty exponent.}
\item{adj}{Integer/logical matrix (N x N).}
\item{check_contiguity}{Logical; set FALSE only for T0 estimation.}
}
\value{
Numeric scalar (0 if any cluster is non-contiguous).
}
\description{
F(x) = (GCI(x) - P(x)) * prod_\{i=1\}^\{m\} C(Mi) The contiguity product makes F = 0 whenever any market is non-contiguous, acting as a hard feasibility constraint. NOTE: all_clusters_contiguous() is defined in 04_contiguity.R.
}
