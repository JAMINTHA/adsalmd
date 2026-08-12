\name{compute_gci}
\alias{compute_gci}
\title{Global Cohesion Index   [Equation 7]}
\usage{
compute_gci(sol, W, row_W, col_W)
}
\arguments{
\item{sol}{Integer vector of assignments (length N).}
\item{W}{Numeric matrix (N x N).}
\item{row_W}{Numeric vector (length N).}
\item{col_W}{Numeric vector (length N).}
}
\value{
Numeric scalar GCI.
}
\description{
GCI = sum_\{g in G\} CI(g) Automatically uses the fast sub-matrix version for N >= 50 (compute_gci_fast, defined later in this file) and the element-wise version for small N where clarity matters more.
}
