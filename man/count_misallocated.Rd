\name{count_misallocated}
\alias{count_misallocated}
\title{Count cohesion-misallocated BGUs in a solution}
\usage{
count_misallocated(sol, W, adj, row_W, col_W, sc_tol = 0.05)
}
\arguments{
\item{sol}{Integer vector.}
\item{W}{Numeric matrix.}
\item{adj}{Integer/logical matrix.}
\item{row_W}{Numeric vector.}
\item{col_W}{Numeric vector.}
\item{sc_tol}{Numeric, tolerated SC reduction.}
}
\value{
Integer count.
}
\description{
See R/08_diagnostics.R for details.
}
