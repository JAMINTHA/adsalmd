\name{build_lma_df}
\alias{build_lma_df}
\title{Build a per-market summary data frame}
\usage{
build_lma_df(sol, W, row_W, col_W)
}
\arguments{
\item{sol}{Integer vector.}
\item{W}{Numeric matrix.}
\item{row_W}{Numeric vector.}
\item{col_W}{Numeric vector.}
}
\value{
data.frame with columns: cluster, pop, scss, scds, sc, n_bgus.
}
\description{
Used internally by the objective function and diagnostics.
}
