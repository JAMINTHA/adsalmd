\name{compute_base_fitness}
\alias{compute_base_fitness}
\title{Base fitness function used for cross-method comparison   [Equation 8]}
\usage{
compute_base_fitness(sol, W, row_W, col_W)
}
\arguments{
\item{sol}{Integer vector.}
\item{W}{Numeric matrix.}
\item{row_W}{Numeric vector.}
\item{col_W}{Numeric vector.}
}
\value{
Numeric scalar.
}
\description{
fbase(x) = m * GCI(x) This is the common performance measure used by GEA and MSA.
}
