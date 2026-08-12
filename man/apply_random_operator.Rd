\name{apply_random_operator}
\alias{apply_random_operator}
\title{Apply one randomly chosen operator (uniform over 1–10)}
\usage{
apply_random_operator(sol, W, adj, row_W, col_W, best_sol = NULL, ...)
}
\arguments{
\item{sol}{Integer vector.}
\item{W}{Numeric matrix.}
\item{adj}{Integer/logical matrix.}
\item{row_W}{Numeric vector.}
\item{col_W}{Numeric vector.}
\item{best_sol}{See Description.}
\item{...}{Additional arguments (see Description / source).}
}
\value{
Perturbed integer vector (may equal sol if move was infeasible).
}
\description{
See R/06_operators.R for details.
}
