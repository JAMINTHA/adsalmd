\name{is_bgu_misallocated}
\alias{is_bgu_misallocated}
\title{Check whether a single BGU is cohesion-misallocated}
\usage{
is_bgu_misallocated(g, sol, W, adj, row_W, col_W, sc_tol = 0.05)
}
\arguments{
\item{g}{Integer, BGU index.}
\item{sol}{Integer vector.}
\item{W}{Numeric matrix.}
\item{adj}{Integer/logical matrix.}
\item{row_W}{Numeric vector.}
\item{col_W}{Numeric vector.}
\item{sc_tol}{Numeric, maximum tolerated SC reduction (default 0.05).}
}
\value{
Logical.
}
\description{
A BGU g is misallocated if there exists a neighbouring LMA such that reassigning g to it would increase CI(g) without causing a substantial loss in self-containment of either the origin or destination market (paper: "without substantial loss in self- containment").
}
