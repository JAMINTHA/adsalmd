\name{compute_sc}
\alias{compute_sc}
\title{Self-containment for a market   [Equation 12]}
\usage{
compute_sc(members, W, row_W, col_W)
}
\arguments{
\item{members}{Integer vector.}
\item{W}{Numeric matrix.}
\item{row_W}{Numeric vector.}
\item{col_W}{Numeric vector.}
}
\value{
Named list: scss, scds, sc.
}
\description{
SC(Ms) = min(SCSS(Ms), SCDS(Ms))
}
