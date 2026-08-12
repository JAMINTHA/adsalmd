\name{compute_scds}
\alias{compute_scds}
\title{Demand-side self-containment for a market   [Equation 11]}
\usage{
compute_scds(members, W, col_W)
}
\arguments{
\item{members}{Integer vector.}
\item{W}{Numeric matrix.}
\item{col_W}{Numeric vector.}
}
\value{
Numeric in [0, 1].
}
\description{
SCDS(Ms) = W(Ms, Ms) / W(G, Ms)
}
