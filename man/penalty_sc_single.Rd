\name{penalty_sc_single}
\alias{penalty_sc_single}
\title{Self-containment penalty for one market   [Equation 13]}
\usage{
penalty_sc_single(sc, SC_min, r)
}
\arguments{
\item{sc}{Numeric, self-containment of market Ms.}
\item{SC_min}{Numeric, minimum SC threshold.}
\item{r}{Integer, steepness exponent (typically 3 or 4).}
}
\value{
Numeric scalar >= 0.
}
\description{
PSC(Ms) = 1 - (SC(Ms) / SC_min)^r   if SC(Ms) < SC_min         = 0                           otherwise
}
