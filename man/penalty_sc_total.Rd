\name{penalty_sc_total}
\alias{penalty_sc_total}
\title{Total self-containment penalty   [Equation 13 summed]}
\usage{
penalty_sc_total(sc_vec, SC_min, r)
}
\arguments{
\item{sc_vec}{Numeric vector, SC values for all markets.}
\item{SC_min}{Numeric.}
\item{r}{Integer.}
}
\value{
Numeric scalar.
}
\description{
PSC(x) = sum_\{Ms in M\} PSC(Ms)
}
