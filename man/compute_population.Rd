\name{compute_population}
\alias{compute_population}
\title{Population of a market  (total outgoing commuters = total residents)}
\usage{
compute_population(members, row_W)
}
\arguments{
\item{members}{Integer vector.}
\item{row_W}{Numeric vector.}
}
\value{
Numeric scalar.
}
\description{
Pop(Ms) = W(Ms, G) = sum_\{i in Ms\} rowSums(W)[i]
}
