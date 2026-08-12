\name{compute_ci}
\alias{compute_ci}
\title{Cohesion Index for a single BGU   [Equation 6]}
\usage{
compute_ci(g, M_g, W, row_W, col_W)
}
\arguments{
\item{g}{Integer, BGU index.}
\item{M_g}{Integer vector, other BGUs in g's market (M \\ \{g\}).}
\item{W}{Numeric matrix, OD matrix.}
\item{row_W}{Numeric vector, rowSums(W).}
\item{col_W}{Numeric vector, colSums(W).}
}
\value{
Numeric scalar >= 0.
}
\description{
CI(g) = W(g, M\\g)^2 / (W(g,G) * W(G, M\\g))       + W(M\\g, g)^2 / (W(M\\g,G) * W(G, g)) Returns 0 for singleton clusters.
}
