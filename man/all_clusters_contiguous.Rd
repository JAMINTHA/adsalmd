\name{all_clusters_contiguous}
\alias{all_clusters_contiguous}
\title{Check that EVERY cluster in a solution is contiguous}
\usage{
all_clusters_contiguous(sol, adj)
}
\arguments{
\item{sol}{Integer vector of assignments.}
\item{adj}{Integer/logical matrix.}
}
\value{
Logical.
}
\description{
Returns FALSE as soon as any cluster fails the BFS test.
}
