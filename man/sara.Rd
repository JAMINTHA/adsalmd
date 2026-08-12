\name{sara}
\alias{sara}
\title{SARA: Spatially Adaptive Random Aggregation  [Algorithm 2]}
\usage{
sara(init_sol, W, adj, row_W, col_W, SC_min, Pop_min, Pop_max, n_low = 3, n_high = 3, max_iter = 10000, active_bgus = NULL)
}
\arguments{
\item{init_sol}{Integer vector, starting partition (e.g. 1:N).}
\item{W}{Numeric matrix (N x N).}
\item{adj}{Integer/logical matrix (N x N).}
\item{row_W}{Numeric vector, rowSums(W).}
\item{col_W}{Numeric vector, colSums(W).}
\item{SC_min}{Numeric, minimum self-containment threshold.}
\item{Pop_min}{Numeric, minimum population size threshold.}
\item{Pop_max}{See Description.}
\item{n_low}{Integer, candidates to sample focal from (default 3).}
\item{n_high}{Integer, candidates to sample target from (default 3).}
\item{max_iter}{Integer, safety cap on iterations.}
\item{active_bgus}{See Description.}
}
\value{
Integer vector, valid relabelled solution.
}
\description{
Starting from an atomic (each BGU = its own cluster) or any other initial partition, iteratively merges the least-valid cluster into one of its best ranked neighbours until every cluster satisfies the validity criterion. Randomness is introduced at two points (matching the paper):   1. The focal cluster is sampled from the n_low lowest-ranked.   2. The merge target is sampled from the n_high highest-ranked      adjacent clusters
}
