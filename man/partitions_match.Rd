\name{partitions_match}
\alias{partitions_match}
\title{Check whether two partitions are identical up to relabelling}
\usage{
partitions_match(sol_a, sol_b)
}
\arguments{
\item{sol_a}{Integer vector, e.g. the true/known partition.}
\item{sol_b}{Integer vector, e.g. the algorithm's output partition.   Must be the same length as sol_a.}
}
\value{
Named list: \\code\{match\} (logical, TRUE if identical up to   relabelling), \\code\{table\} (the cluster-by-cluster contingency   table -- a clean diagonal-only pattern is what a perfect match   looks like; anything off that diagonal shows exactly which BGUs   / clusters disagree).
}
\description{
Two cluster-assignment vectors represent the same grouping if every BGU pair that's co-clustered in one is co-clustered in the other, and vice versa -- regardless of which integer label each cluster happens to have. Use this to check whether an algorithm's output (e.g. \\code\{result$best_sol\}) recovered a known true partition, since cluster labels are arbitrary and won't line up by value.
}
