\name{get_ordered_vec}
\alias{get_ordered_vec}
\alias{getOrderedVec}
\title{Relabel cluster IDs to consecutive integers 1, 2, ..., m}
\usage{
get_ordered_vec(sol)
}
\arguments{
\item{sol}{Integer vector of cluster assignments (length N).}
}
\value{
Integer vector with labels 1 … m.
}
\description{
Replaces your original get_ordered_vec() / assign_custom_order(). Uses match() for O(N) performance instead of a loop.
}
