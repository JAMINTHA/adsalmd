\name{is_cluster_contiguous}
\alias{is_cluster_contiguous}
\title{Check if a single cluster is spatially contiguous using BFS}
\usage{
is_cluster_contiguous(members, adj)
}
\arguments{
\item{members}{Integer vector, BGU indices in the cluster.}
\item{adj}{Integer/logical matrix (N x N), adjacency matrix.}
}
\value{
Logical.
}
\description{
A cluster is contiguous if all its member BGUs are reachable from any one member via adjacency links that stay within the cluster. BFS is O(|members|) and avoids the fragile iteration limit in the original code.
}
