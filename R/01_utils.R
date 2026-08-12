# ============================================================
# FILE: R/01_utils.R
# Pure utility functions — no side effects, no globals.
# ============================================================

#' Relabel cluster IDs to consecutive integers 1, 2, ..., m
#'
#' Replaces your original get_ordered_vec() / assign_custom_order().
#' Uses match() for O(N) performance instead of a loop.
#'
#' @param sol Integer vector of cluster assignments (length N).
#' @return    Integer vector with labels 1 … m.
get_ordered_vec <- function(sol) {
  unique_vals <- sort(unique(sol))
  match(sol, unique_vals)
}

#' Alias for get_ordered_vec()
#'
#' Every caller in this package (sara(), sara_iterated(), adsa_almd(),
#' the operators, sha()) uses the camelCase name; this alias keeps them
#' working without having to touch every call site.
#'
#' @param sol Integer vector of cluster assignments.
#' @return Integer vector with labels 1 ... m.
#' @export
getOrderedVec <- get_ordered_vec

#' Return indices of BGUs belonging to a given cluster
#'
#' @param cluster_id Integer, cluster label.
#' @param sol        Integer vector of assignments.
#' @return Integer vector of indices.
cluster_members <- function(cluster_id, sol) {
  which(sol == cluster_id)
}

#' Sorted unique cluster labels present in a solution
#'
#' @param sol Integer vector.
#' @return Sorted integer vector.
unique_clusters <- function(sol) {
  sort(unique(as.integer(sol)))
}

#' Number of clusters (markets) in a solution
#'
#' @param sol Integer vector.
#' @return Integer scalar m.
n_clusters <- function(sol) {
  length(unique(sol))
}
