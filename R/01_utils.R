# ============================================================
# FILE: R/01_utils.R
# Pure utility functions — no side effects, no globals.
# ============================================================

#' Relabel cluster IDs to consecutive integers 1, 2, ..., m
#'
#' Canonicalises a partition so that its cluster labels are exactly
#' \code{1:m} in order of first appearance of the sorted unique labels. Every
#' function in the package that creates or destroys clusters relabels through
#' this, so that \code{max(sol)} is always the number of markets and a fresh
#' label can be minted with \code{max(sol) + 1L}.
#'
#' @details
#' Implemented as \code{match(sol, sort(unique(sol)))}, which is O(N log N) and
#' vectorised, replacing the per-element loop of the original
#' \code{assign_custom_order()}. The mapping is order-preserving: the smallest
#' original label becomes 1, the next becomes 2, and so on. Gaps are closed and
#' duplicates are untouched, so the \emph{grouping} is never changed — only the
#' names of the groups.
#'
#' Because labels are arbitrary, do not compare two partitions with
#' \code{identical()} even after relabelling if they were built independently;
#' use \code{\link{partitions_match}()}.
#'
#' @param sol Integer vector of cluster assignments (length N). May contain
#'   gaps, negative values or any other integer labels.
#' @return Integer vector of the same length, with labels \code{1 … m} where
#'   \code{m = length(unique(sol))}.
#'
#' @family partition utilities
#' @seealso \code{\link{n_clusters}()}, \code{\link{unique_clusters}()},
#'   \code{\link{partitions_match}()}
#'
#' @examples
#' get_ordered_vec(c(7, 7, 3, 9, 3))
#'
#' ## Gaps left by a merge are closed
#' sol <- c(1, 1, 2, 2, 3, 3)
#' sol[sol == 2] <- 1          # merge market 2 into market 1
#' sol                          # labels 1 and 3 — gap at 2
#' get_ordered_vec(sol)         # labels 1 and 2
#'
#' ## The grouping itself is unchanged
#' identical(
#'   table(get_ordered_vec(c(7, 7, 3, 9, 3))),
#'   unname(table(c(7, 7, 3, 9, 3)))
#' )
#'
#' @export
get_ordered_vec <- function(sol) {
  unique_vals <- sort(unique(sol))
  match(sol, unique_vals)
}

#' Alias for get_ordered_vec()
#'
#' A camelCase copy of \code{\link{get_ordered_vec}()}. Every internal caller
#' (\code{\link{sara}()}, \code{\link{sara_iterated}()},
#' \code{\link{adsa_almd}()}, the operators and \code{\link{sha}()}) uses this
#' name; it is kept so those call sites — and any of your own scripts written
#' against the original code — keep working.
#'
#' @details
#' This is a plain copy of the function object made at build time, not a
#' wrapper, so there is no call overhead and the two names cannot drift apart.
#' New code should prefer \code{\link{get_ordered_vec}()}.
#'
#' @param sol Integer vector of cluster assignments.
#' @return Integer vector with labels \code{1 … m}.
#'
#' @family partition utilities
#' @seealso \code{\link{get_ordered_vec}()} for the full description.
#'
#' @examples
#' getOrderedVec(c(5, 5, 2, 8))
#' identical(getOrderedVec(c(5, 5, 2, 8)), get_ordered_vec(c(5, 5, 2, 8)))
#'
#' @export
getOrderedVec <- get_ordered_vec

#' Return indices of BGUs belonging to a given cluster
#'
#' The package's standard way of going from a market label to its member BGUs.
#' Nearly every metric takes a \code{members} vector produced by this function.
#'
#' @details
#' A thin wrapper on \code{which(sol == cluster_id)}. It returns an empty
#' integer vector — not an error — when \code{cluster_id} is absent from
#' \code{sol}, which is what lets callers loop over stale label lists safely.
#'
#' @param cluster_id Integer, the market label to look up.
#' @param sol Integer vector of assignments (length N).
#' @return Integer vector of BGU indices, in increasing order. Length 0 if no
#'   BGU carries that label.
#'
#' @family partition utilities
#' @seealso \code{\link{unique_clusters}()} to enumerate the labels to pass in.
#'
#' @examples
#' sol <- c(1, 1, 1, 2, 2, 2)
#' cluster_members(1, sol)
#' cluster_members(2, sol)
#' cluster_members(3, sol)   # absent: integer(0)
#'
#' ## Typical use: loop over every market
#' for (cl in unique_clusters(sol)) {
#'   cat("market", cl, "->", cluster_members(cl, sol), "\n")
#' }
#'
#' @export
cluster_members <- function(cluster_id, sol) {
  which(sol == cluster_id)
}

#' Sorted unique cluster labels present in a solution
#'
#' Enumerates the markets in a partition. Coerces to integer first, so a
#' solution stored as doubles (for example after \code{read.csv()}) still
#' yields integer labels that can be compared and used to mint new ones.
#'
#' @param sol Integer vector of assignments.
#' @return Sorted integer vector of the distinct labels in \code{sol}.
#'
#' @family partition utilities
#' @seealso \code{\link{cluster_members}()}, \code{\link{n_clusters}()}
#'
#' @examples
#' unique_clusters(c(3, 1, 1, 2, 3))
#' unique_clusters(c(2.0, 1.0, 1.0))   # coerced to integer
#'
#' @export
unique_clusters <- function(sol) {
  sort(unique(as.integer(sol)))
}

#' Number of clusters (markets) in a solution
#'
#' The \code{m} that multiplies GCI in the base fitness
#' (\code{\link{compute_base_fitness}()}) and that is reported in every run
#' summary.
#'
#' @details
#' Counts distinct labels, so it is correct whether or not \code{sol} has been
#' relabelled by \code{\link{get_ordered_vec}()}. After relabelling it equals
#' \code{max(sol)}, which is why the operators can mint a fresh market with
#' \code{max(sol) + 1L}.
#'
#' @param sol Integer vector of assignments.
#' @return Integer scalar \code{m}, the number of markets.
#'
#' @family partition utilities
#' @seealso \code{\link{compute_base_fitness}()}, which uses this as its
#'   multiplier.
#'
#' @examples
#' n_clusters(c(1, 1, 2, 2, 3))
#' n_clusters(c(7, 7, 3))                       # labels need not be 1:m
#' n_clusters(get_ordered_vec(c(7, 7, 3)))      # same answer after relabelling
#'
#' @export
n_clusters <- function(sol) {
  length(unique(sol))
}
