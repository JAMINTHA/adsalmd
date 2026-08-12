# ============================================================
# FILE: R/12_validation_utils.R
# Helpers for validating the method against known-answer test cases,
# and for choosing SC_min/Pop_min/Pop_max/Pop_tar from a reference
# partition rather than guessing.
# ============================================================

#' Check whether two partitions are identical up to relabelling
#'
#' Two cluster-assignment vectors represent the same grouping if every
#' BGU pair that's co-clustered in one is co-clustered in the other, and
#' vice versa -- regardless of which integer label each cluster happens
#' to have. Use this to check whether an algorithm's output
#' (e.g. \code{result$best_sol}) recovered a known true partition, since
#' cluster labels are arbitrary and won't line up by value.
#'
#' @param sol_a Integer vector, e.g. the true/known partition.
#' @param sol_b Integer vector, e.g. the algorithm's output partition.
#'   Must be the same length as sol_a.
#' @return Named list: \code{match} (logical, TRUE if identical up to
#'   relabelling), \code{table} (the cluster-by-cluster contingency
#'   table -- a clean diagonal-only pattern is what a perfect match
#'   looks like; anything off that diagonal shows exactly which BGUs
#'   / clusters disagree).
#' @export
partitions_match <- function(sol_a, sol_b) {
  stopifnot(length(sol_a) == length(sol_b))
  tab <- table(true = sol_a, predicted = sol_b)
  is_match <- all(rowSums(tab > 0) == 1L) && all(colSums(tab > 0) == 1L)
  list(match = is_match, table = tab)
}

#' Suggest SC_min / Pop_min / Pop_max / Pop_tar from a reference partition
#'
#' If you already know a "good" partition for some data -- a synthetic
#' test case with a known true answer, a previous method's output, or a
#' business rule -- this reads off that partition's own per-market
#' self-containment and population (via \code{\link{build_lma_df}}) and
#' proposes thresholds a small margin below/above the extremes, so real
#' markets in that partition aren't excluded by the very thresholds
#' meant to describe them.
#'
#' This is a starting point, not a substitute for the kind of grid
#' search described in the paper (Table 2) -- treat the output as
#' something to sanity-check and refine, not use blindly.
#'
#' @param sol   Integer vector, a reference partition.
#' @param W     Numeric OD matrix (N x N).
#' @param row_W Numeric vector, rowSums(W).
#' @param col_W Numeric vector, colSums(W).
#' @param margin Numeric in [0, 1), how far below the reference
#'   partition's minimum SC/population (and above its maximum
#'   population) to set the suggested thresholds (default 0.10, i.e.
#'   10%).
#' @return Named list: \code{SC_min}, \code{Pop_min}, \code{Pop_max},
#'   \code{Pop_tar} (the median market population), plus \code{lma_df}
#'   (the per-market table the suggestion was derived from, so you can
#'   see exactly which market drove each number).
#' @export
suggest_thresholds <- function(sol, W, row_W, col_W, margin = 0.10) {
  lma_df <- build_lma_df(sol, W, row_W, col_W)

  list(
    SC_min  = max(0, min(lma_df$sc)  * (1 - margin)),
    Pop_min = max(0, min(lma_df$pop) * (1 - margin)),
    Pop_max = max(lma_df$pop) * (1 + margin),
    Pop_tar = stats::median(lma_df$pop),
    lma_df  = lma_df
  )
}
