# ============================================================
# FILE: R/04_contiguity.R
# Contiguity checking (BFS) and validity score.
# Replaces the original checkClusValidity() / is_cluster_contiguous().
# ============================================================

#' Is a single market spatially contiguous?
#'
#' Tests whether a market is one connected piece on the map rather than two or
#' more disconnected islands sharing a label. Contiguity is a hard constraint
#' throughout the package: every operator checks it before returning a move,
#' and \code{\link{compute_objective}()} zeroes any partition that breaks it.
#'
#' @details
#' Runs a breadth-first search from the first member, expanding only along
#' adjacency links whose far end is \emph{also in the market}. The market is
#' contiguous if the BFS reaches every member.
#'
#' Each wave is computed with one \code{colSums()} over the frontier's rows of
#' \code{adj}, so the whole search is a handful of vectorised passes rather than
#' an element-by-element graph walk. It terminates when the frontier empties —
#' there is no iteration cap and no way for it to give a wrong answer by
#' running out of steps, which was the failure mode of the original
#' fixed-iteration implementation.
#'
#' Markets of 0 or 1 members are contiguous by definition and return
#' immediately.
#'
#' @section Assumptions:
#' \code{adj} is expected to be symmetric. The BFS treats a non-zero entry as a
#' link and only ever reads rows of the frontier, so with an asymmetric matrix
#' it silently tests reachability along directed edges instead — usually not
#' what you want. Non-zero (not necessarily \code{1}) counts as adjacent, so
#' logical, integer and shared-border-length matrices all work.
#'
#' @param members Integer vector of BGU indices in the market.
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @return \code{TRUE} if every member is reachable from every other within the
#'   market, \code{FALSE} otherwise.
#'
#' @family contiguity
#' @seealso \code{\link{all_clusters_contiguous}()} to test a whole partition,
#'   \code{\link{fix_contiguity}()} to repair one by splitting broken markets.
#'
#' @examples
#' ## Six BGUs in a line: 1-2-3-4-5-6
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#'
#' is_cluster_contiguous(c(1, 2, 3), adj)   # a connected run
#' is_cluster_contiguous(c(1, 3), adj)      # a gap at BGU 2
#' is_cluster_contiguous(integer(0), adj)   # trivially TRUE
#' is_cluster_contiguous(4, adj)            # singleton: TRUE
#'
#' @export
is_cluster_contiguous <- function(members, adj) {
  if (length(members) <= 1L) return(TRUE)

  visited  <- members[1L]
  frontier <- members[1L]

  while (length(frontier) > 0L) {
    # All BGUs adjacent to the current frontier AND inside the cluster
    candidates <- which(colSums(adj[frontier, , drop = FALSE]) > 0L)
    new_nodes  <- setdiff(intersect(candidates, members), visited)
    visited    <- c(visited, new_nodes)
    frontier   <- new_nodes
  }

  length(visited) == length(members)
}

#' Is every market in a partition contiguous?
#'
#' Applies \code{\link{is_cluster_contiguous}()} to each market and returns
#' \code{TRUE} only if all of them pass. This is the feasibility gate used by
#' the objective function and by every perturbation operator.
#'
#' @details
#' Short-circuits on the first failure, so cost depends on where the broken
#' market sits in label order — typically far less than a full pass over the
#' partition when a move has just broken something.
#'
#' The operators use it in a generate-and-test loop: build a candidate
#' partition, call this, and fall back to the unchanged solution if it returns
#' \code{FALSE}. That is why the search never has to reason about contiguity
#' explicitly — it simply never accepts a partition that fails here.
#'
#' @param sol Integer vector of market assignments (length N).
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @return \code{TRUE} if every market is contiguous, \code{FALSE} otherwise.
#'
#' @family contiguity
#' @seealso \code{\link{is_cluster_contiguous}()},
#'   \code{\link{compute_objective}()}, \code{\link{fix_contiguity}()}
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#'
#' all_clusters_contiguous(c(1, 1, 1, 2, 2, 2), adj)   # two connected runs
#' all_clusters_contiguous(c(1, 2, 1, 2, 2, 2), adj)   # market 1 is split
#' all_clusters_contiguous(1:6, adj)                   # all singletons
#'
#' ## Find which market is the problem
#' sol <- c(1, 2, 1, 2, 2, 2)
#' vapply(unique_clusters(sol),
#'        function(cl) is_cluster_contiguous(cluster_members(cl, sol), adj),
#'        logical(1))
#'
#' @export
all_clusters_contiguous <- function(sol, adj) {
  clusters <- unique_clusters(sol)
  for (cl in clusters) {
    if (!is_cluster_contiguous(cluster_members(cl, sol), adj)) {
      return(FALSE)
    }
  }
  TRUE
}

#' Validity score of a single market (Equation 1)
#'
#' Scores how close a market is to being a legitimate labour market area:
#' \code{1} when it satisfies every criterion, and a fraction below 1 measuring
#' how far short it falls. This is what \code{\link{sara}()} and
#' \code{\link{sha}()} rank markets by when deciding which one to merge away
#' next.
#'
#' @details
#' \deqn{VS(M_i) = \min\!\left(\frac{SC(M_i)}{SC_{min}}, 1\right)
#'                \cdot \min\!\left(\frac{Pop(M_i)}{Pop_{min}}, 1\right)
#'                \cdot C(M_i)}
#'
#' Each factor is capped at 1, so exceeding a threshold earns no credit and
#' cannot mask a shortfall on the other criterion — a market with twice the
#' required population but half the required self-containment scores 0.5, not
#' 1. The contiguity factor \eqn{C(M_i)} is 0 or 1, so a broken market scores 0
#' outright regardless of its other merits.
#'
#' This multiplicative capped form is what makes \code{VS >= 1} exactly
#' equivalent to "valid on every criterion", which is the stopping test in both
#' aggregation routines.
#'
#' \code{Pop_max} is accepted for signature compatibility with the callers but
#' does not enter the score; the population ceiling is handled separately, by
#' \code{\link{sara_iterated}()} splitting oversized markets and by
#' \code{\link{compute_objective}()} rejecting them.
#'
#' @param members Integer vector of BGU indices in the market.
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param SC_min Numeric in (0, 1], the minimum self-containment threshold.
#' @param Pop_min Numeric > 0, the minimum population threshold.
#' @param Pop_max Numeric, the maximum population threshold (default
#'   \code{Inf}). Accepted but not used (see Details).
#' @return Numeric in \code{[0, 1]}. Exactly \code{1} if and only if the market
#'   is contiguous and meets both thresholds.
#'
#' @family contiguity
#' @seealso \code{\link{sara}()} and \code{\link{sha}()}, which merge the
#'   lowest-scoring markets; \code{\link{compute_sc}()} and
#'   \code{\link{compute_population}()} for the two inputs.
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' ## A real market: valid on every criterion
#' validity_score(1:3, W, row_W, col_W, adj, SC_min = 0.55, Pop_min = 100)
#'
#' ## Population floor it cannot reach
#' validity_score(1:3, W, row_W, col_W, adj, SC_min = 0.55, Pop_min = 1000)
#'
#' ## Non-contiguous: zero regardless of anything else
#' validity_score(c(1, 3), W, row_W, col_W, adj, SC_min = 0.55, Pop_min = 100)
#'
#' ## Score every market in a partition
#' sol <- c(1, 1, 1, 2, 2, 2)
#' vapply(unique_clusters(sol), function(cl)
#'          validity_score(cluster_members(cl, sol), W, row_W, col_W, adj,
#'                         SC_min = 0.55, Pop_min = 100),
#'        numeric(1))
#'
#' @export
validity_score <- function(members, W, row_W, col_W, adj, SC_min, Pop_min, Pop_max = Inf) {
  sc_val  <- compute_sc(members, W, row_W, col_W)$sc
  pop_val <- compute_population(members, row_W)
  contig  <- as.numeric(is_cluster_contiguous(members, adj))

  pop_score <- min(pop_val / Pop_min, 1)
  sc_score <- min(sc_val / SC_min, 1)
  sc_score * pop_score * contig
}
