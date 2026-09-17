# ============================================================
# FILE: R/09_sha.R
# SHA — Subset Hierarchical Aggregation: SARA's merge logic applied
# to a named subset of BGUs, leaving the rest of the map untouched.
# ============================================================

#' SHA: re-aggregate a subset of BGUs into valid markets
#'
#' Rebuilds the partition of a \emph{region} rather than the whole map. Given a
#' set of BGUs — typically just shattered into singletons by an operator — it
#' merges them into valid markets by the same interaction-weighted logic as
#' \code{\link{sara}()}, then writes the result back into the surrounding
#' partition. This is the repair engine behind
#' \code{\link{operator_1}()}, \code{\link{operator_10}()} and
#' \code{\link{operator_11}()}.
#'
#' @details
#' The function works in a \strong{local index space}: \code{bgus} is
#' renumbered \code{1 … length(bgus)} and all merging happens there, with
#' \code{global_idx()} translating back whenever a metric needs real BGU
#' indices. That is what keeps the cost proportional to the subset rather than
#' to \code{N}, and what lets the surrounding markets stay exactly as they were.
#'
#' Each iteration:
#' \enumerate{
#'   \item scores every local cluster with \code{\link{validity_score}()},
#'     evaluated on the \emph{global} members so self-containment is measured
#'     against the whole commuting matrix, not just the subset;
#'   \item stops when every local cluster scores \code{>= 1};
#'   \item samples a focal cluster uniformly from the invalid ones;
#'   \item finds its neighbours inside the subset and merges it into one,
#'     sampled with probability proportional to \code{.cluster_interaction()}
#'     (floored at \code{1e-10}).
#' }
#'
#' @section BGUs with no neighbour left in the subset:
#' This is the case that distinguishes SHA from \code{\link{sara}()}, which can
#' only ever merge within the map it was given. When a focal cluster has no
#' neighbour remaining inside \code{bgus}:
#' \itemize{
#'   \item if it has neighbours \emph{outside} the subset, it is folded
#'     directly into one of those existing markets — chosen with the same
#'     interaction-weighted sampling — and dropped from further local
#'     processing. A shattered region can therefore partly dissolve into the
#'     markets around it rather than being forced to re-form on its own.
#'   \item if it has no neighbours anywhere on the map, it is simply dropped
#'     from the local problem so the rest of the subset can still resolve.
#'     Those BGUs keep whatever label they end up with at the write-back step.
#' }
#'
#' @section Write-back:
#' Whatever local clusters remain when the loop ends are written into
#' \code{result_sol} with fresh labels \code{max(sol) + cl}, guaranteeing no
#' collision with the untouched markets. The whole vector is then relabelled to
#' \code{1 … m}.
#'
#' @section Termination and guarantees:
#' The loop ends on success, on \code{max_iter}, or when the subset empties.
#' Reaching \code{max_iter} returns a perfectly usable partition that simply is
#' not fully valid, so \strong{a result is not a guarantee that every market
#' meets the thresholds}. Nor is contiguity guaranteed: merges here are not
#' contiguity-checked, which is why every caller tests
#' \code{\link{all_clusters_contiguous}()} on the result and discards it if it
#' fails.
#'
#' @param bgus Integer vector of BGU indices to re-aggregate. Usually the
#'   members of one market that has just been broken into singletons.
#' @param sol Integer vector of market assignments for the whole map (length
#'   N), already carrying the shattered labels.
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param adj Integer/logical \code{N x N} adjacency matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @param SC_min Numeric in (0, 1], the minimum self-containment threshold.
#' @param Pop_min Numeric > 0, the minimum market population.
#' @param max_iter Integer, cap on merge iterations (default 6000). Callers use
#'   much lower values — 500 in \code{\link{operator_11}()}, 1000 in
#'   \code{\link{operator_1}()} — to bound the cost per perturbation.
#' @return Integer vector of length N, relabelled to \code{1 … m}: the input
#'   partition with the \code{bgus} region re-aggregated.
#'
#' @family SARA
#' @seealso \code{\link{sara}()} for the whole-map equivalent;
#'   \code{\link{validity_score}()} for the stopping criterion;
#'   \code{\link{operator_1}()}, \code{\link{operator_10}()} and
#'   \code{\link{operator_11}()} for the callers.
#'
#' @examples
#' ## Toy system: six BGUs in a line, forming two 3-BGU markets
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' ## Market 1 has been shattered into singletons; market 2 is intact
#' shattered <- c(3, 4, 5, 2, 2, 2)
#'
#' set.seed(51)
#' out <- sha(bgus = 1:3, sol = shattered, W = W, adj = adj,
#'            row_W = row_W, col_W = col_W,
#'            SC_min = 0.55, Pop_min = 200, max_iter = 500)
#' out
#' build_lma_df(out, W, row_W, col_W)
#'
#' ## The untouched region keeps its grouping
#' identical(out[4:6], rep(out[4], 3))
#'
#' ## Callers must check contiguity themselves -- SHA does not
#' all_clusters_contiguous(out, adj)
#'
#' @export
sha <- function(bgus,
                sol,
                W,
                adj,
                row_W,
                col_W,
                SC_min,
                Pop_min,
                max_iter = 6000L){
  local_sol <- seq_along(bgus)
  names(local_sol) <- as.character(bgus)
  global_idx <- function(local_mbrs){
    bgus[local_mbrs]
  }

  result_sol <- sol
  iter <- 0L

  repeat{
    iter <- iter + 1L
    if(iter > max_iter){
      break
    }

    if(length(bgus) == 0L)
      break

    local_cls <- unique_clusters(local_sol)

    vs <- vapply(local_cls, function(cl){
      local_mbrs <- which(local_sol == cl)
      global_mbrs <- global_idx(local_mbrs)
      validity_score(global_mbrs, W, row_W, col_W, adj, SC_min, Pop_min)
    }, numeric(1L))

    if(all(vs >= 1))
      break

    invalid_cls <- local_cls[vs < 1]

    focal_local <- sample(invalid_cls, 1L)
    focal_local_mbrs <- which(local_sol == focal_local)
    focal_global_mbrs <- global_idx(focal_local_mbrs)

    adj_global_bgus <- setdiff(which(colSums(adj[focal_global_mbrs, ,
                                                 drop = FALSE]) > 0L), focal_global_mbrs)

    adj_subset_bgus <- intersect(adj_global_bgus, bgus)
    if(length(adj_subset_bgus) == 0L) {
      # No neighbour left inside the local `bgus` subset. If it has real
      # neighbours elsewhere on the map, fold it directly into whichever
      # of those existing (outside) clusters it interacts with most.
      if(length(adj_global_bgus) == 0L) {
        # Genuinely isolated on the full map -- drop it from further
        # local processing so the rest of the subset can still resolve.
        keep_idx <- setdiff(seq_along(bgus), focal_local_mbrs)
        bgus <- bgus[keep_idx]
        local_sol <- local_sol[keep_idx]
        next
      }

      outside_cls <- unique(result_sol[adj_global_bgus])
      out_scores <- vapply(outside_cls, function(cl){
        .cluster_interaction(focal_global_mbrs, cluster_members(cl, result_sol), W, row_W, col_W)
      }, numeric(1L))
      out_scores <- pmax(out_scores, 1e-10)

      target_cl <- outside_cls[sample(length(outside_cls),
                                      size = 1L,
                                      prob = out_scores / sum(out_scores))]

      result_sol[focal_global_mbrs] <- target_cl

      keep_idx <- setdiff(seq_along(bgus), focal_local_mbrs)
      bgus <- bgus[keep_idx]
      local_sol <- local_sol[keep_idx]
      next
    }
    adj_local_cls <- setdiff(unique(local_sol[match(adj_subset_bgus, bgus)]), focal_local)

    if(length(adj_local_cls) == 0L)
      next

    scores <- vapply(adj_local_cls, function(cl){
      other_local_mbrs <- which(local_sol == cl)
      other_global_mbrs <- global_idx(other_local_mbrs)

      .cluster_interaction(focal_global_mbrs, other_global_mbrs, W, row_W, col_W)
    }, numeric(1L))

    scores <- pmax(scores, 1e-10)

    target_local <- adj_local_cls[sample(length(adj_local_cls),
                                         size = 1L,
                                         prob = scores / sum(scores))]

    local_sol[local_sol == focal_local] <- target_local
    local_sol <- getOrderedVec(local_sol)
  }

  max_global <- max(sol)

  if(length(bgus) > 0L) {
    local_cls <- unique_clusters(local_sol)

    for(cl in local_cls){
      local_mbrs <- which(local_sol == cl)
      global_mbrs <- global_idx(local_mbrs)
      new_id <- max_global + cl
      result_sol[global_mbrs] <- new_id
    }
  }

  getOrderedVec(result_sol)
}
