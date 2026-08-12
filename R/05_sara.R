#' @param focal_cl Integer, focal cluster ID.
#' @param adj_cls  Integer vector, IDs of adjacent clusters.
#' @param sol      Integer vector.
#' @param W        Numeric matrix.
#' @param row_W    Numeric vector.
#' @param col_W    Numeric vector.
#' @return         Named numeric vector of scores (one per adjacent cluster).
.rank_adjacent_clusters <- function(focal_cl, adj_cls, sol, W, row_W, col_W) {
  focal_mbrs <- cluster_members(focal_cl, sol)
  
  vapply(adj_cls, function(cl) {
    other_mbrs <- cluster_members(cl, sol)
    w_fo <- sum(W[focal_mbrs, other_mbrs])
    w_of <- sum(W[other_mbrs, focal_mbrs])
    r_f <- sum(row_W[focal_mbrs])
    c_f <- sum(col_W[focal_mbrs])
    r_o <- sum(row_W[other_mbrs])
    c_o <- sum(col_W[other_mbrs])
    t1 <- if (r_f > 0 && c_o > 0)
      w_fo ^ 2 / (r_f * c_o)
    else
      0
    t2 <- if (r_o > 0 && c_f > 0)
      w_of ^ 2 / (r_o * c_f)
    else
      0
    t1 + t2
  }, numeric(1L))
}

#' Merge cl_from into cl_into and re_label consecutively
#'
#' @param sol      Integer vector.
#' @param cl_from  Integer, cluster that disappears.
#' @param cl_into  Integer, cluster that absorbs.
#' @return         Integer vector (relabelled 1...m-1)
.merge_clusters <- function(sol, cl_from, cl_into) {
  sol[sol == cl_from] <- cl_into
  getOrderedVec(sol)
}

#' SARA: Spatially Adaptive Random Aggregation  [Algorithm 2]
#'
#' Starting from an atomic (each BGU = its own cluster) or any
#' other initial partition, iteratively merges the least-valid
#' cluster into one of its best ranked neighbours until every
#' cluster satisfies the validity criterion.
#'
#' Randomness is introduced at two points (matching the paper):
#'   1. The focal cluster is sampled from the n_low lowest-ranked.
#'   2. The merge target is sampled from the n_high highest-ranked
#'      adjacent clusters
#'
#' @param init_sol  Integer vector, starting partition (e.g. 1:N).
#' @param W         Numeric matrix (N x N).
#' @param adj       Integer/logical matrix (N x N).
#' @param row_W     Numeric vector, rowSums(W).
#' @param col_W     Numeric vector, colSums(W).
#' @param SC_min    Numeric, minimum self-containment threshold.
#' @param Pop_min   Numeric, minimum population size threshold.
#' @param n_low     Integer, candidates to sample focal from (default 3).
#' @param n_high    Integer, candidates to sample target from (default 3).
#' @param max_iter  Integer, safety cap on iterations.
#' @return          Integer vector, valid relabelled solution.
sara <- function(init_sol,
                 W,
                 adj,
                 row_W,
                 col_W,
                 SC_min,
                 Pop_min,
                 Pop_max,
                 n_low = 3L,
                 n_high = 3L,
                 max_iter = 10000L,
                 active_bgus = NULL,
                 verbose = TRUE) {
  sol <- getOrderedVec(init_sol)
  iter <- 0L
  
  if (is.null(active_bgus)) {
    active_bgus <- seq_len(length(sol))
  }
  
  prev_n <- length(unique(sol[active_bgus]))
  stuck <- 0L
  
  repeat {
    iter <- iter + 1L
    if (iter > max_iter) {
      break
    }
    
    clusters <- unique(sol[active_bgus])
    
    curr_n <- length(clusters)
    
    if (curr_n == prev_n) {
      stuck <- stuck + 1L
    } else {
      stuck <- 0L
    }
    
    prev_n <- curr_n
    
    if (stuck > 100L) {
      if (verbose) message("SARA stuck")
      break
    }
    
    # Validity scores for every current cluster
    vs <- vapply(clusters, function(cl) {
      validity_score(cluster_members(cl, sol),
                     W,
                     row_W,
                     col_W,
                     adj,
                     SC_min,
                     Pop_min, Pop_max)
    }, numeric(1L))
    
    # All valid - done
    if (all(vs >= 1))
      break
    
    # Step 1: sample focal cluster from the n_low least-valid
    n_invalid <- max(1L, sum(vs < 1))
    low_idx <- order(vs)[seq_len(min(n_low, n_invalid))]
    focal_cl <- clusters[sample(low_idx, 1L)]
    focal_mbrs <- cluster_members(focal_cl, sol)
    
    # Step 2: Find adjacent BGUs and their cluster IDs, preferring other
    # active clusters first so the merge stays confined to the region
    # currently being re-partitioned. Falls back to any adjacent cluster
    # if the focal cluster has no active neighbour.
    adj_bgus <- setdiff(intersect(which(colSums(adj[focal_mbrs, , drop = FALSE]) > 0L), active_bgus), focal_mbrs)
    adj_cls <- setdiff(unique(sol[adj_bgus]), focal_cl)

    if (length(adj_cls) == 0L) {
      adj_bgus <- setdiff(which(colSums(adj[focal_mbrs, , drop = FALSE]) > 0L), focal_mbrs)
      adj_cls <- setdiff(unique(sol[adj_bgus]), focal_cl)
    }

    if (length(adj_cls) == 0L) {
      # No neighbours anywhere in the graph -- a genuinely disconnected
      # fragment. Nothing this function can do; stop gracefully.
      break
    }
    
    focal_pop <- compute_population(focal_mbrs, row_W)
    
    adj_pops <- vapply(adj_cls, function(cl) compute_population(cluster_members(cl, sol), row_W), numeric(1L))
    scores <- .rank_adjacent_clusters(focal_cl, adj_cls, sol, W, row_W, col_W)
    
    for (i in seq_along(adj_cls)) {
      
      adj_pop <- adj_pops[i]
      merged_pop <- focal_pop + adj_pop
      
      if (merged_pop > Pop_max) {
        scores[i] <- scores[i] * 0.5
        next
      }
      
    }
    
    top_idx <- order(scores)[seq_len(min(n_high, length(adj_cls)))]
    target_cl <- adj_cls[sample(top_idx, 1L)]
    
    # Step 4: merge focal into target, check contiguity
    candidate <- getOrderedVec(.merge_clusters(sol, focal_cl, target_cl))
    merged_id <- candidate[focal_mbrs[1L]]  # new label of merged cluster
    merged_mbrs <- cluster_members(merged_id, candidate)
    
    if (is_cluster_contiguous(merged_mbrs, adj)) {
      sol <- candidate
    }
    # If non-contiguous, skip this merge and try again next iteration
  }
  getOrderedVec(sol)
}

sara_max_pop_score <- function(sol, W, row_W, col_W, Pop_max) {

  lma_df <- build_lma_df(sol, W, row_W, col_W)
  
  n_violations <- sum(lma_df$pop > Pop_max)
  
  max_pop <- max(lma_df$pop)
  excess <- max(max_pop - Pop_max, 0)
  
  list(
    n_violations = n_violations,
    excess = excess
  )
}

is_better_max_pop <- function(score_a, score_b) {
  
  if (score_a$excess != score_b$excess) {
    return(score_a$excess < score_b$excess)
  }
  
  return(score_a$n_violations <
           score_b$n_violations)
}

sara_iterated <- function(W,
                          adj,
                          row_W,
                          col_W,
                          SC_min,
                          Pop_min,
                          Pop_max,
                          max_it = 5L,
                          verbose = TRUE,
                          temp = NULL) {

  start_sol <- getOrderedVec(seq_len(length(row_W)))
  iter <- 0L
  best_sol <- start_sol

  best_score <- list(
    n_violations = Inf,
    excess = Inf
  )

  # Tracks the true best solution across passes (best_sol below is
  # re-aliased to the current attempt every pass for the split logic).
  overall_best_sol <- start_sol
  overall_best_score <- best_score

  min_accept_slf <- 0.5
  Min_self <- SC_min
  init_min_self <- SC_min

  active_bgus <- start_sol
  prev_mbrs <- active_bgus

  while (Min_self > 0.05) {

    iter <- iter + 1L
    if (iter > max_it) {
      break
    }

    if (!is.null(temp)) {
      sol <- temp
      temp <- NULL
    }
    else {
      sol <- sara(
        init_sol = start_sol,
        W = W,
        adj = adj,
        row_W = row_W,
        col_W = col_W,
        SC_min = Min_self,
        Pop_min = Pop_min,
        Pop_max = Pop_max,
        n_low = 5L,
        n_high = 5L,
        max_iter = 10000L,
        active_bgus = active_bgus,
        verbose = verbose
      )
    }

    current_score <- sara_max_pop_score(sol, W, row_W, col_W, Pop_max)

    if (current_score$excess == 0) {
      return(sol)
    }

    if (is_better_max_pop(current_score, overall_best_score)) {
      overall_best_score <- current_score
      overall_best_sol <- sol
    }

    best_score <- current_score
    best_sol <- sol
    lma_df <- build_lma_df(sol, W, row_W, col_W)
    if (verbose) {
      print(current_score$excess)
      print(lma_df)
    }
    big_cls <- lma_df$cluster[lma_df$pop > Pop_max]
    split_bgus <- integer(0)
    cl <- big_cls[1]
    all_mbrs <- c()
    
    for (cl in big_cls) {
      mbrs <- cluster_members(cl, best_sol)
      all_mbrs <- c(all_mbrs, mbrs)
      split_bgus <- c(split_bgus, mbrs)
      new_ids <- seq(
        max(sol) + 1L,
        max(sol) + length(mbrs))
      sol[mbrs] <- new_ids
    }
    active_bgus <- split_bgus
    if (identical(sort(all_mbrs), sort(prev_mbrs))) {
      Min_self <- Min_self * 0.95
      if (verbose) print(Min_self)
    } else {
      Min_self <- init_min_self
    }
    
    prev_mbrs <- all_mbrs
    start_sol <- getOrderedVec(sol)

  }

  getOrderedVec(overall_best_sol)
}