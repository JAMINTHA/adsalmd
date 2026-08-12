# ============================================================
# FILE: R/06_operators.R
# Ten group-based perturbation operators.
# Implements the operators from Martínez-Bernabeu et al. (2012)
# as adapted for AdSA-ALMD.
#
# All operators:
#   - Accept the full parameter set as arguments (no globals).
#   - Return the original sol unchanged if the move is infeasible.
#   - Never break contiguity of any existing cluster.
# ============================================================

# ---- Internal helpers (not exported) -------------------------

#' Border BGUs: members of a cluster that touch another cluster
.border_bgus <- function(members, sol, adj) {
  all_adj      <- which(colSums(adj[members, , drop = FALSE]) > 0L)
  external_adj <- setdiff(all_adj, members)
  if (length(external_adj) == 0L) 
    return(integer(0))
  members[vapply(members, function(g) 
    any(adj[g, external_adj] > 0L), logical(1L))]
}

#' Cluster IDs adjacent to a given set of members
.adjacent_clusters <- function(members, sol, adj) {
  adj_bgus <- setdiff(which(colSums(adj[members, , drop = FALSE]) > 0L), members)
  unique(sol[adj_bgus])
}

#' Pairwise cluster interaction score (mutual CI)
.cluster_interaction <- function(mbrs_a, mbrs_b, W, row_W, col_W) {
  w_ab <- sum(W[mbrs_a, mbrs_b])
  w_ba <- sum(W[mbrs_b, mbrs_a])
  ra   <- sum(row_W[mbrs_a])
  ca   <- sum(col_W[mbrs_a])
  rb   <- sum(row_W[mbrs_b])
  cb   <- sum(col_W[mbrs_b])
  t1   <- if (rb > 0 && cb > 0) 
    w_ab^2 / (ra * cb)
  else 0
  t2   <- if (rb > 0 && ca > 0) 
    w_ba^2 / (rb * ca) 
  else 0
  result <- t1 + t2
  
  if(!is.finite(result)) return(0)
  result
}



#' Apply a candidate solution only if all clusters remain contiguous
.accept_if_contiguous <- function(candidate, sol, adj) {
  if (all_clusters_contiguous(candidate, adj))
    candidate 
  else 
    sol
}

.tournament_select <- function(values,
                               n_way = 3L,
                               minimize = TRUE){
  n <- length(values)
  if(n == 0L)
    return(NULL)
  n_way <- min(n_way, n)
  if(minimize){
    idx <- order(values)[seq_len(n_way)]
  } else{
    idx <- order(values, decreasing = TRUE) [seq_len(n_way)]
  }
  
  sample(idx, 1L)
}

# ---- Operator 1: Merge least-valid cluster into best neighbour ------
#
# Greedy merge: the market with lowest SC is absorbed by whichever
# adjacent market has the strongest mutual interaction with it.

operator_1 <- function(sol, W, adj, row_W, col_W, ...) {
  lma_df   <- build_lma_df(sol, W, row_W, col_W)
  if(nrow(lma_df) < 2L)
    return(sol)
  
  idx <- .tournament_select(lma_df$sc, n_way = 3L, minimize = TRUE)

  worst_cl <- lma_df$cluster[idx]
  
  w_mbrs   <- cluster_members(worst_cl, sol)
  
  adj_cls  <- .adjacent_clusters(w_mbrs, sol, adj)
  
  if (length(adj_cls) == 0L) 
    return(sol)

  scores  <- vapply(adj_cls, function(cl) {
    .cluster_interaction(w_mbrs, cluster_members(cl, sol), W, row_W, col_W)
  }, numeric(1L))
  
  scores <- pmax(scores, 1e-10)
  
  best_cl <- adj_cls[sample(length(adj_cls),
                            size = 1L,
                            prob = scores / sum(scores))]
  
  cand <- sol
  cand[sol == best_cl] <- worst_cl
  cand <- getOrderedVec(cand)
  
  merged_id <- cand[w_mbrs[1L]]
  
  merged_mbrs <- cluster_members(merged_id, cand)
  
  for(i in seq_along(merged_mbrs)){
    new_id <- max(cand) + 1L
    cand[merged_mbrs[i]] <- new_id
  }
  
  cand <- getOrderedVec(cand)
  
  cand <- sha(
    bgus = merged_mbrs,
    sol = cand,
    W = W,
    adj = adj,
    row_W = row_W,
    col_W = col_W,
    SC_min = params$SC_min,
    Pop_min = params$Pop_min,
    max_iter = 1000
  )
  
  if(all_clusters_contiguous(cand, adj)){
    return(cand)
  }

sol
}

# ---- Operator 2: Expand smallest cluster by absorbing one border BGU ----
#
# The adjacent BGU with the highest CI contribution to the focal cluster
# is annexed.

operator_2 <- function(sol, W, adj, row_W, col_W, ...) {
  lma_df   <- build_lma_df(sol, W, row_W, col_W)
  if (nrow(lma_df) < 2L)
    return(sol)
  
  idx <- .tournament_select(lma_df$pop, n_way = 3L, minimize = TRUE)
  
  small_cl <- lma_df$cluster[idx]
  members  <- cluster_members(small_cl, sol)
  adj_bgus <- setdiff(which(colSums(adj[members, , drop = FALSE]) > 0L), members)
  if (length(adj_bgus) == 0L) 
    return(sol)

  scores <- vapply(adj_bgus, function(g) {
    compute_ci(g, members, W, row_W, col_W)
  }, numeric(1L))
  
  b_idx <- .tournament_select(scores, n_way = 3L, minimize = FALSE)
  best_g <- adj_bgus[b_idx]

  candidate <- getOrderedVec(replace(sol, best_g, small_cl))
  .accept_if_contiguous(candidate, sol, adj)
}



# ---- Operator 3: Dismantle largest cluster — reassign BGUs to neighbours ----
#
# Each member of the biggest cluster is moved to whichever adjacent
# cluster gives it the best CI.

operator_3 <- function(sol, W, adj, row_W, col_W, ...) {
  lma_df  <- build_lma_df(sol, W, row_W, col_W)
  if(nrow(lma_df) < 2L)
    return(sol)
  
  idx <- .tournament_select(lma_df$n_bgus, n_way = 3L, minimize = FALSE)
  
  big_cl  <- lma_df$cluster[idx]
  members <- cluster_members(big_cl, sol)
  if (length(members) <= 1L) 
    return(sol)

  borders <- .border_bgus(members, sol, adj)
  if(length(borders) == 0L)
    return(sol)
  
  for(g in sample(borders)){
    adj_cls <- setdiff(unique(sol[which(adj[g, ] > 0L)]), big_cl)
    if(length(adj_cls) == 0L)
      next
    
    scores <- vapply(adj_cls, function(cl)
      compute_ci(g, cluster_members(cl, sol), W, row_W, col_W), numeric(1L))
    t_idx <- .tournament_select(scores, n_way = 3L, minimize = FALSE)
    best_cl <- adj_cls[t_idx]
    
    cand <- getOrderedVec(replace(sol, g, best_cl))
    
    if(all_clusters_contiguous(cand, adj)) {
      return(cand)
    }
      
  }
  
  sol
}







# ---- Operator 4: Remove a border BGU from its cluster, assign to best neighbour ----
#
# Targets the cluster with the lowest SC.  A randomly chosen border BGU
# from that cluster is reassigned to its best adjacent cluster.

operator_4 <- function(sol, W, adj, row_W, col_W, ...) {
  lma_df   <- build_lma_df(sol, W, row_W, col_W)
  if(nrow(lma_df) < 2L)
    return(sol)
  
  idx <- .tournament_select(lma_df$sc, n_way = 3L, minimize = TRUE)
  focal_cl <- lma_df$cluster[idx]
  members  <- cluster_members(focal_cl, sol)
  borders  <- .border_bgus(members, sol, adj)
  if (length(borders) == 0L) 
    return(sol)
  
  b_scores <- vapply(borders, function(g)
    compute_ci(g, setdiff(members, g), W, row_W, col_W), numeric(1L))
  b_idx <- .tournament_select(b_scores, n_way = 3L, minimize = TRUE)

  g       <- borders[b_idx]
  adj_cls <- setdiff(unique(sol[which(adj[g, ] > 0L)]), focal_cl)
  if (length(adj_cls) == 0L) 
    return(sol)

  scores  <- vapply(adj_cls, function(cl) {
    compute_ci(g, cluster_members(cl, sol), W, row_W, col_W)
  }, numeric(1L))
  t_idx <- .tournament_select(scores, n_way = 3L, minimize = FALSE)
  best_cl <- adj_cls[which.max(scores)]

  candidate <- getOrderedVec(replace(sol, g, best_cl))
  .accept_if_contiguous(candidate, sol, adj)
}

# ---- Operator 5: Include a border BGU from an adjacent cluster ----
#
# The adjacent BGU that would contribute most CI to the focal cluster is
# pulled in. (Inverse of Operator 4.)

operator_5 <- function(sol, W, adj, row_W, col_W, ...) {
  lma_df   <- build_lma_df(sol, W, row_W, col_W)
  if(nrow(lma_df) < 2L)
    return(sol)
  
  idx <- .tournament_select(lma_df$sc, n_way = 3L, minimize = TRUE)
  focal_cl <- lma_df$cluster[idx]
  members  <- cluster_members(focal_cl, sol)
  adj_bgus <- setdiff(which(colSums(adj[members, , drop = FALSE]) > 0L), members)
  if (length(adj_bgus) == 0L) 
    return(sol)

  scores <- vapply(adj_bgus, function(g) {
    compute_ci(g, members, W, row_W, col_W)
  }, numeric(1L))
  
  b_idx <- .tournament_select(scores, n_way = 3L, minimize = FALSE)
  best_g <- adj_bgus[b_idx]

  candidate <- getOrderedVec(replace(sol, best_g, focal_cl))
  .accept_if_contiguous(candidate, sol, adj)
}

# ---- Operator 6: Segregation — split a contiguous subset into a new cluster ----
#
# A random starting BGU within the largest cluster seeds a BFS expansion
# of up to n_split BGUs, which are detached as a new cluster.

operator_6 <- function(sol, W, adj, row_W, col_W, ...) {
  lma_df  <- build_lma_df(sol, W, row_W, col_W)
  if(nrow(lma_df) < 2L)
    return(sol)
  
  idx <- .tournament_select(lma_df$n_bgus, n_way = 3L, minimize = FALSE)
  big_cl  <- lma_df$cluster[idx]
  members <- cluster_members(big_cl, sol)
  if (length(members) < 4L) 
    return(sol)

  ci_scores <- vapply(members, function(g){
    M_g <- setdiff(members, g)
    compute_ci(g, M_g, W, row_W, col_W)
  }, numeric(1L))

  for(n_split in c(1L, 2L, 3L)) {
    if(n_split >= length(members))
      next
    subset_g <- members[order(ci_scores)[seq_len(n_split)]]
    remainder <- setdiff(members, subset_g)

    if(is_cluster_contiguous(subset_g, adj) &&
       is_cluster_contiguous(remainder, adj)){
      new_cl <- max(sol) + 1L
      cand <- getOrderedVec(replace(sol, subset_g, new_cl))
      return(cand)
    }
  }
  sol

}

# ---- Operator 7: Annexation — absorb a weakly connected adjacent cluster ----
#
# A random cluster selects the adjacent cluster with the lowest mutual
# interaction and absorbs it entirely.

operator_7 <- function(sol, W, adj, row_W, col_W, ...) {
  clusters <- unique_clusters(sol)
  if (length(clusters) < 2L) 
    return(sol)

  fcl_idx <- sample(length(clusters), 1L)
  focal_cl  <- fcl_idx
  focal_mbrs <- cluster_members(focal_cl, sol)
  adj_cls   <- .adjacent_clusters(focal_mbrs, sol, adj)
  if (length(adj_cls) == 0L) 
    return(sol)

  scores   <- vapply(adj_cls, function(cl) {
    .cluster_interaction(focal_mbrs, cluster_members(cl, sol), W, row_W, col_W)
  }, numeric(1L))
  weak_cl  <- adj_cls[which.min(scores)]

  candidate <- getOrderedVec(replace(sol, cluster_members(weak_cl, sol), focal_cl))
  .accept_if_contiguous(candidate, sol, adj)
}

# ---- Operator 8: Random border reassignment ----
#
# A random cluster is chosen; one of its border BGUs is moved to a
# random adjacent cluster.

operator_8 <- function(sol, W, adj, row_W, col_W, ...) {
  focal_cl <- sample(unique_clusters(sol), 1L)
  members  <- cluster_members(focal_cl, sol)
  borders  <- .border_bgus(members, sol, adj)
  if (length(borders) == 0L) 
    return(sol)

  g       <- sample(borders, 1L)
  adj_cls <- setdiff(unique(sol[which(adj[g, ] > 0L)]), focal_cl)
  if (length(adj_cls) == 0L) 
    return(sol)

  target_cl <- sample(adj_cls, 1L)
  candidate <- getOrderedVec(replace(sol, g, target_cl))
  .accept_if_contiguous(candidate, sol, adj)
}

# ---- Operator 9: Exchange border BGUs between two adjacent clusters ----
#
# One BGU from cluster A that touches cluster B is swapped with one
# BGU from cluster B that touches cluster A.

operator_9 <- function(sol, W, adj, row_W, col_W, ...) {
  clusters <- unique_clusters(sol)
  if (length(clusters) < 2L) 
    return(sol)

  cl1      <- sample(clusters, 1L)
  adj_cls  <- .adjacent_clusters(cluster_members(cl1, sol), sol, adj)
  if (length(adj_cls) == 0L) 
    return(sol)
  cl2      <- sample(adj_cls, 1L)

  mbrs1    <- cluster_members(cl1, sol)
  mbrs2    <- cluster_members(cl2, sol)

  b1_to_2 <- .border_bgus(mbrs1, sol, adj)
  b1_to_2 <- b1_to_2[vapply(b1_to_2, function(g) 
    any(adj[g, mbrs2] > 0L), logical(1L))]

  b2_to_1 <- .border_bgus(mbrs2, sol, adj)
  b2_to_1 <- b2_to_1[vapply(b2_to_1, function(g) 
    any(adj[g, mbrs1] > 0L), logical(1L))]

  if (length(b1_to_2) == 0L || length(b2_to_1) == 0L) 
    return(sol)

  g1 <- sample(b1_to_2, 1L)
  g2 <- sample(b2_to_1, 1L)

  candidate        <- sol
  candidate[g1]    <- cl2
  candidate[g2]    <- cl1
  candidate        <- getOrderedVec(candidate)
  .accept_if_contiguous(candidate, sol, adj)
}

# ---- Operator 10: Transfer a small subset between adjacent clusters ----
#
# One to three border BGUs from one cluster that touch a neighbouring
# cluster are moved to that neighbour as a group.

# Helper function for operator 10

fix_contiguity <- function(sol, adj, max_iter = 5000L){
  clusters <- unique_clusters(sol)
  
  for(cl in clusters){
    members <- cluster_members(cl, sol)
    if(length(members) <= 1L) next
    if(is_cluster_contiguous(members, adj)) next
    
    remaining <- members
    first <- TRUE
    iter <- 0L
    
    while (length(remaining) > 0L) {
      iter <- iter + 1L
      if(iter > max_iter){
        break
      }
      seed <- remaining[1L]
      visited <- seed
      frontier <- seed
      
      while(length(frontier) > 0L){
        candidates <- which(colSums(adj[frontier, , drop = FALSE]) > 0L)
        new_nodes <- setdiff(intersect(candidates, remaining), visited)
        visited <- c(visited, new_nodes)
        frontier <- new_nodes
      }
      
      if(!first){
        new_cl <- max(sol) + 1L
        sol[visited] <- new_cl
      }
      
      remaining <- setdiff(remaining, visited)
      first <- FALSE
    }
  }
  getOrderedVec(sol)
}

operator_10 <- function(sol, W, adj, row_W, col_W, best_sol = NULL, ...){
  
  low_time_con = TRUE
  
  if(is.null(best_sol) || identical(sol, best_sol)){
    # simple transfer
    clusters <- unique_clusters(sol)
    if(length(clusters) < 2L) return(sol)
    cl1 <- sample(clusters, 1L)
    adj_cls <- .adjacent_clusters(cluster_members(cl1, sol), sol, adj)
    if(length(adj_cls) == 0L) return(sol)
    cl2 <- sample(adj_cls, 1L)
    borders <- .border_bgus(cluster_members(cl1, sol), sol, adj)
    adj_b <- borders[vapply(borders, function(g) any(adj[g, cluster_members(cl2, sol)] > 0L), logical(1L))]
    if(length(adj_b) == 0L) return(sol)
    move_g <- sample(adj_b, sample(1L:min(3L, length(adj_b)), 1L))
    cand <- getOrderedVec(replace(sol, move_g, cl2))
    
    return(.accept_if_contiguous(cand, sol, adj))
  }
  
  # Crossover with sols
  
  N <- length(sol)
  n_cls_2 <- n_clusters(best_sol)
  
  r_min <- max(1L, round(n_cls_2*0.01))
  r_max <- max(1L, round(n_cls_2 * 0.06))
  
  r <- sample(r_min:r_max, 1L)
  
  cls_2 <- unique_clusters(best_sol)
  
  selected <- sample(cls_2, min(r, length(cls_2)))
  
  child <- sol
  
  offset <- max(child) + 1L
  
  for(cl in selected){
    mbrs <- cluster_members(cl, best_sol)
    child[mbrs] <- cl + offset
  }
  
  child <- getOrderedVec(child)
  
  child <- fix_contiguity(child, adj)
  
  all_bgus <- seq_len(N)
  
  affected_bgus <- all_bgus
  
  max_iter <- 3000L

  if(low_time_con && runif(1L) > 0.3){
    max_iter <- 1000L
    affected_bgus <- unique(c(unlist(lapply(selected, function(cl) cluster_members(cl, best_sol))),
                              unlist(lapply(selected, function(cl){
                                mbrs <- cluster_members(cl, best_sol)
                                setdiff(which(colSums(adj[mbrs,,drop = FALSE]) > 0L), mbrs)
                              }))))
  }
  
  child <- sha(
    bgus = affected_bgus,
    sol = child,
    W = W,
    adj = adj,
    row_W = row_W,
    col_W = col_W,
    SC_min = params$SC_min,
    Pop_min = params$Pop_min,
    max_iter = max_iter
  )

  if(all_clusters_contiguous(child, adj)){
    return(child)
  }
  sol
}


operator_11 <- function(sol, W, adj, row_W, col_W, Pop_max = Inf, ...){
  
  lma_df <- build_lma_df(sol, W, row_W, col_W)
  
  big_cls <- lma_df$cluster[lma_df$pop > Pop_max]
  
  if(length(big_cls) == 0L) return(sol)
  
  big_pops <- lma_df$pop[match(big_cls, lma_df$cluster)]
  idx <- .tournament_select(-big_pops, n_way = 3L, minimize  = TRUE)
  
  big_cl <- big_cls[idx]
  big_mbrs <- cluster_members(big_cl, sol)
  
  if(length(big_mbrs) < 3L) return(sol)
  
  cand <- sol
  for(i in seq_along(big_mbrs)){
    new_id <- max(cand) + 1L
    cand[big_mbrs[i]] <- new_id
  }
  
  cand <- getOrderedVec(cand)
  
  cand <- sha(
    bgus = big_mbrs, sol = cand, W = W, adj = adj, row_W = row_W, col_W = col_W, SC_min = params$SC_min,
    Pop_min = params$Pop_min, max_iter = 500L)
  
  if(all_clusters_contiguous(cand, adj)){
    return(cand)
  }
  
  sol
}

init_op_tracker <- function(){
  list(
    called = rep(0L, 11L),
    accepted = rep(0L, 11L),
    improved = rep(0L, 11L),
    unchanged = rep(0L, 11L)
  )
}

update_op_tracker <- function(tracker,
                              op_idx,
                              sol_before,
                              sol_after,
                              F_before,
                              F_after,
                              accepted,
                              improved_best){
  tracker$called[op_idx] <- tracker$called[op_idx] + 1L
  
  if(identical(sol_before, sol_after)){
    tracker$unchanged[op_idx] <- tracker$unchanged[op_idx] + 1L
    
    return(tracker)
  }
  
  if(accepted){
    tracker$accepted[op_idx] <- tracker$accepted[op_idx] + 1L
  }
  
  if(improved_best){
    tracker$improved[op_idx] <- tracker$improved[op_idx] + 1L
  }
  
  tracker
}


print_op_tracker <- function(tracker){
  cat("\n=== Operator Performance Report ===\n")
  cat(
    sprintf(
      "%-12s %8s %8s %8s %8s %8s %8s \n",
      "Operator",
      "Called",
      "Unchanged",
      "Changed",
      "Accepted",
      "Improved",
      "Acc%"
    )
  )
  cat(paste(rep("-", 68), collapse = ""), "\n")
  
  for(i in 1:11){
    called <- tracker$called[i]
    unchanged <- tracker$unchanged[i]
    changed <- called - unchanged
    accepted <- tracker$accepted[i]
    improved <- tracker$improved[i]
    acc_pct <- if (called > 0)
      round(100 * accepted / called, 1)
    else
      0
    
    cat(
      sprintf(
        "%-12s %8s %8s %8s %8s %8s %8s \n",
        paste0("operator_", i),
        called, 
        unchanged,
        changed,
        accepted,
        improved,
        acc_pct
      )
    )
  }
  cat(paste(rep("-", 68),collapse = ""), "\n")
  cat(
    sprintf(
      "%-12s %8s %8s %8s %8s %8s %8s \n",
      "TOTAL",
      sum(tracker$called),
      sum(tracker$unchanged),
      sum(tracker$called) - sum(tracker$unchanged),
      sum(tracker$accepted),
      sum(tracker$improved),
      round(100*sum(tracker$accepted) / max(sum(tracker$called), 1), 1)
    )
  )
  cat("======================================\n\n")
}
# ---- Dispatcher ---------------------------------------------------

#' Apply one randomly chosen operator (uniform over 1–10)
#'
#' @param sol     Integer vector.
#' @param W       Numeric matrix.
#' @param adj     Integer/logical matrix.
#' @param row_W   Numeric vector.
#' @param col_W   Numeric vector.
#' @param SC_min  Numeric.
#' @param Pop_min Numeric.
#' @return Perturbed integer vector (may equal sol if move was infeasible).
apply_random_operator <- function(sol, W, adj, row_W, col_W, best_sol = NULL, ...) {
  op_fns <- list(
    operator_1, 
    operator_2, 
    operator_3, 
    operator_4, 
    operator_5,
    operator_6, 
    operator_7, 
    operator_8, 
    operator_9, 
    operator_10,
    operator_11
  )
  
  weights <- c(1,1,1,1,1,2,1,1,1,0.5,1.5)
  
  op_idx <- sample(11L, 1L, prob = weights / sum(weights))
  
  if(op_idx == 10L){
    new_sol <- operator_10(sol, W, adj, row_W, col_W, best_sol = best_sol, ...)
  }else if(op_idx == 11L){
    new_sol <- operator_11(sol, W, adj, row_W, col_W, Pop_max = params$Pop_max, ...)
  }else {
    new_sol <- op_fns[[op_idx]](sol, W, adj, row_W, col_W, ...)
  }
  list(sol = new_sol, op_idx = op_idx)
}
