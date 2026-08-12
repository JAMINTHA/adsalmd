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
