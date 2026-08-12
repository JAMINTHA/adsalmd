# ============================================================
# FILE: R/03_objective.R
# Penalty functions and the AdSA-ALMD objective function.
# All formulas match the paper exactly.
# ============================================================

# ---- Individual penalty terms -------------------------------

#' Self-containment penalty for one market   [Equation 13]
#'
#' PSC(Ms) = 1 - (SC(Ms) / SC_min)^r   if SC(Ms) < SC_min
#'         = 0                           otherwise
#'
#' @param sc     Numeric, self-containment of market Ms.
#' @param SC_min Numeric, minimum SC threshold.
#' @param r      Integer, steepness exponent (typically 3 or 4).
#' @return Numeric scalar >= 0.
penalty_sc_single <- function(sc, SC_min, r) {
  if (sc >= SC_min) return(0)
  1 - (sc / SC_min)^r
}

penalty_pop_min_single <- function(pop, Pop_min, r){
  if(pop >= Pop_min) return(0)
  1-(pop/Pop_min) ^ r
}

penalty_pop_max_single <- function(pop, Pop_max, r){
  if(pop <= Pop_max*1.5) return(0)
  base_penalty <- (pop/Pop_max - 1)^r
  
  return(base_penalty)
}

penalty_pop_max_total <- function(pop_vec, Pop_max, r){
  sum(vapply(pop_vec, penalty_pop_max_single, numeric(1L), Pop_max = Pop_max, r=r))
}

penalty_concentration <- function(pop_vec, weight = 1.0){
  total <- sum(pop_vec)
  if(total == 0) return(0)
  
  shares <- pop_vec / total
  
  hhi <- sum(shares^2)
  
  m <- length(pop_vec)
  
  perfect <- 1/m
  
  excess <- max(hhi - perfect, 0)
  
  weight * excess
}
  
#' Population penalty for one market   [Equation 14]
#'
#' PPop(Ms) = 1 - (Pop(Ms) / Pop_min)^r   if Pop(Ms) < Pop_min
#'          = 0                             otherwise
#'
#' #' @param pop     Numeric, population of market Ms.
#' #' @param Pop_min Numeric, minimum population threshold.
#' #' @param r       Integer, steepness exponent.
#' #' @return Numeric scalar >= 0.
#' penalty_pop_single <- function(pop, Pop_min, r) {
#'   if (pop >= Pop_min) return(0)
#'   1 - (pop / Pop_min)^r
#' }

# ---- Aggregate penalties ------------------------------------

#' Total self-containment penalty   [Equation 13 summed]
#'
#' PSC(x) = sum_{Ms in M} PSC(Ms)
#'
#' @param sc_vec Numeric vector, SC values for all markets.
#' @param SC_min Numeric.
#' @param r      Integer.
#' @return Numeric scalar.
penalty_sc_total <- function(sc_vec, SC_min, r) {
  sum(vapply(sc_vec, penalty_sc_single, numeric(1L),
             SC_min = SC_min, r = r))
}

#' Total population penalty   [Equation 14 summed]
#'
#' PPop(x) = sum_{Ms in M} PPop(Ms)
#'
#' @param pop Numeric vector, populations for all markets.
#' @param Pop_min Numeric.
#' @param r       Integer.
#' @return Numeric scalar.
penalty_pop_min_total <- function(pop_vec, Pop_min, r = 3L) {
  sum(vapply(pop_vec, penalty_pop_min_single, numeric(1L),
             Pop_min = Pop_min, r = r))
}

#' Total penalty P(x) = PSC(x) + PPop(x)   [Equation 9]
#'
#' @param lma_df  data.frame from build_lma_df().
#' @param SC_min  Numeric.
#' @param Pop_min Numeric.
#' @param r       Integer.
#' @return Numeric scalar.
penalty_total <- function(lma_df, SC_min, Pop_min, Pop_max = Inf, Pop_tar = NULL, r) {
  conc_weight <- 2.0
  p_sc <- penalty_sc_total(lma_df$sc, SC_min, r)
  p_conc <- penalty_concentration(lma_df$pop, weight = conc_weight)
  
  p_sc + p_conc
}

# ---- Objective function -------------------------------------

#' AdSA-ALMD objective function F(x)   [Equation 5]
#'
#' F(x) = (GCI(x) - P(x)) * prod_{i=1}^{m} C(Mi)
#'
#' The contiguity product makes F = 0 whenever any market is
#' non-contiguous, acting as a hard feasibility constraint.
#'
#' NOTE: all_clusters_contiguous() is defined in 04_contiguity.R.
#'
#' @param sol               Integer vector of assignments (length N).
#' @param W                 Numeric matrix (N x N).
#' @param row_W             Numeric vector (length N).
#' @param col_W             Numeric vector (length N).
#' @param SC_min            Numeric.
#' @param Pop_min           Numeric.
#' @param r                 Integer, penalty exponent.
#' @param adj               Integer/logical matrix (N x N).
#' @param check_contiguity  Logical; set FALSE only for T0 estimation.
#' @return Numeric scalar (0 if any cluster is non-contiguous).
compute_objective <- function(sol, W, row_W, col_W,
                               SC_min, Pop_min, Pop_max = Inf, Pop_tar = NULL, r, adj,
                               check_contiguity = TRUE) {
  if (check_contiguity && !all_clusters_contiguous(sol, adj)) {
    return(0)
  }
  lma_df <- build_lma_df(sol, W, row_W, col_W)
  if(is.finite(Pop_max) && any(lma_df$pop > Pop_max)){
    return(0)
  }
  if(any(lma_df$pop < Pop_min)){
    return(0)
  }

  gci    <- compute_gci(sol, W, row_W, col_W)
 
  p      <- penalty_total(lma_df, SC_min, Pop_min, Pop_max, Pop_tar, r)

  gci - p
}
