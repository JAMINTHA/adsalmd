# ============================================================
# FILE: R/02_lma_metrics.R
# LMA-specific metrics.  Every formula is cited to the paper.
#
# Notation (matching paper):
#   W      : OD matrix (N x N).  W[i, j] = commuters from i to j.
#   row_W  : rowSums(W) — precomputed once and passed in.
#   col_W  : colSums(W) — precomputed once and passed in.
#   G      : full BGU index set 1:N.
#   Ms     : set of BGU indices in market s.
#   g      : single BGU index.
#   M\g    : all BGUs in g's market, excluding g.
# ============================================================

# ---- Cohesion -----------------------------------------------

#' Cohesion index of a single BGU within its market (Equation 6)
#'
#' Measures how strongly one BGU \code{g} is tied to the rest of its own
#' market, in both directions: the share of \code{g}'s outgoing commuters that
#' stay inside the market, and the share of the market's commuters that flow
#' into \code{g}. Summing this over every BGU gives the global cohesion index
#' that the whole method maximises.
#'
#' @details
#' \deqn{CI(g) = \frac{W(g, M \setminus g)^2}{W(g, G)\, W(G, M \setminus g)}
#'             + \frac{W(M \setminus g, g)^2}{W(M \setminus g, G)\, W(G, g)}}
#'
#' The first term is supply-side (where \code{g}'s residents work), the second
#' demand-side (where \code{g}'s workers live). Each is a squared flow
#' normalised by the product of the two marginals it could have gone to, so the
#' index is scale-free: doubling every flow in \code{W} leaves it unchanged.
#'
#' Edge cases, all of which occur routinely mid-search:
#' \itemize{
#'   \item a singleton market (\code{M_g} empty) returns 0 — a BGU alone has no
#'     internal cohesion to measure;
#'   \item a term whose denominator is 0 (a BGU with no outgoing or no incoming
#'     commuters at all) contributes 0 rather than \code{NaN}, so the two terms
#'     are handled independently.
#' }
#'
#' \code{M_g} must \emph{exclude} \code{g} itself. Passing the full member
#' vector inflates the result by counting \code{g}'s internal flow
#' \code{W[g, g]}. Callers normally write \code{setdiff(members, g)}.
#'
#' @param g Integer, index of the BGU to score.
#' @param M_g Integer vector, the other BGUs in \code{g}'s market
#'   (\eqn{M \setminus \{g\}}). May be empty.
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @return Numeric scalar \eqn{\ge 0}. Bounded above by 2 (one per direction).
#'
#' @family LMA metrics
#' @seealso \code{\link{compute_gci}()}, which sums this over all BGUs;
#'   \code{\link{is_bgu_misallocated}()}, which compares it across candidate
#'   markets.
#'
#' @examples
#' ## Toy system: six BGUs in a line, forming two 3-BGU markets
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#' sol <- c(1, 1, 1, 2, 2, 2)
#'
#' ## BGU 1 inside its own market
#' compute_ci(1, c(2, 3), W, row_W, col_W)
#'
#' ## The same BGU judged against the far market: much weaker
#' compute_ci(1, c(4, 5, 6), W, row_W, col_W)
#'
#' ## A singleton market has no internal cohesion
#' compute_ci(1, integer(0), W, row_W, col_W)
#'
#' @export
compute_ci <- function(g, M_g, W, row_W, col_W) {
  if (length(M_g) == 0L) return(0)

  W_g_Mg <- sum(W[g,   M_g])        # W(g,  M\g)
  W_Mg_g <- sum(W[M_g, g  ])        # W(M\g, g)
  W_g_G  <- row_W[g]                # W(g,  G)
  W_G_g  <- col_W[g]                # W(G,  g)
  W_Mg_G <- sum(row_W[M_g])         # W(M\g, G)
  W_G_Mg <- sum(col_W[M_g])         # W(G,  M\g)

  term1  <- if (W_g_G > 0 && W_G_Mg > 0) W_g_Mg^2 / (W_g_G  * W_G_Mg) else 0
  term2  <- if (W_Mg_G > 0 && W_G_g > 0) W_Mg_g^2 / (W_Mg_G * W_G_g ) else 0

  term1 + term2
}

#' Global cohesion index of a partition (Equation 7)
#'
#' The quality measure at the heart of the method: the total cohesion of a
#' whole partition, obtained by summing \code{\link{compute_ci}()} over every
#' BGU. \code{\link{compute_objective}()} maximises this minus penalties, and
#' \code{\link{compute_base_fitness}()} scales it by the number of markets for
#' cross-method comparison.
#'
#' @details
#' \deqn{GCI(x) = \sum_{g \in G} CI(g)}
#'
#' Higher is better. The index rewards partitions whose markets absorb their
#' own commuting, but on its own it is trivially maximised by one giant market
#' — which is why the objective function combines it with self-containment,
#' population and contiguity constraints.
#'
#' @section Backends:
#' Two implementations compute the identical value and this function dispatches
#' between them on problem size:
#' \describe{
#'   \item{\code{N < 50}}{the literal element-wise definition above: loop over
#'     markets, loop over members, call \code{\link{compute_ci}()} once each.
#'     Kept because it reads exactly like Equation 7, which matters when
#'     checking the implementation by hand.}
#'   \item{\code{N >= 50}}{delegates to \code{\link{compute_gci_fast}()}, which
#'     extracts the sub-matrix \code{W[Ms, Ms]} once per market and derives
#'     every per-BGU quantity from its row and column sums. Around 15–20x
#'     faster at \code{N = 520}.}
#' }
#' The threshold is on \code{length(sol)}, not on the number of markets. You
#' can call either backend directly to bypass the dispatch — useful when
#' verifying that the two agree on your own data.
#'
#' @param sol Integer vector of market assignments (length N).
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @return Numeric scalar, the global cohesion index \eqn{\ge 0}.
#'
#' @family LMA metrics
#' @seealso \code{\link{compute_gci_fast}()} for the large-N backend,
#'   \code{\link{compute_ci}()} for the per-BGU term,
#'   \code{\link{compute_base_fitness}()} for the comparison measure.
#'
#' @examples
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' ## The natural two-market partition
#' compute_gci(c(1, 1, 1, 2, 2, 2), W, row_W, col_W)
#'
#' ## A partition that splits the markets in the wrong place scores lower
#' compute_gci(c(1, 1, 2, 2, 2, 1), W, row_W, col_W)
#'
#' ## Every BGU alone: no internal cohesion anywhere
#' compute_gci(1:6, W, row_W, col_W)
#'
#' ## Both backends agree; only speed differs
#' all.equal(
#'   compute_gci(c(1, 1, 1, 2, 2, 2), W, row_W, col_W),
#'   compute_gci_fast(c(1, 1, 1, 2, 2, 2), W, row_W, col_W)
#' )
#'
#' @export
compute_gci <- function(sol, W, row_W, col_W) {
  # For N >= 50 use the fast sub-matrix version (~15-20x speedup)
  if (length(sol) >= 50L) return(compute_gci_fast(sol, W, row_W, col_W))

  clusters <- unique_clusters(sol)
  gci <- 0
  for (cl in clusters) {
    members <- cluster_members(cl, sol)
    for (g in members) {
      M_g <- setdiff(members, g)
      gci <- gci + compute_ci(g, M_g, W, row_W, col_W)
    }
  }
  gci
}

# ---- Self-containment ---------------------------------------

#' Supply-side self-containment of a market (Equation 10)
#'
#' The share of a market's resident workers who also work inside it. This is
#' the "do people who live here work here?" half of self-containment.
#'
#' @details
#' \deqn{SCSS(M_s) = \frac{W(M_s, M_s)}{W(M_s, G)}}
#'
#' The denominator is the market's total outgoing commuting, i.e. its working
#' population. A market with no resident workers at all returns 0 rather than
#' \code{NaN}.
#'
#' @param members Integer vector, the BGU indices in the market.
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @return Numeric in \code{[0, 1]}.
#'
#' @family LMA metrics
#' @seealso \code{\link{compute_scds}()} for the demand side,
#'   \code{\link{compute_sc}()} for both at once.
#'
#' @examples
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W)
#'
#' compute_scss(1:3, W, row_W)     # a real market: high
#' compute_scss(3:4, W, row_W)     # straddling the divide: low
#' compute_scss(1, W, row_W)       # a single BGU on its own
#'
#' @export
compute_scss <- function(members, W, row_W) {
  W_Ms_G <- sum(row_W[members])
  if (W_Ms_G == 0) return(0)
  sum(W[members, members]) / W_Ms_G
}

#' Demand-side self-containment of a market (Equation 11)
#'
#' The share of a market's jobs that are filled by its own residents — the
#' mirror image of \code{\link{compute_scss}()}.
#'
#' @details
#' \deqn{SCDS(M_s) = \frac{W(M_s, M_s)}{W(G, M_s)}}
#'
#' The denominator is the market's total incoming commuting, i.e. its
#' employment. A market with no jobs returns 0 rather than \code{NaN}.
#'
#' The two sides differ whenever a market is a net commuting source or sink: a
#' dormitory suburb can fill all of its few local jobs (high SCDS) while
#' sending most of its residents elsewhere to work (low SCSS). Taking the
#' minimum of the two, as \code{\link{compute_sc}()} does, is what stops such a
#' market from passing the validity test on one side alone.
#'
#' @param members Integer vector, the BGU indices in the market.
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @return Numeric in \code{[0, 1]}.
#'
#' @family LMA metrics
#' @seealso \code{\link{compute_scss}()}, \code{\link{compute_sc}()}
#'
#' @examples
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' col_W <- colSums(W)
#'
#' compute_scds(1:3, W, col_W)
#' compute_scds(3:4, W, col_W)
#'
#' @export
compute_scds <- function(members, W, col_W) {
  W_G_Ms <- sum(col_W[members])
  if (W_G_Ms == 0) return(0)
  sum(W[members, members]) / W_G_Ms
}

#' Self-containment of a market, both sides (Equation 12)
#'
#' Computes supply-side and demand-side self-containment and returns them
#' together with their minimum — the single \code{SC} value used by the
#' validity test, the penalty function and every operator that targets "the
#' least self-contained market".
#'
#' @details
#' \deqn{SC(M_s) = \min\bigl(SCSS(M_s),\ SCDS(M_s)\bigr)}
#'
#' Taking the minimum is deliberately strict: a market counts as
#' self-contained only if it retains both its workers and its jobs. See
#' \code{\link{compute_scds}()} for why the two sides can diverge.
#'
#' All three components are returned because the diagnostics table reports
#' \code{scss} and \code{scds} separately (paper Table 3) while the search only
#' ever uses \code{sc}.
#'
#' @param members Integer vector, the BGU indices in the market.
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @return Named list of three numerics in \code{[0, 1]}:
#'   \describe{
#'     \item{\code{scss}}{supply-side self-containment.}
#'     \item{\code{scds}}{demand-side self-containment.}
#'     \item{\code{sc}}{the minimum of the two.}
#'   }
#'
#' @family LMA metrics
#' @seealso \code{\link{validity_score}()} and
#'   \code{\link{penalty_sc_single}()}, the two consumers of \code{sc};
#'   \code{\link{build_lma_df}()} to get all markets at once.
#'
#' @examples
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' compute_sc(1:3, W, row_W, col_W)
#'
#' ## Does this market clear a 0.55 floor?
#' compute_sc(1:3, W, row_W, col_W)$sc >= 0.55
#' compute_sc(3:4, W, row_W, col_W)$sc >= 0.55
#'
#' @export
compute_sc <- function(members, W, row_W, col_W) {
  scss <- compute_scss(members, W, row_W)
  scds <- compute_scds(members, W, col_W)
  list(scss = scss, scds = scds, sc = min(scss, scds))
}

# ---- Population ---------------------------------------------

#' Population of a market
#'
#' The market's working population, taken as its total outgoing commuting
#' \eqn{Pop(M_s) = W(M_s, G)}. Every BGU resident who works anywhere is counted
#' exactly once, so summing over all markets recovers \code{sum(W)}.
#'
#' @details
#' Because population is read off the OD matrix rather than supplied
#' separately, it is automatically consistent with the self-containment
#' figures: there is no way for a market's population and its commuting
#' marginals to disagree.
#'
#' Note this is a \emph{sample} count when \code{W} is a sample (the Queensland
#' study uses a 4 percent sample), which is why \code{Pop_min} and
#' \code{Pop_max} in \code{\link{AdSA_params}} are on that same scale rather
#' than in whole persons.
#'
#' @param members Integer vector, the BGU indices in the market.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @return Numeric scalar \eqn{\ge 0}.
#'
#' @family LMA metrics
#' @seealso \code{\link{penalty_pop_min_single}()},
#'   \code{\link{penalty_pop_max_single}()}, \code{\link{validity_score}()}
#'
#' @examples
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W)
#'
#' compute_population(1:3, row_W)
#' compute_population(1:6, row_W) == sum(W)   # all BGUs: the whole sample
#'
#' @export
compute_population <- function(members, row_W) {
  sum(row_W[members])
}

# ---- Summary data frame -------------------------------------

#' Per-market summary table for a partition
#'
#' Turns a partition into one row per market with its population and
#' self-containment figures. This is the package's standard description of "how
#' good is each market in this solution", and the table you should inspect
#' first when a run produces an unexpected result.
#'
#' @details
#' For each label in \code{\link{unique_clusters}(sol)} the function collects
#' the market's members and evaluates \code{\link{compute_population}()} and
#' \code{\link{compute_sc}()} on them. Rows come out in sorted label order.
#' Markets are never dropped: a singleton BGU appears with
#' \code{n_bgus = 1} and whatever \code{sc} it happens to have.
#'
#' Cost is one pass per market, dominated by the \code{W[members, members]}
#' sub-matrix extractions inside \code{\link{compute_sc}()}. The annealing loop
#' calls this once per candidate solution via
#' \code{\link{compute_objective}()}, so it is on the hot path — build it once
#' and pass the result around rather than recomputing.
#'
#' @param sol Integer vector of market assignments (length N).
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @return A \code{data.frame} with one row per market and columns:
#'   \describe{
#'     \item{\code{cluster}}{integer, the market label.}
#'     \item{\code{pop}}{numeric, market population
#'       (\code{\link{compute_population}()}).}
#'     \item{\code{scss}}{numeric, supply-side self-containment.}
#'     \item{\code{scds}}{numeric, demand-side self-containment.}
#'     \item{\code{sc}}{numeric, \code{min(scss, scds)}.}
#'     \item{\code{n_bgus}}{integer, number of BGUs in the market.}
#'   }
#'   Character columns are never created (\code{stringsAsFactors = FALSE}).
#'
#' @family LMA metrics
#' @seealso \code{\link{compute_objective}()} and
#'   \code{\link{penalty_total}()}, which consume this table;
#'   \code{\link{suggest_thresholds}()}, which reads thresholds off it;
#'   \code{\link{build_results_table}()} for a cross-method summary.
#'
#' @examples
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' build_lma_df(c(1, 1, 1, 2, 2, 2), W, row_W, col_W)
#'
#' ## Which markets fail a 0.55 self-containment floor?
#' lma <- build_lma_df(c(1, 1, 2, 2, 3, 3), W, row_W, col_W)
#' lma[lma$sc < 0.55, ]
#'
#' ## Population always sums to the whole sample, whatever the partition
#' sum(build_lma_df(1:6, W, row_W, col_W)$pop) == sum(W)
#'
#' @export
build_lma_df <- function(sol, W, row_W, col_W) {
  clusters <- unique_clusters(sol)

  rows <- lapply(clusters, function(cl) {
    members <- cluster_members(cl, sol)
    sc_vals <- compute_sc(members, W, row_W, col_W)
    data.frame(
      cluster = cl,
      pop     = compute_population(members, row_W),
      scss    = sc_vals$scss,
      scds    = sc_vals$scds,
      sc      = sc_vals$sc,
      n_bgus  = length(members),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

# ---- Optimised GCI for large N ------------------------------

#' Global cohesion index, sub-matrix backend for large N (Equation 7)
#'
#' Computes exactly the same quantity as \code{\link{compute_gci}()} but with a
#' per-market rather than a per-BGU memory access pattern. This is the backend
#' \code{\link{compute_gci}()} dispatches to when \code{length(sol) >= 50}; call
#' it directly to force the fast path, or to cross-check the two.
#'
#' @details
#' The naive form indexes \code{W[g, M_g]} once per BGU — 520 separate row
#' extractions at \code{N = 520}, each allocating a fresh vector. This version
#' extracts \code{W[Ms, Ms]} once per market and reads every per-BGU quantity
#' off that sub-matrix:
#' \itemize{
#'   \item \eqn{W(g, M \setminus g)} is the row sum of the sub-matrix (the
#'     diagonal term \code{W[g, g]} cancels out of both numerator and the
#'     marginals);
#'   \item \eqn{W(M \setminus g, g)} is the corresponding column sum;
#'   \item \eqn{W(M \setminus g, G)} and \eqn{W(G, M \setminus g)} come from the
#'     precomputed market totals minus \code{row_W[g]} / \code{col_W[g]}, with
#'     no matrix access at all.
#' }
#' Around 15–20x faster than the element-wise version at \code{N = 520}.
#' Singleton markets are skipped outright (their cohesion is 0 by definition).
#'
#' The same zero-denominator guards as \code{\link{compute_ci}()} apply, so
#' results agree term by term and not just in total.
#'
#' @param sol Integer vector of market assignments (length N).
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @return Numeric scalar, the global cohesion index \eqn{\ge 0}. Identical to
#'   \code{\link{compute_gci}()} up to floating-point summation order.
#'
#' @family LMA metrics
#' @seealso \code{\link{compute_gci}()} for the dispatching front end and the
#'   description of both backends.
#'
#' @examples
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#' sol <- c(1, 1, 1, 2, 2, 2)
#'
#' compute_gci_fast(sol, W, row_W, col_W)
#'
#' ## Same answer as the element-wise backend
#' all.equal(
#'   compute_gci_fast(sol, W, row_W, col_W),
#'   compute_gci(sol, W, row_W, col_W)
#' )
#'
#' \donttest{
#' ## Where the speed difference shows up
#' set.seed(1)
#' N <- 200
#' Wb <- matrix(rpois(N * N, 1), N, N)
#' rb <- rowSums(Wb); cb <- colSums(Wb)
#' big <- rep(1:20, each = 10)
#' system.time(compute_gci_fast(big, Wb, rb, cb))
#' }
#'
#' @export
compute_gci_fast <- function(sol, W, row_W, col_W) {
  clusters <- unique_clusters(sol)
  gci      <- 0

  for (cl in clusters) {
    members <- cluster_members(cl, sol)
    n_m     <- length(members)
    if (n_m <= 1L) next

    # --- Extract sub-matrix once per cluster ---
    W_cl   <- W[members, members, drop = FALSE]
    rw_cl  <- rowSums(W_cl)   # W(g, M\g) for each g in order
    cw_cl  <- colSums(W_cl)   # W(M\g, g) for each g in order

    # Cluster-level totals (precomputed)
    cl_row <- sum(row_W[members])   # W(Ms, G)
    cl_col <- sum(col_W[members])   # W(G, Ms)

    for (j in seq_len(n_m)) {
      g       <- members[j]
      W_g_Mg  <- rw_cl[j]                  # W(g,  M\g)
      W_Mg_g  <- cw_cl[j]                  # W(M\g, g)
      W_g_G   <- row_W[g]                  # W(g,  G)
      W_G_g   <- col_W[g]                  # W(G,  g)
      W_Mg_G  <- cl_row - row_W[g]         # W(M\g, G)
      W_G_Mg  <- cl_col - col_W[g]         # W(G,  M\g)

      t1 <- if (W_g_G  > 0 && W_G_Mg > 0) W_g_Mg^2 / (W_g_G  * W_G_Mg) else 0
      t2 <- if (W_Mg_G > 0 && W_G_g  > 0) W_Mg_g^2 / (W_Mg_G * W_G_g ) else 0
      gci <- gci + t1 + t2
    }
  }
  gci
}

# ---- Base fitness -------------------------------------------

#' Base fitness for cross-method comparison (Equation 8)
#'
#' The unpenalised measure used to compare AdSA-ALMD against TTWA, GEA and MSA
#' on equal terms: cohesion multiplied by the number of markets.
#'
#' @details
#' \deqn{f_{base}(x) = m \cdot GCI(x)}
#'
#' The multiplier is what makes the measure discriminating. Cohesion alone is
#' maximised by merging everything into one market; multiplying by \code{m}
#' rewards partitions that achieve high cohesion \emph{while staying split},
#' which is the actual delineation problem.
#'
#' This is deliberately \emph{not} the function the search maximises. The
#' search uses \code{\link{compute_objective}()}, which adds self-containment
#' and population penalties and a hard contiguity constraint. Base fitness is
#' reported afterwards so that methods with different internal objectives can
#' still be ranked on one scale (paper Table 3).
#'
#' @param sol Integer vector of market assignments (length N).
#' @param W Numeric \code{N x N} origin-destination matrix.
#' @param row_W Numeric vector, \code{rowSums(W)}.
#' @param col_W Numeric vector, \code{colSums(W)}.
#' @return Numeric scalar \eqn{\ge 0}.
#'
#' @family LMA metrics
#' @seealso \code{\link{compute_objective}()} for the function actually
#'   optimised; \code{\link{build_results_table}()} and
#'   \code{\link{summarise_runs}()}, which report this value.
#'
#' @examples
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#'
#' compute_base_fitness(c(1, 1, 1, 2, 2, 2), W, row_W, col_W)
#'
#' ## One giant market maximises cohesion but not base fitness
#' compute_gci(rep(1, 6), W, row_W, col_W)
#' compute_base_fitness(rep(1, 6), W, row_W, col_W)
#'
#' @export
compute_base_fitness <- function(sol, W, row_W, col_W) {
  n_clusters(sol) * compute_gci(sol, W, row_W, col_W)
}
