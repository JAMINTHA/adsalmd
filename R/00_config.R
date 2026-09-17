# ============================================================
# FILE: R/00_config.R
# AdSA-ALMD Configuration — All tunable parameters in one place.
# Modify ONLY this file to change algorithm behaviour.
# ============================================================

#' Default AdSA-ALMD configuration
#'
#' The tunable parameters of the method in one place, set to the values used
#' for the Queensland SA2 study (N = 520 BGUs, 4 percent population sample).
#' Every algorithm in the package takes these as ordinary arguments — this list
#' is a convenient bundle of defaults, not a hidden global that functions read.
#'
#' @details
#' Copy and override individual entries rather than editing in place, since the
#' object is shared across everything that reads it:
#'
#' \preformatted{
#' params <- AdSA_params
#' params$SC_min <- 0.60
#' }
#'
#' Entries:
#' \describe{
#'   \item{\code{N}}{Integer. Number of BGUs after merging uninhabited SA2s.
#'     Informational — the algorithms take \code{N} from
#'     \code{length(row_W)}.}
#'   \item{\code{SC_min}}{Numeric in (0, 1]. Minimum self-containment a market
#'     must reach to count as valid. Drives \code{\link{validity_score}()},
#'     \code{\link{penalty_sc_single}()} and the SARA merge loop.}
#'   \item{\code{SC_tar}}{Numeric. Target (not minimum) self-containment, used
#'     when reporting how far above the floor a solution sits.}
#'   \item{\code{Pop_min}, \code{Pop_tar}, \code{Pop_max}}{Numeric population
#'     floor, target and ceiling per market. \code{Pop_min} and \code{Pop_max}
#'     are hard constraints in \code{\link{compute_objective}()} — a partition
#'     violating either scores 0.}
#'   \item{\code{penalty_exp}}{Integer \code{r}, the steepness exponent of the
#'     penalty functions (the paper uses 3 or 4). Larger \code{r} makes small
#'     shortfalls cheaper and large ones much more expensive.}
#'   \item{\code{P0}, \code{T0_samples}}{Target acceptance probability at the
#'     initial temperature, and the number of random perturbations used to
#'     estimate it. See \code{\link{estimate_t0}()}.}
#'   \item{\code{L}, \code{l}, \code{eps}}{Outer trials, inner iterations per
#'     trial, and the convergence threshold of \code{\link{adsa_almd}()}.}
#'   \item{\code{alpha}, \code{beta}}{Adaptive cooling coefficients: how fast
#'     the temperature rises on a rejection streak and falls on an acceptance
#'     streak.}
#'   \item{\code{n_runs}, \code{n_cores}}{Independent runs and worker processes
#'     for \code{\link{run_adsa_pipeline_parallel}()}.}
#' }
#'
#' @format A named list of 16 tunable parameters.
#'
#' @seealso \code{\link{suggest_thresholds}()} to derive \code{SC_min} /
#'   \code{Pop_min} / \code{Pop_max} / \code{Pop_tar} from a partition you
#'   already trust; \code{\link{run_adsa_pipeline}()} for the entry point that
#'   consumes them.
#'
#' @examples
#' str(AdSA_params)
#'
#' ## Override without mutating the shared default
#' params <- AdSA_params
#' params$SC_min <- 0.60
#' params$L <- 50L
#' c(default = AdSA_params$SC_min, mine = params$SC_min)
#'
#' @export
AdSA_params <- list(

  # ---- Data ----
  N            = 520L,    # BGUs after merging uninhabited SA2s

  # ---- Self-containment thresholds (Table 2 grid search optimum) ----
  SC_min       = 0.55,   # Minimum self-containment threshold
  SC_tar       = 0.75,   # Target  self-containment threshold

  # ---- Population thresholds — 4 % sample of QLD (Table 2) ----
  Pop_min      = 500L,   # Minimum population
  Pop_tar      = 2000L,  # Target  population
  Pop_max      = 4000L,
  
  # ---- Penalty function (paper: r = 3 or 4) ----
  penalty_exp  = 3L,     # Exponent r controlling steepness

  # ---- Initial temperature ----
  P0           = 0.9,    # Target acceptance probability at T0
  T0_samples   = 500L,   # Perturbations used to estimate E[|DeltaF|]

  # ---- AdSA-ALMD inner / outer loops (Algorithm 1) ----
  L            = 1000L,  # Number of outer trials
  l            = 25L,    # Inner iterations per trial
  eps          = 0,      # Convergence threshold epsilon

  # ---- Adaptive Cooling Schedule (Equation 15) ----
  alpha        = 0.2,    # Temperature-increase coefficient (consecutive rejects)
  beta         = 0.2,    # Temperature-decrease coefficient (consecutive accepts)

  # ---- Parallel execution ----
  n_runs       = 10L,    # Independent runs (different random seeds)
  n_cores      = 8L      # CPU cores available in ABS DataLab
)
