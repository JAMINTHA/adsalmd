# ============================================================
# FILE: R/adsalmd-package.R
# Package-level documentation (?adsalmd) and shared roxygen
# building blocks.
# ============================================================

#' adsalmd: Adaptive Simulated Annealing Delineation of Local Labour Markets
#'
#' Delineates local labour market areas (LMAs, also called "markets") from an
#' origin-destination commuting matrix and a spatial adjacency matrix. A SARA
#' (Spatially Adaptive Random Aggregation) pass builds a feasible initial
#' partition, AdSA-ALMD simulated annealing refines it, and a deterministic
#' repair pass fixes any remaining cohesion-misallocated BGUs.
#'
#' @section The three-stage workflow:
#' \describe{
#'   \item{1. Initial partition}{\code{\link{sara_iterated}()} repeatedly calls
#'     \code{\link{sara}()}, merging the least-valid cluster into one of its
#'     best-ranked neighbours until every market satisfies the validity
#'     criterion (\code{\link{validity_score}}).}
#'   \item{2. Annealing}{\code{\link{adsa_almd}()} perturbs the partition with
#'     eleven group-based operators (\code{\link{operator_1}} …
#'     \code{\link{operator_11}}), accepting moves under an adaptive cooling
#'     schedule and maximising the objective \code{\link{compute_objective}()}.}
#'   \item{3. Refinement}{\code{\link{refine_misallocated_bgus}()} reassigns any
#'     BGU that would gain cohesion elsewhere without a substantial
#'     self-containment loss.}
#' }
#' \code{\link{run_adsa_pipeline}()} chains all three in one call, and
#' \code{\link{run_adsa_pipeline_parallel}()} runs it once per seed across a
#' PSOCK cluster and keeps the best result.
#'
#' @section Notation used throughout the help pages:
#' \describe{
#'   \item{\code{N}}{number of BGUs (basic geographic units) in the study area.}
#'   \item{\code{W}}{\code{N x N} origin-destination matrix; \code{W[i, j]} is
#'     the number of commuters living in BGU \code{i} and working in BGU
#'     \code{j}.}
#'   \item{\code{adj}}{\code{N x N} symmetric adjacency matrix; non-zero means
#'     the two BGUs share a border.}
#'   \item{\code{row_W}, \code{col_W}}{\code{rowSums(W)} and \code{colSums(W)},
#'     computed once by the caller and threaded through every function that
#'     needs them. Passing stale values silently produces wrong metrics.}
#'   \item{\code{sol}}{integer vector of length \code{N}; \code{sol[i]} is the
#'     market that BGU \code{i} belongs to. Labels are arbitrary — compare two
#'     partitions with \code{\link{partitions_match}()}, not \code{identical()}.}
#'   \item{\code{Ms}}{the set of BGUs in market \code{s}; \code{G} the full BGU
#'     set \code{1:N}; \code{M\\g} market \code{M} excluding BGU \code{g}.}
#'   \item{\code{m}}{number of markets in a partition,
#'     \code{\link{n_clusters}(sol)}.}
#' }
#'
#' @section Choosing thresholds:
#' \code{SC_min}, \code{Pop_min}, \code{Pop_max} and \code{Pop_tar} drive every
#' feasibility test in the package. \code{\link{AdSA_params}} holds the defaults
#' used for the Queensland SA2 study; \code{\link{suggest_thresholds}()} reads
#' plausible values off a partition you already trust.
#'
#' @section Global variables:
#' \code{\link{operator_1}()}, \code{\link{operator_10}()},
#' \code{\link{operator_11}()} and \code{\link{apply_random_operator}()} read a
#' \code{params} list (\code{$SC_min}, \code{$Pop_min}, \code{$Pop_max}) from
#' the global environment rather than taking it as an argument. The pipeline
#' entry points set this up for you; call those operators directly only with a
#' \code{params} object in scope.
#'
#' @references
#' Martínez-Bernabeu, L., Flórez-Revuelta, F. and Casado-Díaz, J. M. (2012)
#' Grouping genetic operators for the delineation of functional areas based on
#' spatial interaction. \emph{Expert Systems with Applications}.
#' (Source of the group-based perturbation operators adapted here.)
#'
#' @keywords internal
#' @importFrom stats median runif sd
#' @importFrom utils write.csv
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ
#'   parLapply detectCores
"_PACKAGE"

#' Toy example data used in the help pages
#'
#' Several examples in this package build the same six-BGU toy system:
#'
#' \preformatted{
#' adj <- matrix(0L, 6, 6); adj[cbind(1:5, 2:6)] <- 1L; adj <- adj + t(adj)
#' W <- matrix(1, 6, 6); W[1:3, 1:3] <- 30; W[4:6, 4:6] <- 30; diag(W) <- 60
#' row_W <- rowSums(W); col_W <- colSums(W)
#' sol <- c(1, 1, 1, 2, 2, 2)
#' }
#'
#' Six BGUs in a line (1-2-3-4-5-6). Commuting is heavy within \code{1:3} and
#' within \code{4:6} and almost absent between them, so \code{sol} is the
#' obvious two-market answer and every metric in the package has an easy-to-
#' check value on it.
#'
#' Larger generated examples (N = 14, 100, 520) with known true partitions live
#' in the \code{examples/} directory of the source tree.
#'
#' @name adsalmd-toy-data
#' @keywords internal
NULL
