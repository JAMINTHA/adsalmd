# =============================================================================
# dev/abs_adjacency_from_shapefile.R
#
# Build the geographic adjacency matrix `adj` expected by run_adsa_pipeline()
# from an ABS SA2 digital boundary shapefile. NOT part of the public adsalmd
# package -- source this file manually:
#
#   source("dev/abs_adjacency_from_shapefile.R")
#
# Requires: sf, spdep. Download boundaries from the ABS (e.g.
# SA2_2016_AUST_GDA2020.shp / later ASGS releases). NOT declared as package
# dependencies -- see dev/README.md.
# =============================================================================

#' Build a binary SA2 adjacency matrix from an ABS boundary shapefile
#'
#' Reads an SA2 shapefile, subsets and reorders it to `sa2_codes`, and derives
#' a queen-contiguity adjacency matrix with \code{spdep::poly2nb()} /
#' \code{spdep::nb2mat()}. The returned matrix's row/column order matches
#' `sa2_codes` exactly, so pass the *same* vector you gave
#' \code{\link{build_od_matrix}} to keep `W` and `adj` aligned.
#'
#' @param shapefile_path Path to the SA2 boundary shapefile (`.shp`), e.g. an
#'   ABS ASGS SA2 digital boundary file.
#' @param sa2_codes Character vector, the canonical, ordered set of SA2 codes
#'   for this run (same vector passed to \code{\link{build_od_matrix}}).
#' @param sa2_id_col Column in the shapefile holding the SA2 code (default
#'   `"SA2_MAIN16"`; use `"SA2_MAIN21"`/`"SA2_CODE21"` etc. for newer ASGS
#'   editions -- check `names(sf::st_read(shapefile_path))`).
#' @param queen Logical (default TRUE); TRUE = queen contiguity (shared edge
#'   or vertex), FALSE = rook contiguity (shared edge only).
#' @return Integer `length(sa2_codes) x length(sa2_codes)` binary matrix,
#'   symmetric, zero diagonal, dimnames set to `sa2_codes` on both margins.
#'   Errors if any `sa2_codes` entry is missing from the shapefile.
build_adjacency_from_shapefile <- function(shapefile_path, sa2_codes,
                                            sa2_id_col = "SA2_MAIN16",
                                            queen = TRUE) {
  stopifnot(requireNamespace("sf", quietly = TRUE))
  stopifnot(requireNamespace("spdep", quietly = TRUE))
  stopifnot(!anyDuplicated(sa2_codes))

  sa2_shp <- sf::st_read(shapefile_path, quiet = TRUE)
  if (!sa2_id_col %in% names(sa2_shp)) {
    stop(sprintf(
      "`%s` not found in shapefile. Available columns: %s",
      sa2_id_col, paste(names(sa2_shp), collapse = ", ")
    ))
  }

  sa2_shp <- sa2_shp[sa2_shp[[sa2_id_col]] %in% sa2_codes, ]
  sa2_shp <- sf::st_make_valid(sa2_shp)

  missing <- setdiff(sa2_codes, sa2_shp[[sa2_id_col]])
  if (length(missing) > 0) {
    stop(sprintf(
      "%d sa2_codes not found in shapefile: %s",
      length(missing), paste(missing, collapse = ", ")
    ))
  }

  # Reorder rows to match sa2_codes exactly before deriving neighbours, so
  # nb2mat()'s output rows/cols line up with sa2_codes by position.
  sa2_shp <- sa2_shp[match(sa2_codes, sa2_shp[[sa2_id_col]]), ]

  nb <- spdep::poly2nb(sa2_shp, queen = queen)
  A  <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)

  dimnames(A) <- list(sa2_codes, sa2_codes)
  storage.mode(A) <- "integer"
  A
}

#' Reattach zero-flow BGUs to the nearest included neighbour's LMA
#'
#' Some SA2s have no residents and no jobs recorded at all (a national park,
#' an industrial or port zone, an uninhabited island) and so get excluded
#' before \code{run_adsa_pipeline()} runs (see \code{USER_GUIDE.md}, step
#' 1.5) -- they have no commuting flow to place them by. This reattaches
#' each one to the labour market area (LMA) of the geographically nearest
#' SA2 that *was* included in the run: first by shared border (queen/rook
#' contiguity, matching \code{\link{build_adjacency_from_shapefile}}), and
#' if none of its borders made it into the run either, by nearest centroid
#' distance among the SA2s that did.
#'
#' @param empty_codes Character vector of SA2 codes to reattach (e.g.
#'   `empty_bgu` from \code{USER_GUIDE.md} step 1.5). Empty input returns an
#'   empty result -- always safe to call, even if nothing was set aside.
#' @param sol Integer vector, the pipeline's partition (e.g.
#'   \code{res$best_sol}), aligned to `sa2_codes`.
#' @param sa2_codes Character vector, the trimmed SA2 codes actually used
#'   for the pipeline run (i.e. `rownames(W)`/`rownames(adj)`), same order
#'   and length as `sol`.
#' @param shapefile_path Path to the SA2 boundary shapefile -- the *full*
#'   one, not subset to `sa2_codes`, since `empty_codes`'s neighbours need
#'   to be found too.
#' @param sa2_id_col Column in the shapefile holding the SA2 code (default
#'   `"SA2_MAIN16"`; match whatever you used for
#'   \code{\link{build_adjacency_from_shapefile}}).
#' @param queen Logical (default TRUE); TRUE = queen contiguity (shared edge
#'   or vertex), FALSE = rook contiguity (shared edge only) -- match whatever
#'   you used to build `adj`.
#' @return Named integer vector, length `length(empty_codes)`: names are
#'   `empty_codes`, values are each one's assigned LMA label (a value from
#'   `sol`). Combine with `sol`/`sa2_codes` (e.g. `c(sol, result)`,
#'   `c(sa2_codes, empty_codes)`) before building the final correspondence
#'   table, so every SA2 in the study area is covered.
reattach_isolated_bgus <- function(empty_codes, sol, sa2_codes,
                                    shapefile_path,
                                    sa2_id_col = "SA2_MAIN16",
                                    queen = TRUE) {
  if (length(empty_codes) == 0) {
    return(stats::setNames(integer(0), character(0)))
  }

  stopifnot(requireNamespace("sf", quietly = TRUE))
  stopifnot(requireNamespace("spdep", quietly = TRUE))
  stopifnot(length(sol) == length(sa2_codes))

  sa2_shp <- sf::st_read(shapefile_path, quiet = TRUE)
  if (!sa2_id_col %in% names(sa2_shp)) {
    stop(sprintf(
      "`%s` not found in shapefile. Available columns: %s",
      sa2_id_col, paste(names(sa2_shp), collapse = ", ")
    ))
  }
  sa2_shp <- sf::st_make_valid(sa2_shp)
  codes_all <- sa2_shp[[sa2_id_col]]

  missing <- setdiff(c(empty_codes, sa2_codes), codes_all)
  if (length(missing) > 0) {
    stop(sprintf(
      "%d code(s) not found in shapefile: %s",
      length(missing), paste(missing, collapse = ", ")
    ))
  }

  nb <- spdep::poly2nb(sa2_shp, queen = queen)

  nearest_included <- function(empty_code) {
    i <- match(empty_code, codes_all)
    bordering  <- codes_all[nb[[i]]]
    candidates <- intersect(bordering, sa2_codes)

    if (length(candidates) == 0) {
      # no bordering SA2 was included in the run either (rare) -- fall back
      # to the nearest centroid among included SA2s
      d <- sf::st_distance(
        sf::st_centroid(sa2_shp[i, ]),
        sf::st_centroid(sa2_shp[match(sa2_codes, codes_all), ])
      )
      candidates <- sa2_codes[which.min(d)]
    }
    candidates[1]
  }

  nearest <- vapply(empty_codes, nearest_included, character(1))
  stats::setNames(sol[match(nearest, sa2_codes)], empty_codes)
}
