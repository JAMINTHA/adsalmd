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
