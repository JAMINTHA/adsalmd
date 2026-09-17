# =============================================================================
# dev/sample_data/generate_sample_data.R
#
# Regenerates the two CSVs in this folder:
#   census_extract_sample.csv      -- synthetic person-level Census extract
#   DZN_SA2_2016_AUST_sample.csv   -- matching DZN -> SA2 lookup
#
# Same 10-SA2 Brisbane scheme as dev/test_census_convert.R, written out as
# plain CSV (base R only) so USER_GUIDE.md section 1 can read.csv() them
# directly -- no invented data needed to try the walkthrough before you have
# a real DataLab extract.
# =============================================================================

set.seed(42)
n <- 1000

sa2_codes <- c(
  "301011001", "301011002", "301011003", "301011004", "301011005",
  "301011006", "301011007", "301011008", "301011009", "301011010"
)

dzn_sa2_lookup <- data.frame(
  DZN_CODE_2016 = c(
    "30101100101", "30101100102",
    "30101100201", "30101100202",
    "30101100301", "30101100302",
    "30101100401",
    "30101100501", "30101100502",
    "30101100601",
    "30101100701", "30101100702",
    "30101100801",
    "30101100901", "30101100902",
    "30101101001"
  ),
  SA2_MAINCODE_2016 = c(
    "301011001", "301011001",
    "301011002", "301011002",
    "301011003", "301011003",
    "301011004",
    "301011005", "301011005",
    "301011006",
    "301011007", "301011007",
    "301011008",
    "301011009", "301011009",
    "301011010"
  ),
  stringsAsFactors = FALSE
)

census_raw <- data.frame(
  ABSPID = sprintf("P%07d", 1:n),
  AGEP   = sample(15:85, n, replace = TRUE),
  LFSP   = sample(c("1", "2", "3", "4", "5", "6", "@"), n, replace = TRUE,
                  prob = c(0.35, 0.20, 0.05, 0.05, 0.03, 0.25, 0.07)),
  PURP   = sample(sa2_codes, n, replace = TRUE,
                  prob = c(0.15, 0.12, 0.10, 0.08, 0.10,
                           0.12, 0.11, 0.08, 0.07, 0.07)),
  POWP   = sample(dzn_sa2_lookup$DZN_CODE_2016, n, replace = TRUE),
  IFPOWP = sample(c(0L, 1L), n, replace = TRUE, prob = c(0.889, 0.111)),
  stringsAsFactors = FALSE
)

write.csv(census_raw, file.path("dev", "sample_data", "census_extract_sample.csv"),
          row.names = FALSE)
write.csv(dzn_sa2_lookup, file.path("dev", "sample_data", "DZN_SA2_2016_AUST_sample.csv"),
          row.names = FALSE)

cat(sprintf("Wrote %d person records and %d DZN->SA2 rows to dev/sample_data/\n",
            nrow(census_raw), nrow(dzn_sa2_lookup)))



