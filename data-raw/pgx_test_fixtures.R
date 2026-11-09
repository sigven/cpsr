## Creates small CPIC reference data fixtures for the pharmacogenomics
## unit tests (tests/testthat/test-pgx.R), extracted from a PCGR/CPSR
## data bundle (variant/tsv/pgx, GRCh38).
##
## Fixtures contain:
##  - pgx_allele: all allele definitions that include any of the test
##    variants below (so that partial/complete allele logic is realistic)
##  - pgx_phenotype: all genotype-to-phenotype mappings (DPYD, TPMT,
##    NUDT15, G6PD)
##  - pgx_guideline: all CPIC recommendations for genes in scope of CPSR
##
## Re-run when the bundle (CPIC release) is updated, and check that
## tests/testthat/test-pgx.R still passes.

bundle_pgx_dir <- file.path(
  "/Users/sigven/project_data/data/data__pcgrdb/bundle_output",
  "20261007/data/grch38/variant/tsv/pgx")
fixture_dir <- "tests/testthat/fixtures/pgx"

read_pgx <- function(fname) {
  as.data.frame(readr::read_tsv(
    file.path(bundle_pgx_dir, fname), na = ".",
    col_types = readr::cols(.default = "c"), show_col_types = FALSE))
}
pgx_allele <- read_pgx("pgx_allele.tsv.gz")

## VCF records (GRCh38) used in the test scenarios
test_var_ids <- c(
  "1_97450058_C_T",   # DPYD c.1905+1G>A (*2A) - no function
  "1_97573863_C_T",   # DPYD c.1236G>A - HapB3 (with c.1129-5923C>G)
  "1_97579893_G_C",   # DPYD c.1129-5923C>G - HapB3 / decreased function
  "1_97883329_A_G",   # DPYD c.85T>C (*9A) - normal function
  "6_18130687_T_C",   # TPMT c.719A>G - *3C / *3A
  "6_18138997_C_T",   # TPMT c.460G>A - *3B / *3A
  "13_48045719_C_T",  # NUDT15 c.415C>T - *3 (+ non-VCF repeat)
  "X_154532389_C_T",  # G6PD c.1361G>A (Andalus) - II/Deficient
  "X_154535277_T_C",  # G6PD c.376A>G (A) - IV/Normal, part of A-
  unique(pgx_allele$VAR_ID[
    pgx_allele$SYMBOL == "DPYD" & pgx_allele$VARIANT_NAME == "c.2846A>T"]),
  unique(pgx_allele$VAR_ID[
    pgx_allele$SYMBOL == "G6PD" & pgx_allele$VARIANT_NAME == "c.202G>A"])
)
stopifnot(all(test_var_ids %in% pgx_allele$VAR_ID))

definition_ids <- unique(
  pgx_allele$ALLELE_DEFINITION_ID[pgx_allele$VAR_ID %in% test_var_ids])

fixtures <- list(
  pgx_allele = pgx_allele |>
    dplyr::filter(.data$ALLELE_DEFINITION_ID %in% definition_ids),
  pgx_phenotype = read_pgx("pgx_phenotype.tsv.gz"),
  pgx_guideline = read_pgx("pgx_guideline.tsv.gz") |>
    dplyr::filter(as.logical(.data$CPSR_CALLABLE))
)

dir.create(fixture_dir, recursive = TRUE, showWarnings = FALSE)
for (tbl in names(fixtures)) {
  readr::write_tsv(
    fixtures[[tbl]],
    file = file.path(fixture_dir, paste0(tbl, ".tsv")),
    na = ".", quote = "needed")
  cat(tbl, ":", NROW(fixtures[[tbl]]), "rows\n")
}
