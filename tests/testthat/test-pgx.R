## Unit tests for pharmacogenomic (CPIC) findings
## (retrieve_pgx_calls, assign_pgx_phenotypes, assign_pgx_recommendations)
##
## Reference data: small CPIC fixtures (tests/testthat/fixtures/pgx),
## created with data-raw/pgx_test_fixtures.R (GRCh38 records)
##
## Coverage goals
## ─────────────────────────────────────────────────────────────────────────────
## Variant level (retrieve_pgx_calls)
##  1. Altered function allele (DPYD *2A)       → reported, function/activity
##  2. Normal function allele (DPYD c.85T>C)    → not reported (unless asked)
##  3. hom_ref / undefined genotypes            → ignored
##  4. 'chr'-prefixed chromosome names          → matched
##  5. Multi-variant allele, one variant (TPMT) → *3C complete, *3A partial
##  6. Variant completely defining a normal allele (G6PD c.376A>G, 'A')
##                                              → partial A- matches ignored
##  7. No PGx variants                          → zero rows
## Phenotypes (assign_pgx_phenotypes)
##  8. DPYD *2A het / hom_alt                   → IM (1.0) / PM (0.0),
##     CPIC-style diplotype (Reference/c.1905+1G>A (*2A)), variant labels
##  9. DPYD *2A + c.2846A>T (het, het)          → PM (0.5), assumed in trans
## 10. DPYD c.1236G>A only / with c.1129-5923C>G → Not determined / IM (1.5)
## 11. More than two altered function alleles  → Not determined
## 12. TPMT *3C hom_alt                         → PM
## 13. TPMT c.460G>A + c.719A>G                  → *3A (IM), phase note
## 14. NUDT15 c.415C>T (repeat not assessable)  → *3, IM
## 15. G6PD class II: male (hom_alt / het), female (het / hom_alt), unknown
##                                              → Deficient / Variable /
##                                                Deficient / Not determined
## 16. G6PD A- (c.202G>A + c.376A>G), male      → Deficient
## 17. Normal function alleles only             → no phenotypes
## Recommendations (assign_pgx_recommendations)
## 18. DPYD IM (1.0)                            → fluoropyrimidines, -50% dose
## 19. TPMT IM only                             → thiopurines, NUDT15 assumed NM
## 20. TPMT IM + NUDT15 IM                      → thiopurines, 20-50% dose
## 21. G6PD Deficient                           → rasburicase, avoid use
## 22. Phenotype not determined                 → no recommendation, with note
## 23. No phenotypes                            → zero rows
## 24. canonical_pgx_key is order-insensitive

# ---------------------------------------------------------------------------
# Fixtures and helpers
# ---------------------------------------------------------------------------
read_pgx_fixture <- function(tbl) {
  as.data.frame(readr::read_tsv(
    testthat::test_path("fixtures", "pgx", paste0(tbl, ".tsv")),
    na = ".", col_types = readr::cols(.default = "c"),
    show_col_types = FALSE))
}

ref_data <- list(variant = list(
  pgx_allele = read_pgx_fixture("pgx_allele"),
  pgx_phenotype = read_pgx_fixture("pgx_phenotype"),
  pgx_guideline = read_pgx_fixture("pgx_guideline")))

pgx_var_id <- function(symbol, variant_name) {
  pa <- ref_data$variant$pgx_allele
  unique(pa$VAR_ID[pa$SYMBOL == symbol & pa$VARIANT_NAME == variant_name &
                     !is.na(pa$VAR_ID)])[1]
}

## GRCh38 VCF records (CHROM_POS_REF_ALT) of test variants
test_vars <- c(
  dpyd_2a = "1_97450058_C_T",
  dpyd_c1236 = "1_97573863_C_T",
  dpyd_c1129 = "1_97579893_G_C",
  dpyd_c85 = "1_97883329_A_G",
  dpyd_c2846 = pgx_var_id("DPYD", "c.2846A>T"),
  tpmt_c719 = "6_18130687_T_C",
  tpmt_c460 = "6_18138997_C_T",
  nudt15_c415 = "13_48045719_C_T",
  g6pd_andalus = "X_154532389_C_T",
  g6pd_c376 = "X_154535277_T_C",
  g6pd_c202 = pgx_var_id("G6PD", "c.202G>A")
)

## sample calls, e.g. make_calls(dpyd_2a = "het", tpmt_c719 = "hom_alt")
make_calls <- function(..., chr_prefix = FALSE) {
  gts <- c(...)
  do.call(rbind, lapply(names(gts), function(v) {
    k <- strsplit(test_vars[[v]], "_")[[1]]
    data.frame(
      CHROM = if (chr_prefix) paste0("chr", k[1]) else k[1],
      POS = as.integer(k[2]), REF = k[3], ALT = k[4],
      GENOTYPE = gts[[v]],
      SYMBOL = toupper(sub("_.*", "", v)),
      stringsAsFactors = FALSE)
  }))
}

phenotype_of <- function(calls, gene, sex = "UNKNOWN") {
  ph <- cpsr::assign_pgx_phenotypes(calls, ref_data = ref_data, sex = sex)
  ph[ph$SYMBOL == gene, ]
}

test_that("fixtures contain all test variants", {
  expect_false(any(is.na(test_vars)))
  expect_true(all(test_vars %in% ref_data$variant$pgx_allele$VAR_ID))
})

# ---------------------------------------------------------------------------
# Variant level - retrieve_pgx_calls
# ---------------------------------------------------------------------------
test_that("1. altered function allele is reported with CPIC annotation", {
  out <- cpsr::retrieve_pgx_calls(make_calls(dpyd_2a = "het"), ref_data)
  expect_equal(NROW(out), 1)
  expect_equal(out$PGX_ALLELE, "c.1905+1G>A (*2A)")
  expect_equal(out$PGX_ALLELE_FUNCTION, "No function")
  expect_equal(out$PGX_ACTIVITY_VALUE, "0.0")
  expect_true(out$PGX_ALLELE_COMPLETE)
})

test_that("2. normal function alleles are only reported on request", {
  calls <- make_calls(dpyd_c85 = "hom_alt")
  expect_equal(NROW(cpsr::retrieve_pgx_calls(calls, ref_data)), 0)
  out <- cpsr::retrieve_pgx_calls(
    calls, ref_data, include_normal_function = TRUE)
  expect_equal(out$PGX_ALLELE_FUNCTION, "Normal function")
})

test_that("3. hom_ref and undefined genotypes are ignored", {
  calls <- make_calls(dpyd_2a = "hom_ref", tpmt_c719 = "undefined")
  expect_equal(NROW(cpsr::retrieve_pgx_calls(calls, ref_data)), 0)
  expect_equal(NROW(cpsr::assign_pgx_phenotypes(calls, ref_data)), 0)
})

test_that("4. chr-prefixed chromosome names are matched", {
  out <- cpsr::retrieve_pgx_calls(
    make_calls(dpyd_2a = "het", chr_prefix = TRUE), ref_data)
  expect_equal(NROW(out), 1)
})

test_that("5. single variant of a multi-variant allele - complete and partial", {
  out <- cpsr::retrieve_pgx_calls(make_calls(tpmt_c719 = "het"), ref_data)
  expect_match(out$PGX_ALLELE, "^\\*3C")
  expect_match(out$PGX_ALLELE, "\\*3A \\(partial - 1/2 defining variants\\)")
  expect_true(out$PGX_ALLELE_COMPLETE)
})

test_that("6. variant completely defining a normal allele is not reported", {
  calls <- make_calls(g6pd_c376 = "hom_alt")
  expect_equal(NROW(cpsr::retrieve_pgx_calls(calls, ref_data)), 0)
  expect_equal(NROW(
    cpsr::assign_pgx_phenotypes(calls, ref_data, sex = "MALE")), 0)
})

test_that("7. no PGx variants gives zero rows", {
  calls <- data.frame(
    CHROM = "17", POS = 43045712L, REF = "T", ALT = "C",
    GENOTYPE = "het", SYMBOL = "BRCA1")
  out <- cpsr::retrieve_pgx_calls(calls, ref_data)
  expect_equal(NROW(out), 0)
})

# ---------------------------------------------------------------------------
# Phenotypes - assign_pgx_phenotypes
# ---------------------------------------------------------------------------
test_that("8. DPYD *2A heterozygous / homozygous", {
  het <- phenotype_of(make_calls(dpyd_2a = "het"), "DPYD")
  expect_equal(het$PGX_PHENOTYPE, "Intermediate Metabolizer")
  expect_equal(het$PGX_ACTIVITY_SCORE, "1.0")
  expect_equal(het$PGX_DIPLOTYPE, "Reference/c.1905+1G>A (*2A)")
  expect_equal(het$PGX_VARIANTS, "c.1905+1G>A")
  expect_equal(het$PGX_GENOTYPES, "heterozygous")
  expect_true(is.na(het$PGX_PHENOTYPE_NOTE))

  hom <- phenotype_of(make_calls(dpyd_2a = "hom_alt"), "DPYD")
  expect_equal(hom$PGX_PHENOTYPE, "Poor Metabolizer")
  expect_equal(hom$PGX_ACTIVITY_SCORE, "0.0")
})

test_that("8b. variant labels use HGVSc from sample annotation", {
  calls <- make_calls(tpmt_c719 = "het")
  calls$HGVSc <- "ENST00000309983.5:c.719A>G"
  ph <- phenotype_of(calls, "TPMT")
  expect_equal(ph$PGX_VARIANTS, "c.719A>G")
  expect_equal(ph$PGX_GENOTYPES, "heterozygous")
  expect_equal(ph$PGX_DIPLOTYPE, "*1/*3C")
})

test_that("9. two different DPYD alleles are assumed in trans", {
  ph <- phenotype_of(make_calls(dpyd_2a = "het", dpyd_c2846 = "het"), "DPYD")
  expect_equal(ph$PGX_PHENOTYPE, "Poor Metabolizer")
  expect_equal(ph$PGX_ACTIVITY_SCORE, "0.5")
  expect_match(ph$PGX_PHENOTYPE_NOTE, "different gene copies \\(in trans\\)")
  expect_equal(ph$PGX_DIPLOTYPE, "c.1905+1G>A (*2A)/c.2846A>T")
  expect_equal(ph$PGX_VARIANTS, "c.1905+1G>A; c.2846A>T")
  expect_equal(ph$PGX_GENOTYPES, "heterozygous; heterozygous")
})

test_that("10. DPYD HapB3 - partial vs. complete", {
  partial <- phenotype_of(make_calls(dpyd_c1236 = "het"), "DPYD")
  expect_equal(partial$PGX_PHENOTYPE, "Not determined")
  expect_match(partial$PGX_PHENOTYPE_NOTE, "HapB3")

  complete <- phenotype_of(
    make_calls(dpyd_c1236 = "het", dpyd_c1129 = "het"), "DPYD")
  expect_equal(complete$PGX_PHENOTYPE, "Intermediate Metabolizer")
  expect_equal(complete$PGX_ACTIVITY_SCORE, "1.5")
  expect_match(complete$PGX_DIPLOTYPE, "HapB3")
})

test_that("11. more than two altered function alleles - not determined", {
  ph <- phenotype_of(
    make_calls(dpyd_2a = "hom_alt", dpyd_c2846 = "het"), "DPYD")
  expect_equal(ph$PGX_PHENOTYPE, "Not determined")
  expect_match(ph$PGX_PHENOTYPE_NOTE, "more altered function alleles")
})

test_that("12. TPMT *3C homozygous", {
  ph <- phenotype_of(make_calls(tpmt_c719 = "hom_alt"), "TPMT")
  expect_equal(ph$PGX_PHENOTYPE, "Poor Metabolizer")
  expect_equal(ph$PGX_DIPLOTYPE, "*3C/*3C")
  expect_equal(ph$PGX_GENOTYPES, "homozygous")
})

test_that("13. TPMT c.460G>A + c.719A>G collapse into *3A, with phase note", {
  ph <- phenotype_of(make_calls(tpmt_c460 = "het", tpmt_c719 = "het"), "TPMT")
  expect_equal(ph$PGX_PHENOTYPE, "Intermediate Metabolizer")
  expect_equal(ph$PGX_DIPLOTYPE, "*1/*3A")
  expect_match(ph$PGX_PHENOTYPE_NOTE, "\\*3B and \\*3C on different gene copies")
})

test_that("14. NUDT15 *3 with non-assessable repeat variant", {
  ph <- phenotype_of(make_calls(nudt15_c415 = "het"), "NUDT15")
  expect_equal(ph$PGX_PHENOTYPE, "Intermediate Metabolizer")
  expect_match(ph$PGX_DIPLOTYPE, "defining variants assessable")
})

test_that("15. G6PD class II allele depends on sex", {
  male_hap <- phenotype_of(make_calls(g6pd_andalus = "hom_alt"), "G6PD", "MALE")
  expect_equal(male_hap$PGX_PHENOTYPE, "Deficient")
  expect_match(male_hap$PGX_DIPLOTYPE, "hemizygous")

  male_het <- phenotype_of(make_calls(g6pd_andalus = "het"), "G6PD", "MALE")
  expect_equal(male_het$PGX_PHENOTYPE, "Deficient")
  expect_match(male_het$PGX_PHENOTYPE_NOTE, "male sample")

  expect_equal(male_hap$PGX_GENOTYPES, "hemizygous")

  female_het <- phenotype_of(make_calls(g6pd_andalus = "het"), "G6PD", "FEMALE")
  expect_equal(female_het$PGX_PHENOTYPE, "Variable")
  expect_equal(female_het$PGX_DIPLOTYPE, "B (reference)/Andalus")

  female_hom <- phenotype_of(
    make_calls(g6pd_andalus = "hom_alt"), "G6PD", "FEMALE")
  expect_equal(female_hom$PGX_PHENOTYPE, "Deficient")

  unknown <- phenotype_of(
    make_calls(g6pd_andalus = "hom_alt"), "G6PD", "UNKNOWN")
  expect_equal(unknown$PGX_PHENOTYPE, "Not determined")
  expect_match(unknown$PGX_PHENOTYPE_NOTE, "sex is unknown")
})

test_that("16. G6PD A- (c.202G>A + c.376A>G) in a male", {
  ph <- phenotype_of(
    make_calls(g6pd_c202 = "hom_alt", g6pd_c376 = "hom_alt"), "G6PD", "MALE")
  expect_equal(ph$PGX_PHENOTYPE, "Deficient")
  expect_match(ph$PGX_DIPLOTYPE, "^A-")
})

test_that("17. normal function alleles only - no phenotypes", {
  calls <- make_calls(dpyd_c85 = "hom_alt", g6pd_c376 = "hom_alt")
  expect_equal(NROW(
    cpsr::assign_pgx_phenotypes(calls, ref_data, sex = "MALE")), 0)
})

# ---------------------------------------------------------------------------
# Recommendations - assign_pgx_recommendations
# ---------------------------------------------------------------------------
recommendations_for <- function(calls, sex = "UNKNOWN") {
  ph <- cpsr::assign_pgx_phenotypes(calls, ref_data = ref_data, sex = sex)
  cpsr::assign_pgx_recommendations(ph, ref_data = ref_data)
}

test_that("18. DPYD Intermediate Metabolizer - reduce fluoropyrimidine dose", {
  rc <- recommendations_for(make_calls(dpyd_2a = "het"))
  expect_setequal(rc$DRUG_NAME, c("capecitabine", "fluorouracil"))
  expect_true(all(rc$LOOKUP_KEY == "DPYD:1.0"))
  expect_true(all(grepl("^Reduce starting dose by 50%", rc$RECOMMENDATION)))
  expect_true(all(rc$CLASSIFICATION == "Strong"))
})

test_that("19. TPMT only - NUDT15 assumed normal metabolizer", {
  rc <- recommendations_for(make_calls(tpmt_c719 = "het"))
  expect_setequal(rc$DRUG_NAME, c("mercaptopurine", "thioguanine"))
  expect_true(all(grepl("TPMT:Intermediate Metabolizer", rc$LOOKUP_KEY)))
  expect_true(all(grepl("NUDT15:Normal Metabolizer", rc$LOOKUP_KEY)))
  expect_true(all(grepl("NUDT15 assumed Normal Metabolizer", rc$PGX_NOTE)))
  expect_true(all(grepl("30-80%", rc$RECOMMENDATION)))
})

test_that("20. TPMT and NUDT15 intermediate metabolizers - 20-50% dose", {
  rc <- recommendations_for(make_calls(tpmt_c719 = "het", nudt15_c415 = "het"))
  mp <- rc[rc$DRUG_NAME == "mercaptopurine", ]
  expect_equal(NROW(mp), 1)
  expect_match(mp$RECOMMENDATION, "20-50%")
  expect_true(is.na(mp$PGX_NOTE))
})

test_that("21. G6PD deficient - avoid rasburicase", {
  rc <- recommendations_for(make_calls(g6pd_andalus = "hom_alt"), sex = "MALE")
  expect_equal(rc$DRUG_NAME, "rasburicase")
  expect_equal(rc$LOOKUP_KEY, "G6PD:Deficient")
  expect_equal(rc$RECOMMENDATION, "Avoid use")
})

test_that("22. phenotype not determined - no recommendation, with note", {
  rc <- recommendations_for(make_calls(dpyd_c1236 = "het"))
  expect_setequal(rc$DRUG_NAME, c("capecitabine", "fluorouracil"))
  expect_true(all(is.na(rc$RECOMMENDATION)))
  expect_true(all(grepl("phenotype not determined for DPYD", rc$PGX_NOTE)))
})

test_that("23. no phenotypes - no recommendations", {
  expect_equal(NROW(
    cpsr::assign_pgx_recommendations(data.frame(), ref_data = ref_data)), 0)
})

test_that("24. canonical_pgx_key is order-insensitive", {
  expect_equal(
    cpsr:::canonical_pgx_key(c("0.5:1", "1.0:1")),
    cpsr:::canonical_pgx_key(c("1.0:1", "0.5:1")))
  expect_equal(
    cpsr:::canonical_pgx_key(c("TPMT:Normal Metabolizer",
                               "NUDT15:Poor Metabolizer")),
    "NUDT15:Poor Metabolizer;TPMT:Normal Metabolizer")
})
