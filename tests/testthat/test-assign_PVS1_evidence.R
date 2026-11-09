## Unit tests for assign_PVS1_evidence
##
## Coverage goals
## ─────────────────────────────────────────────────────────────────────────────
## 1. NULL / stop-gain / frameshift            → PVS1  (no NMD escape)
## 2. NULL_VARIANT + NMD_escaping, domain      → PVS1_STR
## 3. NULL_VARIANT + NMD_escaping, no domain   → PVS1_MOD
## 4. Splice canonical ±1/2, INTRON_POSITION correct → PVS1  (non-last intron)
## 5. Splice canonical ±1/2, INTRON_POSITION WRONG/NA → PVS1 via HGVSc canon
## 6. Multi-term CONSEQUENCE with splice_donor_variant in it → PVS1
## 7. Exon-intron junction spanning deletion (HGVSc exon_pos_intron) → PVS1
## 8. Intron-only positional deletion spanning exon boundary → PVS1
## 9. Splice canonical, last intron             → PVS1_STR (not PVS1)
## 10. Inframe deletion with splice consequence → NOT PVS1 (span excluded for inframe)
## 11. Non-LoF gene                             → nothing triggered
## 12. Non-relevant transcript                  → nothing triggered
## 13. start_lost consequence                   → PVS1_MOD
## 14. Donor +5 intronic variant, MaxEntScan Strong/Moderate → PVS1_MOD
##     (Weak MaxEntScan tier → not PVS1_MOD)
## 15. Splice outside ±2, no HGVSc signal,
##      no span                                 → nothing triggered
##
## Note: EXON_INTRON_JUNCTION_SPAN is set upstream (pcgr/annoutils.py) for
## indels spanning an exon-intron junction; tests 6, 7 and 16 (row 1) mimic
## this.

# ---------------------------------------------------------------------------
# Helper: build a minimal var_calls data.frame
# ---------------------------------------------------------------------------
make_var <- function(
    hgvsc               = NA_character_,
    consequence         = "splice_donor_variant",
    null_variant        = FALSE,
    loss_of_function    = TRUE,
    lof_filter          = NA_character_,
    nmd                 = NA_character_,
    intron_position     = NA_integer_,
    last_intron         = FALSE,
    exonic_status       = "exonic",
    cpg_mod             = "LoF",
    mane_select         = "NM_001234.5",
    mane_select2        = NA_character_,
    mane_plus_clinical  = NA_character_,
    mane_plus_clinical2 = NA_character_,
    refseq_id           = NA_character_,
    protein_rel_pos     = 0.5,
    pfam_domain         = NA_character_,
    maxentscan          = NA_character_,
    exon_intron_junction_span = FALSE
) {
  data.frame(
    VAR_ID                    = "v1",
    HGVSc                     = hgvsc,
    CONSEQUENCE               = consequence,
    NULL_VARIANT              = null_variant,
    LOSS_OF_FUNCTION          = loss_of_function,
    LOF_FILTER                = lof_filter,
    NMD                       = nmd,
    INTRON_POSITION           = intron_position,
    LAST_INTRON               = last_intron,
    EXONIC_STATUS             = exonic_status,
    CPG_MOD                   = cpg_mod,
    MANE_SELECT               = mane_select,
    MANE_SELECT2              = mane_select2,
    MANE_PLUS_CLINICAL        = mane_plus_clinical,
    MANE_PLUS_CLINICAL2       = mane_plus_clinical2,
    REFSEQ_TRANSCRIPT_ID      = ifelse(is.na(refseq_id), NA_character_, refseq_id),
    PROTEIN_RELATIVE_POSITION = protein_rel_pos,
    PFAM_DOMAIN_NAME          = pfam_domain,
    MAXENTSCAN                = maxentscan,
    EXON_INTRON_JUNCTION_SPAN = exon_intron_junction_span,
    stringsAsFactors          = FALSE
  )
}

# ---------------------------------------------------------------------------
# 1. NULL variant (frameshift / stop-gain), no NMD escape → PVS1
# ---------------------------------------------------------------------------
test_that("NULL_VARIANT with no NMD escape → ACMG_PVS1", {
  v <- make_var(
    hgvsc        = "ENST00000375499.8:c.500del",
    consequence  = "frameshift_variant",
    null_variant = TRUE
  )
  out <- cpsr::assign_PVS1_evidence(v)
  expect_true(out$ACMG_PVS1)
  expect_false(out$ACMG_PVS1_STR)
  expect_false(out$ACMG_PVS1_MOD)
})

# ---------------------------------------------------------------------------
# 2. NULL variant + NMD-escaping + PFAM domain → PVS1_STR
# ---------------------------------------------------------------------------
test_that("NULL_VARIANT + NMD_escaping + PFAM domain → ACMG_PVS1_STR", {
  v <- make_var(
    hgvsc        = "ENST00000375499.8:c.2900del",
    consequence  = "frameshift_variant",
    null_variant = TRUE,
    nmd          = "NMD_escaping_variant",
    pfam_domain  = "Helicase_C"
  )
  out <- cpsr::assign_PVS1_evidence(v)
  expect_false(out$ACMG_PVS1)
  expect_true(out$ACMG_PVS1_STR)
  expect_false(out$ACMG_PVS1_MOD)
})

# ---------------------------------------------------------------------------
# 3. NULL variant + NMD-escaping, no domain, late position → PVS1_MOD
# ---------------------------------------------------------------------------
test_that("NULL_VARIANT + NMD_escaping + no domain + late pos → ACMG_PVS1_MOD", {
  v <- make_var(
    hgvsc           = "ENST00000375499.8:c.2900del",
    consequence     = "frameshift_variant",
    null_variant    = TRUE,
    nmd             = "NMD_escaping_variant",
    pfam_domain     = NA_character_,
    protein_rel_pos = 0.95
  )
  out <- cpsr::assign_PVS1_evidence(v)
  expect_false(out$ACMG_PVS1)
  expect_false(out$ACMG_PVS1_STR)
  expect_true(out$ACMG_PVS1_MOD)
})

# ---------------------------------------------------------------------------
# 4. Canonical splice +1, INTRON_POSITION correct → PVS1 (non-last intron)
# ---------------------------------------------------------------------------
test_that("Canonical +1 splice with correct INTRON_POSITION → ACMG_PVS1", {
  v <- make_var(
    hgvsc           = "ENST00000375499.8:c.849_849+1del",
    consequence     = "coding_sequence_variant&splice_donor_variant",
    null_variant    = FALSE,
    intron_position = 1L,
    last_intron     = FALSE
  )
  out <- cpsr::assign_PVS1_evidence(v)
  expect_true(out$ACMG_PVS1)
  expect_false(out$ACMG_PVS1_STR)
})

# ---------------------------------------------------------------------------
# 5. Canonical +1 splice, INTRON_POSITION NA (buggy upstream) → PVS1 via HGVSc
# ---------------------------------------------------------------------------
test_that("Canonical +1 splice, INTRON_POSITION NA → PVS1 rescued by HGVSc", {
  v <- make_var(
    hgvsc           = "ENST00000456914.7:c.849_849+1del",
    consequence     = "coding_sequence_variant&splice_donor_variant",
    null_variant    = FALSE,
    intron_position = NA_integer_,   # buggy Python output
    last_intron     = FALSE
  )
  out <- cpsr::assign_PVS1_evidence(v)
  expect_true(out$ACMG_PVS1)
})

# ---------------------------------------------------------------------------
# 6. Multi-term CONSEQUENCE containing splice_donor_variant +
#    HGVSc exon→intron junction del → PVS1
#    (Row 1 in the missed TSV: c.743_765+36del)
# ---------------------------------------------------------------------------
test_that("Multi-term CONSEQUENCE + exon-to-intron junction del → PVS1", {
  v <- make_var(
    hgvsc = "ENST00000375499.8:c.743_765+36del",
    consequence =
      "coding_sequence_variant&intron_variant&splice_donor_5th_base_variant&splice_donor_variant",
    null_variant    = FALSE,
    intron_position = 36L,   # not canonical ±1/2 → previously missed
    last_intron     = FALSE,
    exon_intron_junction_span = TRUE   # set upstream for junction-spanning indels
  )
  out <- cpsr::assign_PVS1_evidence(v)
  expect_true(out$ACMG_PVS1)
})

# ---------------------------------------------------------------------------
# 7. Exon+intron spanning del with splice_acceptor in consequence
#    (Row 6: c.1103-32_1115del)
# ---------------------------------------------------------------------------
test_that("Splice_acceptor + intron-to-exon span del → PVS1", {
  v <- make_var(
    hgvsc = "ENST00000456914.7:c.1103-32_1115del",
    consequence =
      "coding_sequence_variant&intron_variant&splice_acceptor_variant",
    null_variant    = FALSE,
    intron_position = 0L,   # position 0 was not caught by old code
    last_intron     = FALSE,
    exon_intron_junction_span = TRUE   # set upstream for junction-spanning indels
  )
  out <- cpsr::assign_PVS1_evidence(v)
  expect_true(out$ACMG_PVS1)
})

# ---------------------------------------------------------------------------
# 8. Intronic-only del spanning into the next exon via HGVSc
#    (Row 15: c.555+1_555+9del — pure intron but disrupts canonical site)
# ---------------------------------------------------------------------------
test_that("Intronic del spanning canonical donor (+1_+9) → PVS1", {
  v <- make_var(
    hgvsc = "ENST00000366560.4:c.555+1_555+9del",
    consequence =
      "intron_variant&splice_donor_5th_base_variant&splice_donor_variant",
    null_variant    = FALSE,
    intron_position = 9L,   # not ±1/2 as returned by upstream tool
    last_intron     = FALSE
  )
  out <- cpsr::assign_PVS1_evidence(v)
  expect_true(out$ACMG_PVS1)
})

# ---------------------------------------------------------------------------
# 9. Canonical -2 splice, last intron → PVS1_STR (not PVS1)
# ---------------------------------------------------------------------------
test_that("Canonical -2 splice in last intron → ACMG_PVS1_STR only", {
  v <- make_var(
    hgvsc           = "ENST00000366560.4:c.1391-2del",
    consequence     = "splice_acceptor_variant",
    null_variant    = FALSE,
    intron_position = -2L,
    last_intron     = TRUE
  )
  out <- cpsr::assign_PVS1_evidence(v)
  expect_false(out$ACMG_PVS1)
  expect_true(out$ACMG_PVS1_STR)
})

# ---------------------------------------------------------------------------
# 10. Inframe deletion with a splice consequence → span excluded → not PVS1
#     (inframe_deletion precludes PVS1 via the span route)
# ---------------------------------------------------------------------------
test_that("Inframe del with splice consequence → ACMG_PVS1 FALSE", {
  v <- make_var(
    hgvsc = "ENST00000675843.1:c.8851_9171del",
    consequence =
      "inframe_deletion&intron_variant&splice_acceptor_variant&splice_donor_5th_base_variant&splice_donor_variant",
    null_variant    = FALSE,
    intron_position = 0L,
    last_intron     = FALSE
  )
  out <- cpsr::assign_PVS1_evidence(v)
  expect_false(out$ACMG_PVS1)
})

# ---------------------------------------------------------------------------
# 11. Non-LoF gene (CPG_MOD != "LoF") → nothing triggered
# ---------------------------------------------------------------------------
test_that("Non-LoF gene → all PVS1 flags FALSE", {
  v <- make_var(
    hgvsc        = "ENST00000375499.8:c.500_500+1del",
    consequence  = "coding_sequence_variant&splice_donor_variant",
    null_variant = FALSE,
    cpg_mod      = "GoF"
  )
  out <- cpsr::assign_PVS1_evidence(v)
  expect_false(out$ACMG_PVS1)
  expect_false(out$ACMG_PVS1_STR)
  expect_false(out$ACMG_PVS1_MOD)
})

# ---------------------------------------------------------------------------
# 12. Non-relevant transcript (no MANE, not in curated set) → nothing
# ---------------------------------------------------------------------------
test_that("Non-relevant transcript → all PVS1 flags FALSE", {
  v <- make_var(
    hgvsc               = "ENST00000999999.1:c.100_100+1del",
    consequence         = "coding_sequence_variant&splice_donor_variant",
    null_variant        = FALSE,
    mane_select         = NA_character_,
    mane_select2        = NA_character_,
    mane_plus_clinical  = NA_character_,
    mane_plus_clinical2 = NA_character_,
    refseq_id           = NA_character_
  )
  out <- cpsr::assign_PVS1_evidence(v)
  expect_false(out$ACMG_PVS1)
  expect_false(out$ACMG_PVS1_STR)
  expect_false(out$ACMG_PVS1_MOD)
})

# ---------------------------------------------------------------------------
# 13. start_lost consequence → PVS1_MOD
# ---------------------------------------------------------------------------
test_that("start_lost → ACMG_PVS1_MOD", {
  v <- make_var(
    hgvsc       = "ENST00000375499.8:c.1A>G",
    consequence = "start_lost",
    null_variant = FALSE
  )
  out <- cpsr::assign_PVS1_evidence(v)
  expect_true(out$ACMG_PVS1_MOD)
})

# ---------------------------------------------------------------------------
# 14. intron_variant&splice_donor_5th_base_variant (no splice_donor_variant)
#     → PVS1_MOD when MaxEntScan predicts a Strong/Moderate donor +5 effect
#     (MAXENTSCAN: '<score>|<stratum>|<tier>'); Weak tier → not PVS1_MOD
# ---------------------------------------------------------------------------
test_that("Donor +5 intronic variant with MaxEntScan Strong → ACMG_PVS1_MOD", {
  v <- make_var(
    hgvsc       = "ENST00000375499.8:c.267+5G>A",
    consequence = "intron_variant&splice_donor_5th_base_variant",
    null_variant = FALSE,
    intron_position = 5L,
    maxentscan  = "-8.5|Donor_+5|Strong"
  )
  out <- cpsr::assign_PVS1_evidence(v)
  expect_true(out$ACMG_PVS1_MOD)
  expect_false(out$ACMG_PVS1)
})

test_that("Donor +5 intronic variant with MaxEntScan Weak → no ACMG_PVS1_MOD", {
  v <- make_var(
    hgvsc       = "ENST00000375499.8:c.267+5G>A",
    consequence = "intron_variant&splice_donor_5th_base_variant",
    null_variant = FALSE,
    intron_position = 5L,
    maxentscan  = "-1.2|Donor_+5|Weak"
  )
  out <- cpsr::assign_PVS1_evidence(v)
  expect_false(out$ACMG_PVS1_MOD)
})

# ---------------------------------------------------------------------------
# 15. Deep intronic splice-region variant, no HGVSc span, no canonical pos
#     → all flags FALSE
# ---------------------------------------------------------------------------
test_that("Deep intronic splice_region_variant → all PVS1 FALSE", {
  v <- make_var(
    hgvsc           = "ENST00000375499.8:c.500+50G>T",
    consequence     = "intron_variant&splice_region_variant",
    null_variant    = FALSE,
    intron_position = 50L,
    last_intron     = FALSE
  )
  out <- cpsr::assign_PVS1_evidence(v)
  expect_false(out$ACMG_PVS1)
  expect_false(out$ACMG_PVS1_STR)
  expect_false(out$ACMG_PVS1_MOD)
})

# ---------------------------------------------------------------------------
# 16. Vectorised: multiple rows processed correctly in a single call
# ---------------------------------------------------------------------------
test_that("Vectorised call assigns PVS1 correctly per row", {
  rows <- rbind(
    ## row 1 – should fire PVS1
    make_var(
      hgvsc           = "ENST00000375499.8:c.761_765+33del",
      consequence     = "coding_sequence_variant&intron_variant&splice_donor_5th_base_variant&splice_donor_variant",
      null_variant    = FALSE,
      intron_position = 33L,
      last_intron     = FALSE,
      exon_intron_junction_span = TRUE
    ),
    ## row 2 – deep intronic, should NOT fire
    make_var(
      hgvsc           = "ENST00000375499.8:c.500+80T>A",
      consequence     = "intron_variant",
      null_variant    = FALSE,
      intron_position = 80L,
      last_intron     = FALSE
    )
  )
  out <- cpsr::assign_PVS1_evidence(rows)
  expect_equal(nrow(out), 2L)
  expect_true(out$ACMG_PVS1[1])
  expect_false(out$ACMG_PVS1[2])
})

# ---------------------------------------------------------------------------
# Helper columns must NOT appear in returned data.frame
# ---------------------------------------------------------------------------
test_that("Temporary PVS1 helper columns cleaned up on return", {
  v <- make_var(
    hgvsc           = "ENST00000375499.8:c.500_500+1del",
    consequence     = "coding_sequence_variant&splice_donor_variant",
    null_variant    = FALSE,
    intron_position = 1L
  )
  out <- cpsr::assign_PVS1_evidence(v)
  helper_cols <- c(
    "tmp_MES", "MES_STRATUM", "MES_TIER",
    "HGVSC_LOCAL", "PVS1_RELEVANT_TRANSCRIPT", "PVS1_LOF_GENE",
    "PVS1_PURELY_INTRONIC", "PVS1_SPLICE_CONSEQUENCE",
    "PVS1_INFRAME_CONSEQUENCE", "PVS1_HGVSC_CANONICAL_SITE",
    "PVS1_SPLICE_HIGH_IMPACT"
  )
  expect_false(any(helper_cols %in% colnames(out)))
})
