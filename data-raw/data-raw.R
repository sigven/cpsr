
insilico_path_predictors <-
  c("DBNSFP_ALPHA_MISSENSE",
    "DBNSFP_BAYESDEL_ADDAF",
    "DBNSFP_CLINPRED",
    "DBNSFP_DEOGEN2",
    "DBNSFP_ESM1B",
    "DBNSFP_FATHMM_XF",
    "DBNSFP_LIST_S2",
    "DBNSFP_M_CAP",
    "DBNSFP_META_RNN",
    "DBNSFP_MUTFORMER",
    "DBNSFP_MUTATIONASSESSOR",
    "DBNSFP_MUTATIONTASTER",
    "DBNSFP_PHACTBOOST",
    "DBNSFP_PRIMATEAI",
    "DBNSFP_PROVEAN",
    "DBNSFP_SIFT",
    "DBNSFP_POLYPHEN2_HVAR",
    "DBNSFP_SPLICE_SITE_ADA",
    "DBNSFP_SPLICE_SITE_RF")

usethis::use_data(insilico_path_predictors, overwrite = T)

col_format_output <- list()

###--- column formats ----####

###--- HTML tier ---####
col_format_output[['report_tbl_classification']] <-
  c(
    "SYMBOL",
    "GENENAME",
    "CONSEQUENCE",
    "ALTERATION",
    "GENOTYPE",
    "CLASSIFICATION",
    "PROTEIN_DOMAIN",
    "DP_CONTROL",
    "PROTEIN_CHANGE",
    "HGVSp",
    "HGVSc",
    "CODING_STATUS",
    "HGVSc_RefSeq",
    "ENSEMBL_GENE_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "ONCOGENE",
    "TUMOR_SUPPRESSOR",
    "CDS_CHANGE",
    "MUTATION_HOTSPOT",
    "RMSK_HIT",
    "PREDICTED_EFFECT",
    "SPLICE_EFFECT",
    "miRNA_TARGET_HIT",
    "miRNA_TARGET_HIT_PREDICTION",
    "TF_BINDING_SITE_VARIANT",
    "TF_BINDING_SITE_VARIANT_INFO",
    "GERP_SCORE",
    "LOSS_OF_FUNCTION",
    "NMD",
    "EXON_INTRON_JUNCTION_SPAN",
    "LOF_FILTER",
    "DBSNP_RSID",
    "CLINVAR",
    "CLINVAR_CLASSIFICATION",
    "CLINVAR_GOLD_STARS",
    "CLINVAR_CONFLICTED",
    "CLINVAR_VARIANT_ORIGIN",
    "CLINVAR_PHENOTYPE",
    "ASSERTION_AUTHORITY",
    "ASSERTION_RATIONALE",
    "CPSR_CLASSIFICATION",
    "CPSR_PATHOGENICITY_SCORE",
    "ACMG_DOC",
    "ACMG_CODE",
    "gnomADe_AF",
    "gnomADg_AF",
    "GENOMIC_CHANGE",
    "GENOME_VERSION"
  )

## define tags/variables to display in data tables (secondary findings)
col_format_output[['report_tbl_sf']] <-
  c(
    "SYMBOL",
    "ALTERATION",
    "CLINVAR_CLASSIFICATION",
    "CLINVAR_PHENOTYPE",
    "GENOTYPE",
    "CONSEQUENCE",
    "GENENAME",
    "PROTEIN_DOMAIN",
    "DP_CONTROL",
    "PROTEIN_CHANGE",
    "HGVSp",
    "HGVSc",
    "HGVSc_RefSeq",
    "ENSEMBL_GENE_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "CDS_CHANGE",
    "PREDICTED_EFFECT",
    "SPLICE_EFFECT",
    "LOSS_OF_FUNCTION",
    "LOF_FILTER",
    "DBSNP_RSID",
    "CLINVAR",
    "CLINVAR_GOLD_STARS",
    "CLINVAR_CONFLICTED",
    "ASSERTION_AUTHORITY",
    "gnomADe_AF",
    "gnomADg_AF",
    "GENOMIC_CHANGE",
    "GENOME_VERSION"
  )

## define tags/variables to display in data tables (PGx findings)
col_format_output[['report_tbl_pgx']] <-
  c(
    "SYMBOL",
    "ALTERATION",
    "CLINVAR_CLASSIFICATION",
    "CLINVAR_PHENOTYPE",
    "GENOTYPE",
    "CONSEQUENCE",
    "GENENAME",
    "PROTEIN_DOMAIN",
    "DP_CONTROL",
    "PROTEIN_CHANGE",
    "HGVSp",
    "HGVSc",
    "HGVSc_RefSeq",
    "ENSEMBL_GENE_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "CDS_CHANGE",
    "PREDICTED_EFFECT",
    "SPLICE_EFFECT",
    "LOSS_OF_FUNCTION",
    "LOF_FILTER",
    "DBSNP_RSID",
    "CLINVAR",
    "CLINVAR_GOLD_STARS",
    "CLINVAR_CONFLICTED",
    "gnomADe_AF",
    "gnomADg_AF",
    "GENOMIC_CHANGE",
    "GENOME_VERSION"
  )


## define tags/variables to display in data tables (GWAS findings)
col_format_output[['report_tbl_gwas']] <-
  c(
    "SYMBOL",
    "CONSEQUENCE",
    "GWAS_CITATION",
    "ALTERATION",
    "GENOTYPE",
    "GENENAME",
    "LOSS_OF_FUNCTION",
    "PROTEIN_DOMAIN",
    "DP_CONTROL",
    "PROTEIN_CHANGE",
    "GWAS_PHENOTYPE",
    "HGVSp",
    "HGVSc",
    "HGVSc_RefSeq",
    "ENSEMBL_GENE_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "CDS_CHANGE",
    "CODING_STATUS",
    "miRNA_TARGET_HIT",
    "miRNA_TARGET_HIT_PREDICTION",
    "TF_BINDING_SITE_VARIANT",
    "TF_BINDING_SITE_VARIANT_INFO",
    "GERP_SCORE",
    "PREDICTED_EFFECT",
    "DBSNP_RSID",
    "gnomADe_AF",
    "gnomADg_AF",
    "GENOMIC_CHANGE",
    "GENOME_VERSION"
  )

## define tags/variables to display in output TSV
col_format_output[['tsv']] <-
  c("SAMPLE_ID",
    "GENOMIC_CHANGE",
    "VAR_ID",
    "GENOME_VERSION",
    "GENOTYPE",
    "DP_CONTROL",
    "ASSERTION_AUTHORITY",
    "ASSERTION_RATIONALE",
    "VARIANT_CLASS",
    "CODING_STATUS",
    "SYMBOL",
    "GENENAME",
    "CCDS",
    "ENTREZGENE",
    "UNIPROT_ID",
    "ENSEMBL_GENE_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "MANE_SELECT",
    "MANE_SELECT2",
    "ONCOGENE",
    "TUMOR_SUPPRESSOR",
    "CPG_MOD",
    "CPG_MOI",
    "CONSEQUENCE",
    "ALTERATION",
    "PROTEIN_CHANGE",
    "PFAM_DOMAIN",
    "PFAM_DOMAIN_NAME",
    "HGVSp",
    "GRANTHAM_DISTANCE",
    "HGVSc",
    "HGVSc_RefSeq",
    "CDS_CHANGE",
    "LAST_EXON",
    "LAST_INTRON",
    "EXON",
    "EXON_AFFECTED",
    "EXON_POSITION",
    "INTRON_POSITION",
    "NMD",
    "EXON_INTRON_JUNCTION_SPAN",
    "EXONIC_STATUS",
    "PROTEIN_RELATIVE_POSITION",
    "MUTATION_HOTSPOT",
    "RMSK_HIT",
    "EFFECT_PREDICTIONS",
    "MAXENTSCAN",
    "SPLICE_EFFECT",
    "LOSS_OF_FUNCTION",
    "LOF_FILTER",
    "NULL_VARIANT",
    "DBMTS",
    "REGULATORY_ANNOTATION",
    "TF_BINDING_SITE_VARIANT",
    "TF_BINDING_SITE_VARIANT_INFO",
    "VEP_ALL_CSQ",
    "GERP_SCORE",
    "DBSNP_RSID",
    "CLINVAR_CLASSIFICATION",
    "CLINVAR_MSID",
    "CLINVAR_VARIANT_ORIGIN",
    "CLINVAR_CONFLICTED",
    "CLINVAR_PHENOTYPE",
    "CLINVAR_PHENOTYPE_CANCER",
    "CLINVAR_GOLD_STARS",
    "N_INSILICO_CALLED",
    "N_INSILICO_DAMAGING",
    "N_INSILICO_TOLERATED",
    "N_INSILICO_SPLICING_NEUTRAL",
    "N_INSILICO_SPLICING_AFFECTED",
    "gnomADe_AF",
    "gnomADg_AF",
    "CLASSIFICATION",
    "CPSR_CLASSIFICATION",
    "CPSR_PATHOGENICITY_SCORE",
    "ACMG_CODE"
  )

col_format_output[['report_tbl_biomarker']] <-
  c('VAR_ID',
    'VARIANT_CLASS',
    'ENTREZGENE',
    'GENOTYPE',
    'SYMBOL',
    'ALTERATION',
    'SAMPLE_ALTERATION',
    'CONSEQUENCE',
    'PROTEIN_CHANGE',
    'ASSERTION_AUTHORITY',
    'ASSERTION_RATIONALE',
    'CLASSIFICATION',
    'BM_CANCER_TYPE',
    'BM_EVIDENCE_LEVEL',
    'BM_CONTEXT',
    'BM_MOLECULAR_PROFILE',
    'BM_REFERENCE',
    'BM_EVIDENCE_DESCRIPTION',
    'BM_CLINICAL_SIGNIFICANCE',
    'BM_EVIDENCE_TYPE',
    'BM_THERAPEUTIC_CONTEXT',
    'BM_RATING',
    'BM_EVIDENCE_DIRECTION',
    'BM_SOURCE_DB',
    'BM_EVIDENCE_ID',
    'BM_DISEASE_ONTOLOGY_ID',
    'BM_PRIMARY_SITE',
    'BM_MAPPING_CONFIDENCE',
    'BM_RESOLUTION')

## define tags/variables to display in output Excel
col_format_output[['xlsx_panel']] <-
  c("SYMBOL",
    "ENSEMBL_GENE_ID",
    "ENTREZGENE",
    "GENE_BIOTYPE",
    "PRIMARY_TARGET",
    "MOI",
    "MOD",
    "ID",
    "PANEL_NAME",
    "PANEL_URL",
    "PANEL_VERSION",
    "CONFIDENCE_LEVEL",
    "CPG_SOURCE"
    )

## define tags/variables to display in output Excel
col_format_output[['xlsx_classification']] <-
  c("SAMPLE_ID",
    "GENOMIC_CHANGE",
    "GENOTYPE",
    "DP_CONTROL",
    "GENOME_VERSION",
    "VARIANT_CLASS",
    "CODING_STATUS",
    "SYMBOL",
    "GENENAME",
    "CONSEQUENCE",
    "ALTERATION",
    "PROTEIN_CHANGE",
    "CLASSIFICATION",
    "ASSERTION_AUTHORITY",
    "ASSERTION_RATIONALE",
    "CLINVAR_CLASSIFICATION",
    "CPSR_CLASSIFICATION",
    "ACMG_CODE",
    "CPSR_PATHOGENICITY_SCORE",
    "ENTREZGENE",
    "ENSEMBL_GENE_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "ONCOGENE",
    "TUMOR_SUPPRESSOR",
    "CPG_MOD",
    "CPG_MOI",
    "PFAM_DOMAIN_NAME",
    "HGVSp",
    "GRANTHAM_DISTANCE",
    "HGVSc",
    "HGVSc_RefSeq",
    "CDS_CHANGE",
    "MUTATION_HOTSPOT",
    "EFFECT_PREDICTIONS",
    "SPLICE_EFFECT",
    "MAXENTSCAN",
    "LOSS_OF_FUNCTION",
    "LOF_FILTER",
    "NULL_VARIANT",
    "DBSNP_RSID",
    "CLINVAR_MSID",
    "CLINVAR_VARIANT_ORIGIN",
    "CLINVAR_CONFLICTED",
    "CLINVAR_PHENOTYPE",
    "CLINVAR_GOLD_STARS",
    "gnomADe_AF",
    "gnomADg_AF"
  )

## define tags/variables to display in output Excel
col_format_output[['xlsx_secondary']] <-
  c("SAMPLE_ID",
    "GENOMIC_CHANGE",
    "GENOTYPE",
    "DP_CONTROL",
    "GENOME_VERSION",
    "VARIANT_CLASS",
    "SYMBOL",
    "GENENAME",
    "CONSEQUENCE",
    "PROTEIN_CHANGE",
    "ASSERTION_AUTHORITY",
    "ASSERTION_RATIONALE",
    "CLASSIFICATION",
    "CPSR_CLASSIFICATION",
    "CPSR_PATHOGENICITY_SCORE",
    "ACMG_CODE",
    "CLINVAR_CLASSIFICATION",
    "CLINVAR_MSID",
    "CLINVAR_VARIANT_ORIGIN",
    "CLINVAR_CONFLICTED",
    "CLINVAR_PHENOTYPE",
    "CLINVAR_GOLD_STARS",
    "ENSEMBL_GENE_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "PFAM_DOMAIN_NAME",
    "HGVSp",
    "HGVSc",
    "HGVSc_RefSeq",
    "CDS_CHANGE",
    "CODING_STATUS",
    "MUTATION_HOTSPOT",
    "EFFECT_PREDICTIONS",
    "SPLICE_EFFECT",
    "MAXENTSCAN",
    "LOSS_OF_FUNCTION",
    "LOF_FILTER",
    "NULL_VARIANT",
    "DBSNP_RSID",
    "gnomADe_AF",
    "gnomADg_AF"
  )

## define tags/variables to display in output Excel
col_format_output[['xlsx_pgx']] <-
  c("SAMPLE_ID",
    "GENOMIC_CHANGE",
    "GENOTYPE",
    "DP_CONTROL",
    "GENOME_VERSION",
    "VARIANT_CLASS",
    "SYMBOL",
    "GENENAME",
    "CONSEQUENCE",
    "PROTEIN_CHANGE",
    "CLINVAR_CLASSIFICATION",
    "CLINVAR_MSID",
    "CLINVAR_VARIANT_ORIGIN",
    "CLINVAR_CONFLICTED",
    "CLINVAR_PHENOTYPE",
    "CLINVAR_GOLD_STARS",
    "ENSEMBL_GENE_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "PFAM_DOMAIN_NAME",
    "HGVSp",
    "HGVSc",
    "HGVSc_RefSeq",
    "CDS_CHANGE",
    "CODING_STATUS",
    "MUTATION_HOTSPOT",
    "EFFECT_PREDICTIONS",
    "SPLICE_EFFECT",
    "LOSS_OF_FUNCTION",
    "LOF_FILTER",
    "NULL_VARIANT",
    "DBSNP_RSID",
    "gnomADe_AF",
    "gnomADg_AF"
  )


## define tags/variables to display in output Excel
col_format_output[['xlsx_biomarker']] <-
  c("SAMPLE_ID",
    "GENOMIC_CHANGE",
    "GENOTYPE",
    "DP_CONTROL",
    "GENOME_VERSION",
    "VARIANT_CLASS",
    "SYMBOL",
    "GENENAME",
    "CONSEQUENCE",
    "PROTEIN_CHANGE",
    "ASSERTION_AUTHORITY",
    "ASSERTION_RATIONALE",
    "CLASSIFICATION",
    "CPSR_CLASSIFICATION",
    "CPSR_PATHOGENICITY_SCORE",
    "ACMG_CODE",
    "CLINVAR_CLASSIFICATION",
    "BM_CANCER_TYPE",
    "BM_DISEASE_ONTOLOGY_ID",
    "BM_PRIMARY_SITE",
    "BM_CLINICAL_SIGNIFICANCE",
    "BM_THERAPEUTIC_CONTEXT",
    "BM_CITATION",
    "BM_RATING",
    "BM_MOLECULAR_PROFILE",
    "BM_EVIDENCE_TYPE",
    "BM_EVIDENCE_LEVEL",
    "BM_EVIDENCE_DIRECTION",
    "BM_EVIDENCE_DESCRIPTION",
    "BM_SOURCE_DB",
    "BM_EVIDENCE_ID",
    "BM_VARIANT_ORIGIN",
    "BM_MATCH",
    "BM_RESOLUTION"
  )


#---- acmg ----#
acmg <- list()
acmg[['gnomAD_pops']] <-
  c("GLOBAL","NFE","AMR","AFR","SAS","EAS","FIN")

acmg[["score2tier"]] <- data.frame()
acmg[["evidence_codes"]] <-
  utils::read.table(file = "data-raw/acmg_evidence.tsv",
                    header = T, stringsAsFactors = F,
                    comment.char = "", na.strings = c("NA"),
                    sep = "\t")
acmg[["pathogenic_range_gnomad"]] <- list()
acmg[["pathogenic_range_gnomad"]][["af"]] <- 0.00005
acmg[["pathogenic_range_gnomad"]][["min_an"]] <- 4000
acmg[["insilico_pred_min_majority"]] <- 8
acmg[["insilico_pred_max_minority"]] <- 2

acmg[['score_thresholds']] <- list()
acmg[['score_thresholds']][['p_lower']] <- 4
acmg[['score_thresholds']][['lp_upper']] <- 3
acmg[['score_thresholds']][['lp_lower']] <- 2
acmg[['score_thresholds']][['vus_upper']] <- 1
acmg[['score_thresholds']][['vus_lower']] <- 0
acmg[['score_thresholds']][['lb_upper']] <- -1
acmg[['score_thresholds']][['lb_lower']] <- -3
acmg[['score_thresholds']][['b_upper']] <- -4


acmg[["score2tier"]] <-
  data.frame("CPSR_CLASSIFICATION" = "Pathogenic",
             "CPSR_PATHOGENICITY_SCORE" = "<b>[4, ]</b>")
acmg[["score2tier"]] <-
  dplyr::bind_rows(
    acmg[["score2tier"]],
    data.frame("CPSR_CLASSIFICATION" = "Likely Pathogenic",
               "CPSR_PATHOGENICITY_SCORE" = "<b>[2, 3]</b>"))
acmg[["score2tier"]] <-
  dplyr::bind_rows(
    acmg[["score2tier"]],
    data.frame("CPSR_CLASSIFICATION" = "VUS",
               "CPSR_PATHOGENICITY_SCORE" = "<b>[0, 1]</b>"))
acmg[["score2tier"]] <-
  dplyr::bind_rows(
    acmg[["score2tier"]],
    data.frame("CPSR_CLASSIFICATION" = "Likely Benign",
               "CPSR_PATHOGENICITY_SCORE" = "<b>[-3, -1]</b>"))
acmg[["score2tier"]] <-
  dplyr::bind_rows(
    acmg[["score2tier"]],
    data.frame("CPSR_CLASSIFICATION" = "Benign",
               "CPSR_PATHOGENICITY_SCORE" = "<b>[, -4]</b>"))


color_palette <- list()
color_palette[['report']] <- "#007a74"
color_palette[['none']] <- "#8B8989"
color_palette[['warning']] <- "#FFA500"
color_palette[['genotypes']] <- list()
color_palette[['genotypes']][['bgcolor_values']] <-
  c('#bdbdbd','#bdbdbd','#007a74','#004c47')
color_palette[['genotypes']][['color_values']] <-
  c('#2c313c','#2c313c','#ffffff','#ffffff')
color_palette[['genotypes']][['levels']] <-
  c('undefined','hom_ref','het','hom_alt')
color_palette[['review_status_clinvar']] <- list()
color_palette[['review_status_clinvar']][['bgcolor_values']] <-
  c('#ffffff','#73bbb4','#469a94','#007a74','#004c47')
color_palette[['review_status_clinvar']][['color_values']] <-
  c('#000000','#ffffff','#ffffff','#ffffff','#ffffff')
color_palette[['review_status_clinvar']][['levels']] <-
  c(0, 1, 2, 3, 4)

color_palette[['pathogenicity']] <- list()
color_palette[['pathogenicity']][['levels']] <-
  c("Benign",
    "Likely Benign",
    "VUS",
    "Likely Pathogenic",
    "Pathogenic")
color_palette[['pathogenicity']][['values']] <-
  c("#077009",
    "#6FB572",
    "#2c313c",
    "#9C3948",
    "#9E0142")
color_palette[['pathogenicity_score']] <- list()
color_palette[['pathogenicity_score']][['levels']] <- -9:9
color_palette[['pathogenicity_score']][['bgcolor_values']] <-
  c("#077009",
    "#208321",
    "#39963A",
    "#52A953",
    "#6BBC6C",
    "#6FB572",
    "#5E9464",
    "#4D7357",
    "#3C5249",
    "#2C313C",
    "#513340",
    "#773644",
    "#9C3948",
    "#C23C4C",
    "#CE374D",
    "#C2294A",
    "#B61C47",
    "#AA0E44",
    "#9E0142")

color_palette[['biomarker_types']] <- list()
color_palette[['biomarker_types']][['levels']] <-
  c("Sensitivity/Response",
    "Poor Outcome",
    "Better Outcome",
    "Resistance/Non-response",
    "Predisposition",
    "Diagnostic",
    "Toxicity",
    "Adverse Response")
color_palette[['biomarker_types']][['values']] <-
  c("#1565C0",
    "#E65100",
    "#2E7D32",
    "#C62828",
    color_palette[['report']],
    "#6A1B9A",
    "#868686",
    "#868686")



## Numeric scale
#scores <- -9:9

# Your 5 category colors (benign → pathogenic)
# category_cols <- c(
#   "Benign"            = "#077009",  # -4 to -9
#   "Likely Benign"     = "#78C679",  # -1 to -3
#   "VUS"               = "#2c313c",  # 0, 1
#   "Likely Pathogenic" = "#D53E4F",  # 2, 3
#   "Pathogenic"        = "#9E0142"   # 4 to 9
# )
#
# # Assign one anchor per category midpoint
# anchors <- c(-6.5, -2, 0.5, 2.5, 6.5)
#
# # Create color ramp
# ramp_fun <- colorRampPalette(category_cols)
#
# # Interpolate colors across full scale
# expanded_colors <- ramp_fun(length(scores))
#
# names(expanded_colors) <- scores
# expanded_colors <- unname(expanded_colors)

usethis::use_data(color_palette, overwrite = T)
usethis::use_data(acmg, overwrite = T)
usethis::use_data(col_format_output, overwrite = T)

# #---- create CPSR curated transcripts ----#
# MANE Select + MANE Plus Clinical NM_ accessions for all CPG genes (grch38
# bundle), plus manually curated older accessions used by VEP/ClinVar on
# GRCh37 where MANE annotations are absent.
curated_transcripts <- utils::read.table(
  file = "data-raw/curated_transcripts.tsv",
  header = TRUE, stringsAsFactors = FALSE,
  sep = "\t"
)
usethis::use_data(curated_transcripts, overwrite = T)

# my_log4r_layout <- function(level, ...) {
#   paste0(format(Sys.time()), " - cpsr-report-generation - ",
#          level, " - ", ..., "\n", collapse = "")
# }
#
# log4r_logger <-
#   log4r::logger(
#     threshold = "INFO", appenders = log4r::console_appender(my_log4r_layout))
#
# # this gets passed on to all the log4r_* functions inside the pkg
# options("PCGRR_LOG4R_LOGGER" = log4r_logger)
#
# panel_zero <- list()
# for(build in c('grch37','grch38')){
#   ref_data <- pcgrr::load_reference_data(
#     pcgr_db_assembly_dir =
#       file.path(
#         "/Users/sigven/project_data/data/data__pcgrdb/dev/pcgrdb",
#         "20260426/data",
#         build),
#     genome_assembly = build
#   )
#
#   set1 <- ref_data$gene$cpg |>
#     dplyr::filter(CPG_SOURCE != "ACMG_SF" &
#                     CPG_SOURCE != "CPIC_PGX_ONCOLOGY") |>
#     dplyr::filter(!is.na(ENSEMBL_GENE_ID)) |>
#     dplyr::inner_join(
#       dplyr::select(ref_data$gene$gene_xref,
#                     ENTREZGENE,
#                     ENSEMBL_GENE_ID,
#                     GENE_BIOTYPE,
#                     GENENAME,
#                     TSG,
#                     TSG_SUPPORT,
#                     ONCOGENE,
#                     ONCOGENE_SUPPORT),
#       by = c("ENTREZGENE","ENSEMBL_GENE_ID")
#     )
#
#   set2 <- ref_data$gene$cpg |>
#     dplyr::filter(CPG_SOURCE != "ACMG_SF" &
#                     CPG_SOURCE != "CPIC_PGX_ONCOLOGY") |>
#     dplyr::filter(is.na(ENSEMBL_GENE_ID)) |>
#     dplyr::select(-c("ENSEMBL_GENE_ID")) |>
#     dplyr::inner_join(
#       dplyr::select(ref_data$gene$gene_xref,
#                     ENTREZGENE,
#                     ENSEMBL_GENE_ID,
#                     GENE_BIOTYPE,
#                     GENENAME,
#                     TSG,
#                     TSG_SUPPORT,
#                     ONCOGENE,
#                     ONCOGENE_SUPPORT),
#       by = c("ENTREZGENE")
#     )
#
#   panel_zero[[build]] <- dplyr::bind_rows(set1, set2) |>
#     dplyr::mutate(
#       PANEL_NAME = "CPSR superpanel of cancer predisposition genes",
#       PANEL_VERSION = "v2026_05") |>
#     dplyr::rename(
#       TUMOR_SUPPRESSOR = TSG,
#       TUMOR_SUPPRESSOR_SUPPORT = TSG_SUPPORT
#     ) |>
#     dplyr::left_join(
#       dplyr::select(
#         dplyr::filter(
#           ref_data$variant$clinvar_gene_stats,
#           .data$CONFIDENCE == "min2goldstars"),
#         c("ENTREZGENE",
#           "N_TRUNC_PATH",
#           "N_NONTRUNC_PATH",
#           "N_MISSENSE_PATH",
#           "N_MISSENSE_BENIGN",
#           "BENIGN_MISSENSE_FRAC",
#           "PATH_TRUNC_FRAC")
#       ), by = "ENTREZGENE"
#     ) |>
#     dplyr::select(
#       dplyr::any_of(
#         c("ENTREZGENE",
#           "SYMBOL",
#           "GENENAME",
#           "GENE_BIOTYPE",
#           "ENSEMBL_GENE_ID",
#           "TUMOR_SUPPRESSOR",
#           "TUMOR_SUPPRESSOR_SUPPORT",
#           "ONCOGENE",
#           "ONCOGENE_SUPPORT",
#           "CPG_SOURCE",
#           "CPG_MOD",
#           "CPG_MOI",
#           "CPG_PHENOTYPES",
#           "CPG_CANCER_CUI",
#           "CPG_SYNDROME_CUI")
#       ),
#       dplyr::everything()
#     ) |>
#     dplyr::distinct()
# }
# #
# workbook <- openxlsx2::wb_workbook() |>
#   openxlsx2::wb_add_worksheet(sheet = "CPSR_SUPERPANEL.GRCH37") |>
#   openxlsx2::wb_add_worksheet(sheet = "CPSR_SUPERPANEL.GRCH38") |>
#   openxlsx2::wb_add_data_table(
#     sheet = "CPSR_SUPERPANEL.GRCH37",
#     x = panel_zero[['grch37']],
#     start_row = 1,
#     start_col = 1,
#     col_names = TRUE,
#     na.strings = "NA",
#     table_style = "TableStyleMedium15") |>
#   openxlsx2::wb_add_data_table(
#     sheet = "CPSR_SUPERPANEL.GRCH38",
#     x = panel_zero[['grch38']],
#     start_row = 1,
#     start_col = 1,
#     col_names = TRUE,
#     na.strings = "NA",
#     table_style = "TableStyleMedium16")
#
# openxlsx2::wb_save(
#   wb = workbook,
#   "pkgdown/assets/cpsr_superpanel_2026_04.xlsx",
#   overwrite = TRUE)
#
# #
# panel_zero_display <- panel_zero$grch38 |>
#   dplyr::select(
#     c("ENTREZGENE",
#       "SYMBOL",
#       "ENTREZGENE",
#       "ENSEMBL_GENE_ID",
#       "GENENAME",
#       "CPG_PHENOTYPES",
#       "CPG_MOI",
#       "CPG_MOD",
#       "CPG_SOURCE",
#     )
#   ) |>
#   dplyr::mutate(
#     CPG_SOURCE = stringr::str_replace_all(
#       CPG_SOURCE, "&", ", "
#     )) |>
#   dplyr::mutate(
#     CPG_SOURCE = stringr::str_replace_all(
#       CPG_SOURCE, "ACMG_SF", ""
#     )
#   ) |>
#   dplyr::mutate(
#     CPG_SOURCE = stringr::str_replace_all(
#       CPG_SOURCE, "^, ", ""
#     )
#   ) |>
#   dplyr::mutate(
#     GENE = paste0(
#       "<a href='https://www.ncbi.nlm.nih.gov/gene/",
#       .data$ENTREZGENE,
#       "' target='_blank'>",
#       .data$SYMBOL, "</a>"
#     )
#   ) |>
#   dplyr::select(
#     c("GENE","ENTREZGENE","ENSEMBL_GENE_ID",
#       "CPG_MOD", "CPG_MOI", "GENENAME",
#       "CPG_SOURCE", "CPG_PHENOTYPES")
#   )
#
# readr::write_tsv(
#   panel_zero_display, file = "inst/extdata/panel_zero.tsv.gz",
#   na = "NA", col_names = T,quote = "none"
# )
#
