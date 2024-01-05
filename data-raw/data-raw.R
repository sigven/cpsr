
col_format_output <- list()

col_format_output[['html_tier']] <-
  c(
    "SYMBOL",
    "CLINVAR_PHENOTYPE",
    "CONSEQUENCE",
    "PROTEIN_CHANGE",
    "GENOTYPE",
    "OFFICIAL_GENENAME",
    "PROTEIN_DOMAIN",
    "HGVSp",
    "HGVSc",
    "REFSEQ",
    "CDS_CHANGE",
    "MUTATION_HOTSPOT",
    "RMSK_HIT",
    "PREDICTED_EFFECT",
    "miRNA_TARGET_HIT",
    "miRNA_TARGET_HIT_PREDICTION",
    "TF_BINDING_SITE_VARIANT",
    "TF_BINDING_SITE_VARIANT_INFO",
    "GERP_SCORE",
    "LOSS_OF_FUNCTION",
    "DBSNP",
    "CLINVAR",
    "CLINVAR_CLASSIFICATION",
    "CLINVAR_REVIEW_STATUS_STARS",
    "CLINVAR_CONFLICTED",
    "CLINVAR_VARIANT_ORIGIN",
    "CPSR_CLASSIFICATION_SOURCE",
    "CPSR_CLASSIFICATION",
    "CPSR_PATHOGENICITY_SCORE",
    "CPSR_CLASSIFICATION_DOC",
    "CPSR_CLASSIFICATION_CODE",
    "ONCOGENE",
    "TUMOR_SUPPRESSOR",
    "CANCERGENE_EVIDENCE",
    "gnomADe_AF",
    "GENOMIC_CHANGE",
    "GENOME_VERSION"
  )

## define tags/variables to display in data tables (secondary findings)
col_format_output[['html_sf']] <-
  c(
    "SYMBOL",
    "CONSEQUENCE",
    "CLINVAR_CLASSIFICATION",
    "CLINVAR_PHENOTYPE",
    "PROTEIN_CHANGE",
    "GENOTYPE",
    "OFFICIAL_GENENAME",
    "PROTEIN_DOMAIN",
    "HGVSp",
    "HGVSc",
    "REFSEQ",
    "CDS_CHANGE",
    "REFSEQ_TRANSCRIPT_ID",
    "MUTATION_HOTSPOT",
    "RMSK_HIT",
    "PREDICTED_EFFECT",
    "LOSS_OF_FUNCTION",
    "DBSNP",
    "CLINVAR",
    "CLINVAR_REVIEW_STATUS_STARS",
    "CLINVAR_CONFLICTED",
    "CLINVAR_VARIANT_ORIGIN",
    "GWAS_CITATION",
    "ONCOGENE",
    "TUMOR_SUPPRESSOR",
    "CANCERGENE_EVIDENCE",
    "gnomADe_AF",
    "GENOMIC_CHANGE",
    "GENOME_VERSION"
  )

## define tags/variables to display in data tables (GWAS findings)
col_format_output[['html_gwas']] <-
  c(
    "SYMBOL",
    "CONSEQUENCE",
    "GWAS_CITATION",
    "PROTEIN_CHANGE",
    "GENOTYPE",
    "LOSS_OF_FUNCTION",
    "OFFICIAL_GENENAME",
    "GWAS_PHENOTYPE",
    "PROTEIN_DOMAIN",
    "HGVSp",
    "HGVSc",
    "REFSEQ",
    "CDS_CHANGE",
    "CODING_STATUS",
    "miRNA_TARGET_HIT",
    "miRNA_TARGET_HIT_PREDICTION",
    "TF_BINDING_SITE_VARIANT",
    "TF_BINDING_SITE_VARIANT_INFO",
    "GERP_SCORE",
    "PREDICTED_EFFECT",
    "DBSNP",
    "gnomADe_AF",
    "GENOMIC_CHANGE",
    "GENOME_VERSION"
  )

## define tags/variables to display in output TSV
col_format_output[['tsv']] <-
  c("SAMPLE_ID",
    "GENOMIC_CHANGE",
    "VAR_ID",
    "GENOTYPE",
    "GENOME_VERSION",
    "VARIANT_CLASS",
    "SYMBOL",
    "GENENAME",
    "CCDS",
    "ENTREZGENE",
    "UNIPROT_ID",
    "ENSEMBL_GENE_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "ONCOGENE",
    "TUMOR_SUPPRESSOR",
    "CONSEQUENCE",
    "PROTEIN_CHANGE",
    "PFAM_DOMAIN_NAME",
    "HGVSp",
    "HGVSc",
    "CDS_CHANGE",
    "LAST_EXON",
    "EXON",
    "EXON_AFFECTED",
    "CODING_STATUS",
    "EXON_POSITION",
    "INTRON_POSITION",
    "VEP_ALL_CSQ",
    "CANCER_PHENOTYPE",
    "MUTATION_HOTSPOT",
    "RMSK_HIT",
    "EFFECT_PREDICTIONS",
    "LOSS_OF_FUNCTION",
    "NULL_VARIANT",
    "DBMTS",
    "miRNA_TARGET_HIT",
    "miRNA_TARGET_HIT_PREDICTION",
    "REGULATORY_ANNOTATION",
    "TF_BINDING_SITE_VARIANT",
    "TF_BINDING_SITE_VARIANT_INFO",
    "GERP_SCORE",
    "DBSNP",
    "CLINVAR_CLASSIFICATION",
    "CLINVAR_MSID",
    "CLINVAR_VARIANT_ORIGIN",
    "CLINVAR_CONFLICTED",
    "CLINVAR_PHENOTYPE",
    "CLINVAR_REVIEW_STATUS_STARS",
    "N_INSILICO_CALLED",
    "N_INSILICO_DAMAGING",
    "N_INSILICO_TOLERATED",
    "N_INSILICO_SPLICING_NEUTRAL",
    "N_INSILICO_SPLICING_AFFECTED",
    "gnomADe_AF"
  )

## define tags/variables to display in output Excel
col_format_output[['xlsx_classification']] <-
  c("SAMPLE_ID",
    "GENOMIC_CHANGE",
    "GENOTYPE",
    "GENOME_VERSION",
    "VARIANT_CLASS",
    "SYMBOL",
    "GENENAME",
    "CONSEQUENCE",
    "PROTEIN_CHANGE",
    "FINAL_CLASSIFICATION",
    "CPSR_CLASSIFICATION_SOURCE",
    "CLINVAR_CLASSIFICATION",
    "CPSR_CLASSIFICATION",
    "CPSR_CLASSIFICATION_CODE",
    "CPSR_PATHOGENICITY_SCORE",
    "ENTREZGENE",
    "ENSEMBL_GENE_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "ONCOGENE",
    "TUMOR_SUPPRESSOR",
    "PFAM_DOMAIN_NAME",
    "HGVSp",
    "HGVSc",
    "CDS_CHANGE",
    "CODING_STATUS",
    "MUTATION_HOTSPOT",
    "EFFECT_PREDICTIONS",
    "LOSS_OF_FUNCTION",
    "NULL_VARIANT",
    "DBSNP",
    "CLINVAR_MSID",
    "CLINVAR_VARIANT_ORIGIN",
    "CLINVAR_CONFLICTED",
    "CLINVAR_PHENOTYPE",
    "CLINVAR_REVIEW_STATUS_STARS",
    "gnomADe_AF"
  )

## define tags/variables to display in output Excel
col_format_output[['xlsx_secondary']] <-
  c("SAMPLE_ID",
    "GENOMIC_CHANGE",
    "GENOTYPE",
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
    "CLINVAR_REVIEW_STATUS_STARS",
    "ENSEMBL_GENE_ID",
    "ENSEMBL_TRANSCRIPT_ID",
    "REFSEQ_TRANSCRIPT_ID",
    "PFAM_DOMAIN_NAME",
    "HGVSp",
    "HGVSc",
    "CDS_CHANGE",
    "CODING_STATUS",
    "MUTATION_HOTSPOT",
    "EFFECT_PREDICTIONS",
    "LOSS_OF_FUNCTION",
    "NULL_VARIANT",
    "DBSNP",
    "gnomADe_AF"
  )


## define tags/variables to display in output Excel
col_format_output[['xlsx_biomarker']] <-
  c("SAMPLE_ID",
    "GENOMIC_CHANGE",
    "GENOTYPE",
    "GENOME_VERSION",
    "VARIANT_CLASS",
    "SYMBOL",
    "GENENAME",
    "CONSEQUENCE",
    "PROTEIN_CHANGE",
    "CLINICAL_SIGNIFICANCE",
    "BIOMARKER_SOURCE",
    "BIOMARKER_RESOLUTION",
    "BIOMARKER_MATCH",
    "EVIDENCE_TYPE",
    "EVIDENCE_LEVEL",
    "EVIDENCE_DIRECTION",
    "EVIDENCE_ID",
    "THERAPEUTIC_CONTEXT",
    "EVIDENCE_DESCRIPTION",
    "MOLECULAR_PROFILE_NAME"
  )


#---- acmg ----#
acmg <- list()
acmg[["score2tier"]] <- data.frame()
acmg[["evidence_codes"]] <-
  utils::read.table(file = "data-raw/acmg_evidence.tsv",
                    header = T, stringsAsFactors = F,
                    comment.char = "", na.strings = c("NA"),
                    sep = "\t")
acmg[["pathogenic_range_gnomad"]] <- list()
acmg[["pathogenic_range_gnomad"]][["af"]] <- 0.0005
acmg[["pathogenic_range_gnomad"]][["min_an"]] <- 12000
acmg[["insilico_pred_min_majority"]] <- 8
acmg[["insilico_pred_max_minority"]] <- 2

acmg[["score2tier"]] <-
  data.frame("CPSR_CLASSIFICATION" = "Pathogenic",
             "CPSR_PATHOGENICITY_SCORE" = "<b>[5, ]</b>")
acmg[["score2tier"]] <-
  dplyr::bind_rows(
    acmg[["score2tier"]],
    data.frame("CPSR_CLASSIFICATION" = "Likely Pathogenic",
               "CPSR_PATHOGENICITY_SCORE" = "<b>[2.5, 4.5]</b>"))
acmg[["score2tier"]] <-
  dplyr::bind_rows(
    acmg[["score2tier"]],
    data.frame("CPSR_CLASSIFICATION" = "VUS",
               "CPSR_PATHOGENICITY_SCORE" = "<b>[-1.0, 2.0]</b>"))
acmg[["score2tier"]] <-
  dplyr::bind_rows(
    acmg[["score2tier"]],
    data.frame("CPSR_CLASSIFICATION" = "Likely Benign",
               "CPSR_PATHOGENICITY_SCORE" = "<b>[-4.5, -1.5]</b>"))
acmg[["score2tier"]] <-
  dplyr::bind_rows(
    acmg[["score2tier"]],
    data.frame("CPSR_CLASSIFICATION" = "Benign",
               "CPSR_PATHOGENICITY_SCORE" = "<b>[, -5]</b>"))


usethis::use_data(acmg, overwrite = T)
usethis::use_data(col_format_output, overwrite = T)
