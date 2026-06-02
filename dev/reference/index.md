# Package index

## All functions

- [`acmg`](https://sigven.github.io/cpsr/dev/reference/acmg.md) : Scores
  and documentation of ACMG evidence criteria used for variant
  classification in CPSR
- [`append_cpg_properties()`](https://sigven.github.io/cpsr/dev/reference/append_cpg_properties.md)
  : Function that appends mechanism-of-inheritance (MOI) and mechanism
  of disease (MOD) annotations wrt cancer predisposition genes (cpg's),
  as well as estimated fractions of truncation vs. non-truncating
  variants etc per predisposition gene (from ClinVar), and pathogenic AF
  ranges based on gnomAD data, to the cpg_calls data frame.
- [`assign_BA1_evidence()`](https://sigven.github.io/cpsr/dev/reference/assign_BA1_evidence.md)
  : Function that assigns ACMG evidence indicators for BA1
- [`assign_BP1_evidence()`](https://sigven.github.io/cpsr/dev/reference/assign_BP1_evidence.md)
  : Function that assigns ACMG evidence indicators for BP1
- [`assign_BP4_evidence()`](https://sigven.github.io/cpsr/dev/reference/assign_BP4_evidence.md)
  : Function that assigns ACMG evidence indicators for BP4
- [`assign_BP7_evidence()`](https://sigven.github.io/cpsr/dev/reference/assign_BP7_evidence.md)
  : Function that assigns ACMG evidence indicators for BP7
- [`assign_BS1_evidence()`](https://sigven.github.io/cpsr/dev/reference/assign_BS1_evidence.md)
  : Function that assigns ACMG evidence indicators for BS1_ST and
  BS1_SUP
- [`assign_PM1_evidence()`](https://sigven.github.io/cpsr/dev/reference/assign_PM1_evidence.md)
  : Function that assigns ACMG evidence indicators for PM1 ('ACMG_PM1'
  and 'ACMG_PM1_SUPP')
- [`assign_PM2_evidence()`](https://sigven.github.io/cpsr/dev/reference/assign_PM2_evidence.md)
  : Function that assigns ACMG evidence indicators for PM2_supporting
- [`assign_PM4_evidence()`](https://sigven.github.io/cpsr/dev/reference/assign_PM4_evidence.md)
  : Function that assigns ACMG evidence indicators for PM4
- [`assign_PM5_evidence()`](https://sigven.github.io/cpsr/dev/reference/assign_PM5_evidence.md)
  : Function that assigns ACMG evidence indicators for PM5
- [`assign_PP2_evidence()`](https://sigven.github.io/cpsr/dev/reference/assign_PP2_evidence.md)
  : Function that assigns ACMG evidence indicators for PP2
- [`assign_PP3_evidence()`](https://sigven.github.io/cpsr/dev/reference/assign_PP3_evidence.md)
  : Function that assigns ACMG evidence indicators for PP3
- [`assign_PS1_evidence()`](https://sigven.github.io/cpsr/dev/reference/assign_PS1_evidence.md)
  : Function that assigns ACMG evidence indicators for PS1
  ('ACMG_PS1') - Same amino acid change as a previously established
  pathogenic variant regardless of nucleotide change
- [`assign_PVS1_evidence()`](https://sigven.github.io/cpsr/dev/reference/assign_PVS1_evidence.md)
  : Function that assigns ACMG evidence indicators for PVS1
  ('ACMG_PVS1', 'ACMG_PVS1_STR', 'ACMG_PVS1_MOD')
- [`assign_acmg_consensus()`](https://sigven.github.io/cpsr/dev/reference/assign_acmg_consensus.md)
  : Function that assigns a summary string of all ACMG evidence codes
  met for a given variant (e.g. "PVS1\|PS1\|PM2_supporting\|PP3")
- [`assign_acmg_evidence()`](https://sigven.github.io/cpsr/dev/reference/assign_acmg_evidence.md)
  : Function that assigns variant pathogenicity evidence based on ACMG
  guidelines
- [`assign_classification_authority()`](https://sigven.github.io/cpsr/dev/reference/assign_classification_authority.md)
  : Function that combines classifications of novel and pre-classified
  variants
- [`bs_icon2()`](https://sigven.github.io/cpsr/dev/reference/bs_icon2.md)
  : Use Bootstrap icons (as inline SVG)
- [`classification_thresholds_intro_doc_note()`](https://sigven.github.io/cpsr/dev/reference/classification_thresholds_intro_doc_note.md)
  : Get documentation string for classification threshold calibration
  (intro)
- [`classification_thresholds_outro_doc_note()`](https://sigven.github.io/cpsr/dev/reference/classification_thresholds_outro_doc_note.md)
  : Get documentation string for classification threshold calibration
  (outro)
- [`col_format_output`](https://sigven.github.io/cpsr/dev/reference/col_format_output.md)
  : Format of CPSR output data frames (HTML, TSV)
- [`color_palette`](https://sigven.github.io/cpsr/dev/reference/color_palette.md)
  : CPSR color palette
- [`cpsr_intro_doc_note()`](https://sigven.github.io/cpsr/dev/reference/cpsr_intro_doc_note.md)
  : Get documentation string for the CPSR report introduction
- [`create_clinvar_reactable()`](https://sigven.github.io/cpsr/dev/reference/create_clinvar_reactable.md)
  : Create a ClinVar variant reactable with expandable detail rows
- [`create_unified_variant_filters()`](https://sigven.github.io/cpsr/dev/reference/create_unified_variant_filters.md)
  : Create crosstalk-linked unified filters for all variants
- [`create_unified_variant_reactable()`](https://sigven.github.io/cpsr/dev/reference/create_unified_variant_reactable.md)
  : Create the complete unified variant reactable
- [`create_variant_table_theme()`](https://sigven.github.io/cpsr/dev/reference/create_variant_table_theme.md)
  : Applies consistent styling using CPSR color palette.
- [`curated_transcripts`](https://sigven.github.io/cpsr/dev/reference/curated_transcripts.md)
  : Curated transcripts
- [`exclude_vars_by_maf()`](https://sigven.github.io/cpsr/dev/reference/exclude_vars_by_maf.md)
  : Function that filters novel (non-ClinVar) variants based on a MAF
  threshold in gnomAD
- [`exclude_vars_non_cancer()`](https://sigven.github.io/cpsr/dev/reference/exclude_vars_non_cancer.md)
  : Function that ignores/exclude variants attributed to non-cancer
  related phenotypes in ClinVar
- [`exclude_vars_noncoding()`](https://sigven.github.io/cpsr/dev/reference/exclude_vars_noncoding.md)
  : Function that filters novel (non-ClinVar) variants based on a MAF
  threshold in gnomAD
- [`generate_cpsr_report()`](https://sigven.github.io/cpsr/dev/reference/generate_cpsr_report.md)
  : Function that generates variant predisposition report - CPSR
- [`get_insilico_prediction_statistics()`](https://sigven.github.io/cpsr/dev/reference/get_insilico_prediction_statistics.md)
  : Function that counts insilico predictions of variant effects (i.e.
  damaging/tolerated) from dbNSFP
- [`get_max_rows_pr_datatable()`](https://sigven.github.io/cpsr/dev/reference/get_max_rows_pr_datatable.md)
  : Function that gets the maximum number of rows across different tier
  data frames in CPSR report
- [`gnomad_germline_doc_note()`](https://sigven.github.io/cpsr/dev/reference/gnomad_germline_doc_note.md)
  : Get documentation string for gnomAD usage in germline variant
  classification
- [`insilico_path_predictors`](https://sigven.github.io/cpsr/dev/reference/insilico_path_predictors.md)
  : Insilico protein impact prediction methods
- [`is_clinvar_cancer_phenotype()`](https://sigven.github.io/cpsr/dev/reference/is_clinvar_cancer_phenotype.md)
  : Function that retrieves variants in cancer predisposition genes
  linked to cancer-related conditions according to ClinVar
- [`is_purely_intronic_deletion()`](https://sigven.github.io/cpsr/dev/reference/is_purely_intronic_deletion.md)
  : Function that checks if a variant is a purely intronic deletion
  (i.e. deletion of intronic sequence that does not extend into exonic
  sequence at either end)
- [`known_path_benign_sites()`](https://sigven.github.io/cpsr/dev/reference/known_path_benign_sites.md)
  : Function that extracts known benign peptide changes from ClinVar
  (\>= 2 gold stars), as well as known pathogenic nucleotide and codon
  sites
- [`load_germline_snv_indel()`](https://sigven.github.io/cpsr/dev/reference/load_germline_snv_indel.md)
  : Function that reads and validates an annotated germline SNV/InDel
  file from CPSR pre-reporting pipeline
- [`min_distance_to_pathogenic()`](https://sigven.github.io/cpsr/dev/reference/min_distance_to_pathogenic.md)
  : Function that calculates the closest distance to a known pathogenic
  variant in the same gene
- [`plot_summary_statistics()`](https://sigven.github.io/cpsr/dev/reference/plot_summary_statistics.md)
  : Function that makes a piechart showing the number of variants at
  each significance level
- [`plot_virtual_panels()`](https://sigven.github.io/cpsr/dev/reference/plot_virtual_panels.md)
  : Function that makes a HTML display of virtual gene panel
- [`prep_biomarker_tbl()`](https://sigven.github.io/cpsr/dev/reference/prep_biomarker_tbl.md)
  : Function that gathers data table on biomarker variants for display
  in germline report
- [`prepare_unified_all_variants()`](https://sigven.github.io/cpsr/dev/reference/prepare_unified_all_variants.md)
  : Prepare unified variant dataset for single-table display
- [`render_actble_bm_table()`](https://sigven.github.io/cpsr/dev/reference/render_actble_bm_table.md)
  : Build biomarker reactable with category-aware styling
- [`report_color`](https://sigven.github.io/cpsr/dev/reference/report_color.md)
  : CPSR report color
- [`retrieve_pgx_calls()`](https://sigven.github.io/cpsr/dev/reference/retrieve_pgx_calls.md)
  : Function that retrieves variants in genes recommended for secondary
  findings
- [`retrieve_secondary_calls()`](https://sigven.github.io/cpsr/dev/reference/retrieve_secondary_calls.md)
  : Function that retrieves variants in genes recommended for secondary
  findings
- [`rt_cell_bm_significance()`](https://sigven.github.io/cpsr/dev/reference/rt_cell_bm_significance.md)
  : Cell renderer factory for biomarker clinical significance Returns a
  function that renders biomarker clinical significance values as
  colored pills based on the provided color palette.
- [`rt_cell_classification()`](https://sigven.github.io/cpsr/dev/reference/rt_cell_classification.md)
  : Cell renderer factory for classification column Returns a function
  that renders classification values as colored pills based on the
  provided color palette.
- [`rt_cell_classification_rank()`](https://sigven.github.io/cpsr/dev/reference/rt_cell_classification_rank.md)
  : Cell renderer factory for classification rank column Maps a numeric
  rank (1=Benign ... 5=Pathogenic) to a labelled colored pill, enabling
  correct sort order via JavaScript.
- [`rt_cell_genotype()`](https://sigven.github.io/cpsr/dev/reference/rt_cell_genotype.md)
  : Cell renderer factory for genotype column Returns a function that
  renders genotype values as colored pills based on the provided color
  palette.
- [`variant_classification_synopsis_doc_note()`](https://sigven.github.io/cpsr/dev/reference/variant_classification_synopsis_doc_note.md)
  : Get documentation string for variant classification synopsis
- [`write_cpsr_output()`](https://sigven.github.io/cpsr/dev/reference/write_cpsr_output.md)
  : Function that writes contents of CPSR report object to various
  output formats (quarto HTML reports, TSV, XLSX workbooks etc)
