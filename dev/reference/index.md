# Package index

## All functions

- [`acmg`](https://sigven.github.io/cpsr/dev/reference/acmg.md) : Scores
  and documentation of ACMG evidence criteria used for variant
  classification in CPSR
- [`append_cpg_properties()`](https://sigven.github.io/cpsr/dev/reference/append_cpg_properties.md)
  : Function that appends mechanism-of-inheritance (MOI) and mechanism
  of disease (MOD) annotations to cancer predisposition genes (cpg's),
  as well as estimated fractions of truncation vs. non-truncating
  variants etc per predisposition gene (from ClinVar)
- [`assign_classification()`](https://sigven.github.io/cpsr/dev/reference/assign_classification.md)
  : Function that assigns final pathogenicity classification (B, LB,
  VUS, P, LP) based on accumulated scores from different ACMG criteria
  and pre-defined cutoffs (calibrated against ClinVar)
- [`assign_pathogenicity_evidence()`](https://sigven.github.io/cpsr/dev/reference/assign_pathogenicity_evidence.md)
  : Function that assigns variant pathogenicity evidence based on ACMG
  guidelines
- [`bs_icon2()`](https://sigven.github.io/cpsr/dev/reference/bs_icon2.md)
  : Use Bootstrap icons (as inline SVG)
- [`check_variant2cancer_phenotype()`](https://sigven.github.io/cpsr/dev/reference/check_variant2cancer_phenotype.md)
  : Function that retrieves variants in cancer predisposition genes
  linked to cancer-related conditions according to ClinVar
- [`col_format_output`](https://sigven.github.io/cpsr/dev/reference/col_format_output.md)
  : Format of CPSR output data frames (HTML, TSV)
- [`color_palette`](https://sigven.github.io/cpsr/dev/reference/color_palette.md)
  : CPSR color palette
- [`combine_novel_and_preclassified()`](https://sigven.github.io/cpsr/dev/reference/combine_novel_and_preclassified.md)
  : Function that combines classifications of novel and pre-classified
  variants
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
- [`load_germline_snv_indel()`](https://sigven.github.io/cpsr/dev/reference/load_germline_snv_indel.md)
  : Function that reads and validates an annotated germline SNV/InDel
  file from CPSR pre-reporting pipeline
- [`plot_summary_statistics()`](https://sigven.github.io/cpsr/dev/reference/plot_summary_statistics.md)
  : Function that makes a piechart showing the number of variants at
  each significance level
- [`plot_virtual_panels()`](https://sigven.github.io/cpsr/dev/reference/plot_virtual_panels.md)
  : Function that makes a HTML display of virtual gene panel
- [`report_color`](https://sigven.github.io/cpsr/dev/reference/report_color.md)
  : CPSR report color
- [`retrieve_pgx_calls()`](https://sigven.github.io/cpsr/dev/reference/retrieve_pgx_calls.md)
  : Function that retrieves variants in genes recommended for secondary
  findings
- [`retrieve_secondary_calls()`](https://sigven.github.io/cpsr/dev/reference/retrieve_secondary_calls.md)
  : Function that retrieves variants in genes recommended for secondary
  findings
- [`write_cpsr_output()`](https://sigven.github.io/cpsr/dev/reference/write_cpsr_output.md)
  : Function that writes contents of CPSR report object to various
  output formats (quarto HTML reports, TSV, XLSX workbooks etc)
