Variants in the cancer predisposition gene set are classified through two complementary
assertion authorities:

- <b style="color:#e65100;">CPSR</b> — All coding variants are classified *de novo* by
  CPSR according to a <i>five-level pathogenicity scheme</i> (*CPSR_CLASSIFICATION*):
  pathogenic / likely pathogenic / VUS / likely benign / benign. The classification is
  rule-based, implementing a key subset of published **ACMG/AMP criteria**. Gene-level
  attributes relevant to pathogenicity assessment — including mode of inheritance and
  mechanism of disease (loss-of-function vs. gain-of-function) — are sourced from
  [Genomics England PanelApp](https://panelapp.genomicsengland.co.uk/),
  [Maxwell et al., Am J Hum Genet, 2016](https://www.ncbi.nlm.nih.gov/pubmed/27153395),
  and [Huang et al., Cell, 2018](https://www.ncbi.nlm.nih.gov/pubmed/29625052).
  Importantly, note that some ACMG/AMP criteria have been accommodated with
  gene-specific recommendations:
     * BA1/BS1/PM2 allele-specific thresholds (as specified by
       [ClinGen VCEPs](https://cspec.genome.network/cspec/ui/svi/) for key
       cancer-predisposing genes, i.e. BRCA1/2, MMR genes, PALB2, APC, ATM, PTEN, TP53).
     * PM4 - exceptions for TP53 and ATM
     * PM1 - regional considerations for matching against TP53 and PTEN, and exceptions
       for PALB2, ATM, APC, MMR genes and BRCA1/2

- <b style="color:#0277bd;">ClinVar</b> — For variants with an existing record in
  [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), the ClinVar interpretation may
  override the CPSR-computed classification as the final reported verdict. Whether this
  override applies depends on the *clinvar_trust_level* setting: at higher trust levels,
  only classifications backed by expert panel review or multiple submitters with no
  conflicts are accepted; at lower trust levels, a broader set of ClinVar submissions is
  considered sufficient to override CPSR. If a variant's ClinVar record does not meet
  the configured trust threshold, the CPSR rule-based classification is retained as the
  final call.

The ACMG/AMP criteria listed in {criteria_ref} form the basis for the
*CPSR_CLASSIFICATION* variable. The <i>score</i> column indicates how much each evidence
item contributes to either of the two pathogenicity poles (positive values indicate
pathogenic support, negative values indicate benign support). Scores along each pole
('B' and 'P') are aggregated and represented through the *CPSR_PATHOGENICITY_SCORE*
variable. See {threshold_ref} on how CPSR establishes optimal thresholds for converting
pathogenicity scores to categorical classification levels.
