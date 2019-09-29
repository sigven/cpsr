
## CHANGELOG

#### 0.5.0 - September 23rd 2019

##### Fixed
  * Bug in implementation of ACMG criteria; genes without a known loss-of-function mechanism were handled inappropriately
  * Bug in assignment of heterozygous/homozygous states (input VCF)
  * Bug in implementation of ACMG_PS1 - Same amino acid change as previously pathogenic variant
  * Improved consequence prioritisation for variants with transcript consequences in multiple, distinct cancer predisposition genes
  * Upper MAF threshold (as given by user) only applied for unclassified (i.e. non-ClinVar variants)
  * Handling of non-coding variants (synonymous, upstream_variants) in the report, no longer excluded

##### Added
  * Section on _genomic biomarkers_; indicating which variants in the query VCF that overlaps with existing germline biomarkers (CIViC)
