
## CHANGELOG

#### 0.6.0rc - September 24th 2020

- Data updates: ClinVar, GWAS catalog, GENCODE, CIViC, CancerMine, UniProt KB, dbNSFP, Pfam, KEGG, Open Targets Platform, Genomics England PanelApp
- Software updates: VEP 101

##### Fixed
  * Duplicated entries in incidental findings

##### Changed
  * All arguments to `cpsr.py` are now non-positional
  * Arguments to `cpsr.py` are divided into two groups: _required_ and _optional_
  * `secondary_findings` is now coined `incidental_findings`
  * Option ___gwas:gwas_hits___ in CPSR configuration file is now optional argument `--gwas_findings` in `cpsr.py`
  * Option ___classification:clinvar_cpsr___ in CPSR configuration file is now optional argument `--classify_all` in `cpsr.py`
  * Option ___maf_imits:maf_gnomad___ in CPSR configuration file is now optional argument `--maf_upper_threshold` in `cpsr.py`
  * Option ___secondary_findings:show_sf___ in CPSR configuration file is now optional argument `--incidental_findings` in `cpsr.py`
  * Virtual panels is now displayed through HTML (previously static ggplot plot)
  * __Settings__ section of report is now divived into three:
	  * Sample metadata
	  * Report configuration
	  * Virtual panel
  * Classifications of genes as tumor suppressors/oncogenes are now based on a combination of CancerMine citation count and presence in Network of Cancer Genes

##### Added
  * Missing ACMG criterion for classification of silent and intronic variants outside of splice regions (_ACMG_BP7_)
  * Missing ACMG criterion for classification of variants in promoter and untranslated regions (_ACMG_BP3_)
  * Possibility to create custom virtual panel - any combination of genes from panel 0 provided as a single-column text file with argument `--custom_list`
  * Ensured that non-empty datatables pr. tier (__ClinVar__ and __Non-ClinVar__) are set as the active tab
  * Improved documentation of variant classification in the __References__ section
  * DOIs available for all references

#### 0.5.2 - November 18th 2019

##### Changed
  * Definition of pathogenic range (wrt to variant frequency) takes into account population- and position-specific allele numbers - no longer defined only by allele counts (i.e. *AC*) but by *AC* and *AN*
  * Moved virtual panel identifier from positional argument to optional argument (`--panel_id`) in `cpsr.py`

##### Added
  * Ability to analyze custom panels, provided through option `--custom_panel`. Needs to be defined as BED file with four columns, i.e. chromosome, start, stop, genesymbol

#### 0.5.1 - October 14th 2019

##### Fixed
  * Bug in `cpsr_validate_input.py`, [GitHub Issue](https://github.com/sigven/cpsr/issues/18)
  * Bug when there are zero variants with a 'PASS' status in VCF - omitting report generation

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
