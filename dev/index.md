  

## Cancer Predisposition Sequencing Reporter [![cpsr logo](reference/figures/logo.png)](https://sigven.github.io/cpsr/)

  
  

The *Cancer Predisposition Sequencing Reporter (CPSR)* is a
computational workflow that **interprets DNA sequence variants**
identified from next-generation sequencing **in the context of cancer
predisposition**.

*CPSR* accepts a query file with *germline* variant calls (SNVs/InDels)
from a single sample (i.e. cancer patient), encoded in the [VCF
format](https://samtools.github.io/hts-specs/VCFv4.2.pdf). Through
comprehensive gene and variant annotation procedures, CPSR offers the
following functionalities to the user:

1.  Flexible **selection of cancer predisposition genes** that restricts
    variant classification and reporting - through the use of virtual
    gene panels
2.  **Variant classification** (*Pathogenic* to *Benign*) through a
    dedicated and well-performing implementation of [ACMG/AMP
    guidelines](https://pubmed.ncbi.nlm.nih.gov/25741868/)
3.  **Detection of biomarkers** - variants with prognostic, diagnostic,
    or drug sensitivity/resistance implications in cancer, as well as
    optional detection of variants related to adverse events/toxicity
    for common chemotherapies  
4.  Optional reporting of **secondary/incidental findings** ([ACMG
    recommendations](https://pubmed.ncbi.nlm.nih.gov/37347242/))
5.  **Interactive HTML output report** with detailed variant
    information, gene annotations, and external links to relevant
    databases

The CPSR workflow is integrated with the framework that underlies
[Personal Cancer Genome Reporter -
PCGR](https://github.com/sigven/pcgr). While *PCGR* is intended for
reporting and analysis of somatic variants detected in a tumor, *CPSR*
is intended for reporting and ranking of germline variants in
protein-coding genes that are implicated in cancer predisposition and
inherited cancer syndromes.

Four snapshots of sections in the [quarto](https://quarto.org)-based
cancer predisposition genome report (artificial sample, with more
findings than usual):

  

![Virtual gene panel & summary of findings](img/cpsr_panel_summary.png)

![Variant classification](img/cpsr_variant_classification.png)

![Genomic biomarkers](img/cpsr_biomarkers.png)

![Pharmacogenetic & secondary findings](img/cpsr_pgx_secondary.png)

  

### News

- *September 17th 2025*: **2.2.5 release**
  - patch - safeguard against missing data in gnomAD non-cancer variant
    data
- *September 8th 2025*: **2.2.4 release**
  - patch to avoid duplicate matching of PVS1 criteria
- *March 23rd 2025*: **2.2.1 release**
  - patch to fix bug with non-standard ClinVar significance levels (Drug
    Response, Risk Factor)
- *March 22nd 2025*: **2.2.0 release**
  - more predisposition genes in panel zero
  - optional reporting of pharmacogenomics-related variants (TPMT, DPYD,
    NUDT15)
  - [CHANGELOG](https://sigven.github.io/cpsr/articles/CHANGELOG.html)
- *October 2024*: **2.1.2 release**
  - cosmetic fixes in HTML report
  - fix for VEP consequence pick exception
  - [CHANGELOG](https://sigven.github.io/cpsr/articles/CHANGELOG.html)
- *September 2024*: **2.1.0 release**
  - data bundle upgrade
  - re-calibration of classification tresholds
  - [CHANGELOG](https://sigven.github.io/cpsr/articles/CHANGELOG.html)
- *June 2024*: **2.0.0 release**
  - New HTML report generation and layout with
    [quarto](https://quarto.org/)
  - Excel output supported
  - Data bundle update
  - Singularity/Apptainer support
  - [CHANGELOG](http://cpsr.readthedocs.io/en/latest/CHANGELOG.md)
- *November 2022*: **1.0.1 release**
  - Added CPSR logo (designed by [Hal Nakken](https://halvetica.net))

### Example report

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17117326.svg)](https://doi.org/10.5281/zenodo.17117326)

### Getting started

- [Installation
  instructions](https://sigven.github.io/cpsr/dev/articles/installation.md)

- [Run through an
  example](https://sigven.github.io/cpsr/dev/articles/running.html#example-run)

- Learn more about

  1.  Details regarding [CPSR input
      files](https://sigven.github.io/cpsr/dev/articles/input.md), and
      how they should be formatted
  2.  The types and contents of [CPSR output
      files](https://sigven.github.io/cpsr/dev/articles/output.md)
  3.  [ACMG variant classification
      procedure](https://sigven.github.io/cpsr/dev/articles/variant_classification.md)
      used in CPSR
  4.  The list of [virtual gene
      panels](https://sigven.github.io/cpsr/dev/articles/virtual_panels.md)
      available in CPSR

### Contact

<sigven@ifi.uio.no>
