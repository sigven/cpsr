<br>

## Cancer Predisposition Sequencing Reporter <a href="https://sigven.github.io/cpsr/"><img src="man/figures/logo.png" align="right" height="106" width="90"/></a>

<br><br>

The *Cancer Predisposition Sequencing Reporter (CPSR)* is a computational workflow that **interprets DNA sequence variants** identified from next-generation sequencing **in the context of cancer predisposition**. 

*CPSR* accepts a query file with _germline_ variant calls (SNVs/InDels) from a single sample (i.e. cancer patient), encoded in the [VCF format ](https://samtools.github.io/hts-specs/VCFv4.2.pdf). Through comprehensive gene and variant annotation procedures, CPSR offers the following main functionalities to the user:

1) Flexible **selection of cancer predisposition genes** subject to analysis and reporting - through the use of virtual gene panels
2) **Variant classification** (*Pathogenic* to _Benign_) through a dedicated implementation of [ACMG/AMP guidelines](https://pubmed.ncbi.nlm.nih.gov/25741868/)
3) **Detection of germline biomarkers** - for prognosis, diagnosis, or drug sensitivity/resistance in cancer
4) Reporting of **secondary/incidental findings** ([ACMG recommendations](https://pubmed.ncbi.nlm.nih.gov/37347242/))
5) **Interactive HTML output report** with detailed variant information, gene annotations, and external links to relevant databases


The CPSR workflow is integrated with the framework that underlies [Personal Cancer Genome Reporter - PCGR ](https://github.com/sigven/pcgr). While *PCGR* is intended for reporting and analysis of somatic variants detected in a tumor, *CPSR* is intended for reporting and ranking of germline variants in protein-coding genes that are implicated in cancer predisposition and inherited cancer syndromes.

Snapshots of sections in the [quarto](https://quarto.org)-based cancer predisposition genome report (artificial sample, with more findings than usual):

![](img/cpsr_sc.png)

<br>

### News

*  *October 2024*: **2.1.2 release**
    -  cosmetic fixes in HTML report
    -  fix for VEP consequence pick exception
    -  [CHANGELOG](https://sigven.github.io/cpsr/articles/CHANGELOG.html)
    
* *September 2024*: **2.1.0 release**
  - data bundle upgrade
  - re-calibration of classification tresholds
  - [CHANGELOG](https://sigven.github.io/cpsr/articles/CHANGELOG.html)
    
* *August 2024*: **2.0.3 release**
  - patch with bug fix for copying in quarto templates ([pr62](https://github.com/sigven/cpsr/pull/62))
  - [CHANGELOG](https://sigven.github.io/cpsr/articles/CHANGELOG.html)
    
* *July 2024*: **2.0.1 release**
  - patch with bug fix for mitochondrial input variants ([pr245](https://github.com/sigven/pcgr/pull/245))
  - [CHANGELOG](https://sigven.github.io/cpsr/articles/CHANGELOG.html)

* *June 2024*: **2.0.0 release**
  - New HTML report generation and layout with [quarto](https://quarto.org/)
  - Excel output supported
  - Data bundle update
  - Singularity/Apptainer support
  - [CHANGELOG](http://cpsr.readthedocs.io/en/latest/CHANGELOG.html)

* *November 2022*: **1.0.1 release**
  * Added CPSR logo (designed by [Hal Nakken](https://halvetica.net))


### Example report

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13855978.svg)](https://doi.org/10.5281/zenodo.13855978)

### Getting started

- [Installation instructions](articles/installation.html)
- [Run through an example](articles/running.html#example-run)
- Learn more about

   1) Details regarding [CPSR input files](articles/input.html), and how they should be formatted 
   2) The types and contents of [CPSR output files](articles/output.html)
   3) [ACMG variant classification procedure](articles/variant_classification.html) used in CPSR
   4) The list of [virtual gene panels](articles/virtual_panels.html) available in CPSR


### Contact

sigven@ifi.uio.no
