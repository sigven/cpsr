<br>

## Cancer Predisposition Sequencing Reporter <a href="https://sigven.github.io/cpsr/"><img src="man/figures/logo.png" align="right" height="106" width="90"/></a>

<br><br>

The *Cancer Predisposition Sequencing Reporter (CPSR)* is a computational workflow that **interprets germline variants** identified from next-generation sequencing **in the context of cancer predisposition**. 

*CPSR* accepts a query file with raw germline variant calls (SNVs/InDels) from a single sample (cancer patient), encoded in the [VCF format ](https://samtools.github.io/hts-specs/VCFv4.2.pdf). CPSR conducts comprehensive gene and variant annotation on the input calls, and generates a dedicated _variant HTML report_, that provides the following main functionality:

1) Flexible **selection of cancer predisposition genes** subject to analysis and reporting
2) **Variant classification** (*Pathogenic* to _Benign_) through implementation of ACMG guidelines
3) **Biomarker matching** of sample variants (prognosis, diagnosis, drug sensitivity/resistance)
4) Reporting of **secondary/incidental findings** (ACMG recommendations)


The workflow is integrated with the framework that underlies [Personal Cancer Genome Reporter - PCGR ](https://github.com/sigven/pcgr). While *PCGR* is intended for reporting and analysis of somatic variants detected in a tumor, *CPSR* is intended for reporting and ranking of germline variants in protein-coding genes that are implicated in cancer predisposition and inherited cancer syndromes.

Snapshots of sections in the [quarto](https://quarto.org)-based cancer predisposition genome report ((artificial sample, with more findings than usual):

![](img/cpsr_sc.png)

<br>

### News

* *May 2024*: **2.x.x release**
  - New HTML report generation and layout with [quarto](https://quarto.org/)
  - Excel output supported
  - Updated virtual gene panels (Genomics England PanelApp, Cancer Gene Census)
  - Reference data updates, most importantly including 
    - ClinVar - June 2024
    - CIViC - June 2024
    - GENCODE - v46
  - Software updates - VEP 112
  - Extensive code clean-up and re-structuring

* *November 2022*: **1.0.1 release**
  * Added CPSR logo (designed by [Hal Nakken](https://halvetica.net))

* *February 2022*: **1.0.0 release**
  * Complete restructure of code and Conda installation routines, contributed largely by the great [@pdiakumis](https://github.com/pdiakumis)
  * Updated bundle (ClinVar, CancerMine, UniprotKB, PanelApp, CIViC, GWAS catalog)
  * Software upgrade (VEP, R/BioConductor)
  * New documentation site (https://sigven.github.io/cpsr)

* *June 30th 2021*: **0.6.2 release**
  * Updated bundle (ClinVar, CancerMine, UniprotKB, PanelApp, CIViC, GWAS catalog)
  * Software upgrade (VEP, R/BioConductor)
  * [CHANGELOG](http://cpsr.readthedocs.io/en/latest/CHANGELOG.html)
* *November 30th 2020*: **0.6.1 release**
  * Updated bundle (ClinVar, CancerMine, UniprotKB, CIViC, GWAS catalog)
  * [CHANGELOG](http://cpsr.readthedocs.io/en/latest/CHANGELOG.html)


### Example report

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11401491.svg)](https://doi.org/10.5281/zenodo.11401491)

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
