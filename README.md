# Cancer Predisposition Sequencing Reporter <a href="https://sigven.github.io/cpsr/"><img src="man/figures/logo.png" align="right" height="118" width="100"/></a>

The *Cancer Predisposition Sequencing Reporter (CPSR)* is a computational workflow that **interprets and classifies the clinical significance of germline DNA variants** identified from next-generation sequencing **in the context of cancer predisposition and inherited cancer syndromes**. The workflow can also report **incidental findings (ACMG v3.2)**.

The CPSR workflow is integrated with the framework that underlies the [Personal Cancer Genome Reporter - PCGR](https://github.com/sigven/pcgr). While *PCGR* is intended for reporting and analysis of somatic variants detected in a tumor, *CPSR* is intended for reporting and ranking of germline variants in protein-coding genes that are implicated in cancer predisposition and inherited cancer syndromes. 

*CPSR* accepts a query file from a single case/patient, containing raw germline variant calls encoded in the [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format (i.e. SNVs/InDels). A comprehensive set of **virtual cancer predisposition gene panels** harvested from the [Genomics England PanelApp](https://panelapp.genomicsengland.co.uk/) allows the user to flexibly put a restriction on which genes and findings are displayed in the cancer predisposition report.

Snapshots of sections in the [quarto](https://quarto.org)-based cancer predisposition genome report (artificial sample, with more findings than usual) are shown below:

![CPSR views](pkgdown/assets/img/cpsr_sc.png)

## News

-   *July 2024*: **2.0.1 release**
    - patch with bug fix for mitochondrial input variants ([pr245](https://github.com/sigven/pcgr/pull/245))
    - [CHANGELOG](https://sigven.github.io/cpsr/articles/CHANGELOG.html)

-   *June 2024*: **2.0.0 release**
    - New HTML report generation and layout with [quarto](https://quarto.org/)
    - Excel output supported
    - Data bundle update
    - Singularity/Apptainer support
    - [CHANGELOG](https://sigven.github.io/cpsr/articles/CHANGELOG.html)

-   *November 2022*: **1.0.1 release**
    -   Added CPSR logo (designed by [Hal Nakken](https://halvetica.net))
    
-   *February 2022*: **1.0.0 release**
    -   Complete restructure of code and Conda installation routines, contributed largely by the great [@pdiakumis](https://github.com/pdiakumis)
    -   Updated data bundle
        - ClinVar - Feb 2022
        - CancerMine - Dec 2021
        - UniprotKB - Nov 2021
        - CIViC - Feb 2022
        - GWAS catalog - Dec 2021
    -   Software upgrade (VEP 105, R/BioConductor)
    -   New documentation site ([https://sigven.github.io/cpsr](https://sigven.github.io/cpsr))


## Example report

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11401491.svg)](https://doi.org/10.5281/zenodo.11401491)

## Getting started

-   [Installation instructions](https://sigven.github.io/cpsr/articles/installation.html)
-   [Run through an example](https://sigven.github.io/cpsr/articles/running.html#example-run)
-   Learn more about
    * Details regarding the [CPSR input file](https://sigven.github.io/cpsr/articles/input.html), and how it should be formatted
    * The types and contents of [CPSR output files](https://sigven.github.io/cpsr/articles/output.html)
    * [ACMG variant classification procedure](https://sigven.github.io/cpsr/articles/variant_classification.html) used in CPSR
    * The list of [virtual gene panels](https://sigven.github.io/cpsr/articles/virtual_panels.html) available in CPSR

## Citation

If you use CPSR, please cite the following publication:

Sigve Nakken, Vladislav Saveliev, Oliver Hofmann, Pål Møller, Ola Myklebost, and Eivind Hovig. **Cancer Predisposition Sequencing Reporter (CPSR): a flexible variant report engine for high-throughput germline screening in cancer** (2021). *Int J Cancer*. [doi:[10.1002/ijc.33749](doi:%5B10.1002/ijc.33749)](https://doi.org/10.1002/ijc.33749)

## Contact

[sigven\@ifi.uio.no](mailto:sigven@ifi.uio.no)
