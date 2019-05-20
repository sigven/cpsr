## About

###  What is the Cancer Predisposition Sequencing Preporter (CPSR)?

The *Cancer Predisposition Sequencing Reporter (CPSR)* is a computational workflow that **interprets germline variants** identified from next-generation sequencing **in the context of cancer predisposition**. The workflow is integrated with the framework that underlies the [Personal Cancer Genome Reporter (PCGR)](https://github.com/sigven), utilizing the Docker environment for encapsulation of code and software dependencies. While *PCGR* is intended for reporting and analysis of somatic variants detected in a tumor, *CPSR* is intended for reporting and ranking of germline variants in protein-coding genes that are implicated in cancer predisposition and inherited cancer syndromes.

*CPSR* accepts a query file with raw germline variant calls encoded in the [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format (i.e. analyzing SNVs/InDels). Furthermore, through the use several different _virtual cancer predisposition gene panels_ harvested from the [Genomics England PanelApp](https://panelapp.genomicsengland.co.uk/), the user can flexibly put a restriction on which genes and findings are displayed in the cancer predisposition report.

The software performs extensive variant annotation on the selected geneset and produces an interactive HTML report, in which the user can investigate four types of variants:

1. __ClinVar variants__ - pre-classified variants according to a five-level tier scheme (Pathogenic to Benign)
2. __Non-ClinVar variants__ - classified by CPSR according to a five-level tier scheme (Pathogenic to Benign)
3. __Secondary findings (optional)__ - pathogenic ClinVar variants in the ACMG recommended list for reporting of incidental findings
4. __GWAS hits (optional)__ - variants overlapping with previously identified hits in genome-wide association studies (GWAS) of cancer phenotypes (i.e. low to moderate risk conferring alleles), using [NHGRI-EBI Catalog of published genome-wide association studies](https://www.ebi.ac.uk/gwas/) as the underlying source.

The variant sets can be interactively explored and filtered further through different types of filters (phenotypes, genes, variant consequences, population MAF etc.). Importantly, the unclassified non-ClinVar variants are assigned a *pathogenicity level* based on the aggregation of scores according to previously established [ACMG criteria](https://www.ncbi.nlm.nih.gov/pubmed/25741868). The ACMG criteria includes cancer-specific criteria, as outlined and specified in several previous studies ([Huang et al., *Cell*, 2018](https://www.ncbi.nlm.nih.gov/pubmed/29625052); [Nykamp et al., *Genet Med.*, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28492532); [Maxwell et al., *Am J Hum Genet.*, 2016](https://www.ncbi.nlm.nih.gov/pubmed/27153395); [Amendola et al., *Am J Hum Genet.*,  2016](https://www.ncbi.nlm.nih.gov/pubmed/27181684)). See also [*Related work*](https://github.com/sigven/cpsr#related-work) below).

##### Cancer predisposition genes

The cancer predisposition report can show variants found in a number of well-known cancer predisposition genes, and the specific set of genes can be customized by the user by choosing any of the following __virtual gene panels (0 - 33)__:

  * **Panel 0 (default)** is a comprehensive gene panel assembled through known sources on cancer predisposition:
	* A list of 152 genes that were curated and established within TCGA’s pan-cancer study ([Huang et al., *Cell*, 2018](https://www.ncbi.nlm.nih.gov/pubmed/29625052))
	* A list of 107 protein-coding genes that has been manually curated in COSMIC’s [Cancer Gene Census v88](https://cancer.sanger.ac.uk/census),
	* A list of 148 protein-coding genes established by experts within the Norwegian Cancer Genomics Consortium (http://cancergenomics.no)

	The combination of the three sources resulted in a non-redundant set of [209 protein-coding genes](https://github.com/sigven/cpsr/blob/master/predisposition.md) of relevance for predisposition to tumor development.

  * **Panels 1 - 33** are collected from cancer predisposition panels assembled within the [Genomics England PanelApp](https://panelapp.genomicsengland.co.uk/):
	  * [1 = Adult solid tumours cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/245/)
	  * [2 = Adult solid tumours for rare disease](https://panelapp.genomicsengland.co.uk/panels/391/)
	  * [3 = Bladder cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/208/)
	  * [4 = Brain cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/166/)
	  * [5 = Breast cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/55/)
	  * [6 = Childhood solid tumours cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/259/)
	  * [7 = Colorectal cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/244/)
	  * [8 = Endometrial cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/271/)
	  * [9 = Familial Tumours Syndromes of the central & peripheral Nervous system](https://panelapp.genomicsengland.co.uk/panels/167/)
	  * [10 = Familial breast cancer](https://panelapp.genomicsengland.co.uk/panels/158/)
	  * [11 = Familial melanoma](https://panelapp.genomicsengland.co.uk/panels/522/)
	  * [12 = Familial prostate cancer](https://panelapp.genomicsengland.co.uk/panels/318/)
	  * [13 = Familial rhabdomyosarcoma](https://panelapp.genomicsengland.co.uk/panels/290/)
	  * [14 = Haematological malignancies cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/59/)
	  * [15 = Head and neck cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/115/)
	  * [16 = Inherited colorectal cancer (with or without polyposis)](https://panelapp.genomicsengland.co.uk/panels/254/)
	  * [17 = Inherited non-medullary thyroid cancer](https://panelapp.genomicsengland.co.uk/panels/171/)
	  * [18 = Inherited ovarian cancer (without breast cancer)](https://panelapp.genomicsengland.co.uk/panels/143/)
	  * [19 = Inherited pancreatic cancer](https://panelapp.genomicsengland.co.uk/panels/524/)
	  * [20 = Inherited renal cancer](https://panelapp.genomicsengland.co.uk/panels/521/)
	  * [21 = Inherited phaeochromocytoma and paraganglioma](https://panelapp.genomicsengland.co.uk/panels/97/)
	  * [22 = Melanoma pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/133/)
	  * [23 = Multiple endocrine tumours](https://panelapp.genomicsengland.co.uk/panels/36/)
	  * [24 = Neuroendocrine cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/183/)
	  * [25 = Ovarian cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/117/)
	  * [26 = Parathyroid Cancer](https://panelapp.genomicsengland.co.uk/panels/86/)
	  * [27 = Prostate cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/17/)
	  * [28 = Renal cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/154/)
	  * [29 = Rhabdoid tumour predisposition](https://panelapp.genomicsengland.co.uk/panels/600/)
	  * [30 = Sarcoma cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/217/)
	  * [31 = Thyroid cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/421/)
	  * [32 = Tumour predisposition - childhood onset](https://panelapp.genomicsengland.co.uk/panels/243/)
	  * [33 = Upper gastrointestinal cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/273/)


### Example report

* [Cancer predisposition sequencing report](http://folk.uio.no/sigven/example.cpsr.grch37.html)

### Docker-based technology

The CPSR workflow is developed using the [Docker technology](https://www.docker.com/what-docker). The software is thus packaged into isolated containers, in which the installation of all software libraries/tools and required dependencies have been taken care of. In addition to the bundled software, in the form of a Docker image, the workflow needs to be attached with an [annotation data bundle ](annotation_resources.html).

![](docker-logo50.png)

### Contact

sigven@ifi.uio.no
