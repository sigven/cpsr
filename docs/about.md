## About

###  What is the Cancer Predisposition Sequencing Reporter (CPSR)?

The *Cancer Predisposition Sequencing Reporter (CPSR)* is a computational workflow that **interprets germline variants** identified from next-generation sequencing **in the context of cancer predisposition**. The workflow is integrated with the framework that underlies the [Personal Cancer Genome Reporter (PCGR)](https://github.com/sigven), utilizing the Docker environment for encapsulation of code and software dependencies. While *PCGR* is intended for reporting and analysis of somatic variants detected in a tumor, *CPSR* is intended for reporting and ranking of germline variants in protein-coding genes that are implicated in cancer predisposition and inherited cancer syndromes.

*CPSR* accepts a query file with raw germline variant calls encoded in the [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format (i.e. analyzing SNVs/InDels). Furthermore, through the use several different _virtual cancer predisposition gene panels_ harvested from the [Genomics England PanelApp](https://panelapp.genomicsengland.co.uk/), the user can flexibly put a restriction on which genes and findings are displayed in the cancer predisposition report.


Snapshots of sections in the cancer predisposition genome report:

![CPSR views](cpsr_views.png)


The software performs extensive variant annotation on the selected geneset and produces an interactive HTML report, in which the user can investigate:

* __ClinVar variants__ - pre-classified variants according to a five-level tier scheme in ClinVar (Pathogenic to Benign)
* __Non-ClinVar variants__ - classified by CPSR through ACMG criteria (variant frequency levels and functional effects) into to a five-level tier scheme (Pathogenic to Benign)
* __Genomic biomarkers__ - cancer predisposition variants with reported implications for prognosis, diagnosis or therapeutic regimens
* __Secondary findings (optional)__ - pathogenic ClinVar variants in the [ACMG recommended list for reporting of incidental findings](https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/)
* __GWAS hits (optional)__ - variants overlapping with previously identified hits in genome-wide association studies (GWAS) of cancer phenotypes (i.e. low to moderate risk conferring alleles), using [NHGRI-EBI Catalog of published genome-wide association studies](https://www.ebi.ac.uk/gwas/) as the underlying source.

The variant sets can be interactively explored and filtered further through different types of filters (phenotypes, genes, variant consequences, population MAF etc.). Importantly, the unclassified non-ClinVar variants are assigned a *pathogenicity level* based on the aggregation of scores according to previously established [ACMG criteria](https://www.ncbi.nlm.nih.gov/pubmed/25741868). The ACMG criteria includes cancer-specific criteria, as outlined and specified in several previous studies ([Huang et al., *Cell*, 2018](https://www.ncbi.nlm.nih.gov/pubmed/29625052); [Nykamp et al., *Genet Med.*, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28492532); [Maxwell et al., *Am J Hum Genet.*, 2016](https://www.ncbi.nlm.nih.gov/pubmed/27153395); [Amendola et al., *Am J Hum Genet.*,  2016](https://www.ncbi.nlm.nih.gov/pubmed/27181684)). See also [*Related work*](https://github.com/sigven/cpsr#related-work)

##### Cancer predisposition genes

The cancer predisposition report can show variants found in a number of well-known cancer predisposition genes, and the specific set of genes can be customized by the user by choosing any of the following __virtual gene panels (0 - 42)__:

  * **Panel 0 (default)** is a comprehensive, research-based gene panel assembled through known sources on cancer predisposition:
	* A list of 152 genes that were curated and established within TCGA’s pan-cancer study ([Huang et al., *Cell*, 2018](https://www.ncbi.nlm.nih.gov/pubmed/29625052))
	* A list of 107 protein-coding genes that has been manually curated in COSMIC’s [Cancer Gene Census v91](https://cancer.sanger.ac.uk/census),
	* A list of 148 protein-coding genes established by experts within the Norwegian Cancer Genomics Consortium (http://cancergenomics.no)
	* Additional genes deemed relevant for cancer predisposition (contributed by the CPSR user community)

	The combination of the above sources resulted in a non-redundant set of [216 protein-coding genes](https://github.com/sigven/cpsr/blob/master/predisposition.md) which carry variants that may predispose to tumor development.

	* **Panels 1 - 42** are panels for inherited cancer syndromes and cancer predisposition assembled within the [Genomics England PanelApp](https://panelapp.genomicsengland.co.uk/):
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
       * [14 = GI tract tumours](https://panelapp.genomicsengland.co.uk/panels/254/)
       * [15 = Genodermatoses with malignancies](https://panelapp.genomicsengland.co.uk/panels/201/)
       * [16 = Haematological malignancies cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/59/)
       * [17 = Haematological malignancies for rare disease](https://panelapp.genomicsengland.co.uk/panels/407/)
       * [18 = Head and neck cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/115/)
	  * [19 = Inherited MMR deficiency (Lynch syndrome)](https://panelapp.genomicsengland.co.uk/panels/503/)
       * [20 = Inherited non-medullary thyroid cancer](https://panelapp.genomicsengland.co.uk/panels/171/)
       * [21 = Inherited ovarian cancer (without breast cancer)](https://panelapp.genomicsengland.co.uk/panels/143/)
       * [22 = Inherited pancreatic cancer](https://panelapp.genomicsengland.co.uk/panels/524/)
	  * [23 = Inherited polyposis](https://panelapp.genomicsengland.co.uk/panels/504/)
	  * [24 = Inherited predisposition to acute myeloid leukaemia (AML)](https://panelapp.genomicsengland.co.uk/panels/525/)
	  * [25 = Inherited predisposition to GIST](https://panelapp.genomicsengland.co.uk/panels/523/)
       * [26 = Inherited renal cancer](https://panelapp.genomicsengland.co.uk/panels/521/)
       * [27 = Inherited phaeochromocytoma and paraganglioma](https://panelapp.genomicsengland.co.uk/panels/97/)
       * [28 = Melanoma pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/133/)
       * [29 = Multiple endocrine tumours](https://panelapp.genomicsengland.co.uk/panels/36/)
       * [30 = Multiple monogenic benign skin tumours](https://panelapp.genomicsengland.co.uk/panels/558/)
       * [31 = Neuroendocrine cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/183/)
       * [32 - Neurofibromatosis Type 1](https://panelapp.genomicsengland.co.uk/panels/255/)
       * [33 = Ovarian cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/117/)
       * [34 = Parathyroid Cancer](https://panelapp.genomicsengland.co.uk/panels/86/)
       * [35 = Prostate cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/17/)
       * [36 = Renal cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/154/)
       * [37 = Rhabdoid tumour predisposition](https://panelapp.genomicsengland.co.uk/panels/600/)
       * [38 = Sarcoma cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/217/)
       * [39 = Sarcoma susceptibility](https://panelapp.genomicsengland.co.uk/panels/734/)
       * [40 = Thyroid cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/421/)
       * [41 = Tumour predisposition - childhood onset](https://panelapp.genomicsengland.co.uk/panels/243/)
       * [42 = Upper gastrointestinal cancer pertinent cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/273/)

### Example report

* [Cancer predisposition genome report](http://insilico.hpc.uio.no/pcgr/example_reports/cpsr/0.6.0rc/SAMPLE-001.cpsr.grch37.html)

### Docker-based technology

The CPSR workflow is developed using the [Docker technology](https://www.docker.com/what-docker). The software is thus packaged into isolated containers, in which the installation of all software libraries/tools and required dependencies have been taken care of. In addition to the bundled software, in the form of a Docker image, the workflow needs to be attached with an [annotation data bundle ](annotation_resources.html).

![](docker-logo50.png)

### Contact

sigven@ifi.uio.no
