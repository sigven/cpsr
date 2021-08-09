## Cancer Predisposition Sequencing Reporter (CPSR)

### Contents

- [Overview](#overview)
- [News](#news)
- [Example report](#example-report)
- [Annotation resources](#annotation-resources-included-in-cpsr---0.6.2)
- [CPSR Documentation](#cpsr-documentation)
- [Getting started](#getting-started)
- [Related work](#related-work)
- [Contact](#contact)


### Overview

The *Cancer Predisposition Sequencing Reporter (CPSR)* is a computational workflow that **interprets germline variants** identified from next-generation sequencing **in the context of cancer predisposition**. The workflow is integrated with the framework that underlies the [Personal Cancer Genome Reporter (PCGR)](https://github.com/sigven/pcgr), utilizing the Docker environment for encapsulation of code and software dependencies. While *PCGR* is intended for reporting and analysis of somatic variants detected in a tumor, *CPSR* is intended for reporting and ranking of germline variants in protein-coding genes that are implicated in cancer predisposition and inherited cancer syndromes.

*CPSR* accepts a query file with raw germline variant calls encoded in the [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format (i.e. analyzing SNVs/InDels). Furthermore, through the use several different _virtual cancer predisposition gene panels_ harvested from the [Genomics England PanelApp](https://panelapp.genomicsengland.co.uk/), the user can flexibly put a restriction on which genes and findings are displayed in the cancer predisposition report.


Snapshots of sections in the cancer predisposition genome report:

![CPSR views](cpsr_views.png)


The software performs extensive variant annotation on the selected geneset and produces an interactive HTML report, in which the user can investigate:

* __ClinVar variants__ - pre-classified variants according to a five-level tier scheme in ClinVar (Pathogenic to Benign)
* __Non-ClinVar variants__ - classified by CPSR through ACMG criteria (variant frequency levels and functional effects) into to a five-level tier scheme (Pathogenic to Benign)
* __Variant biomarkers__ - cancer predisposition variants with reported implications for prognosis, diagnosis or therapeutic regimens
* __Secondary findings (optional)__ - pathogenic ClinVar variants in the [ACMG recommended list for reporting of secondary findings](https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/)
* __GWAS hits (optional)__ - variants overlapping with previously identified hits in genome-wide association studies (GWAS) of cancer phenotypes (i.e. low to moderate risk conferring alleles), using [NHGRI-EBI Catalog of published genome-wide association studies](https://www.ebi.ac.uk/gwas/) as the underlying source.

The variant sets can be interactively explored and filtered further through different types of filters (phenotypes, genes, variant consequences, population MAF etc.). Importantly, the unclassified (i.e. non-ClinVar) variants are assigned a *pathogenicity score* based on the aggregation of scores according to previously established [ACMG criteria](https://www.ncbi.nlm.nih.gov/pubmed/25741868). The ACMG criteria includes cancer-specific criteria, as outlined and specified in several previous studies ([Huang et al., *Cell*, 2018](https://www.ncbi.nlm.nih.gov/pubmed/29625052); [Nykamp et al., *Genet Med.*, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28492532); [Maxwell et al., *Am J Hum Genet.*, 2016](https://www.ncbi.nlm.nih.gov/pubmed/27153395); [Amendola et al., *Am J Hum Genet.*,  2016](https://www.ncbi.nlm.nih.gov/pubmed/27181684)). See also [*Related work*](https://github.com/sigven/cpsr#related-work) below).

##### Cancer predisposition genes

The cancer predisposition report can show variants found in a number of well-known cancer predisposition genes, and the specific set of genes can be customized by the user by choosing any of the following __virtual gene panels (0 - 42)__:

  * **Panel 0** is a comprehensive, research-based _superpanel_ assembled through known sources on cancer predisposition:
	* A list of 152 genes that were curated and established within TCGA’s pan-cancer study ([Huang et al., *Cell*, 2018](https://www.ncbi.nlm.nih.gov/pubmed/29625052))
	* A list of 107 protein-coding genes that has been manually curated in COSMIC’s [Cancer Gene Census v91](https://cancer.sanger.ac.uk/census),
	* Genes from all [Genomics England PanelApp](https://panelapp.genomicsengland.co.uk/ panels for inherited cancers and tumor syndromes (detailed below)
	* Additional genes deemed relevant for cancer predisposition (contributed by the CPSR user community)

	The combination of the above sources resulted in a [non-redundant set of n = 433
	genes](https://cpsr.readthedocs.io/en/latest/superpanel.html) of relevance for cancer predisposition.

	Do you miss other genes of relevance for cancer predisposition/inherited tumor syndromes? Please forward a list of gene identifiers, preferably also with mode of inheritance and literature support to sigven AT ifi.uio.no, so we can include them in Panel 0.

* **Panels 1 - 42** are panels for inherited cancers and tumor syndromes assembled within the [Genomics England PanelApp](https://panelapp.genomicsengland.co.uk/):
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

### News
* *June 30th 2021*: **0.6.2 release**
  * Updated bundle (ClinVar, CancerMine, UniprotKB, PanelApp, CIViC, GWAS catalog)
  * Software upgrade (VEP, R/BioConductor)
  * [CHANGELOG](http://cpsr.readthedocs.io/en/latest/CHANGELOG.html)
* *November 30th 2020*: **0.6.1 release**
  * Updated bundle (ClinVar, CancerMine, UniprotKB, CIViC, GWAS catalog)
  * [CHANGELOG](http://cpsr.readthedocs.io/en/latest/CHANGELOG.html)


### Example report

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5035960.svg)](https://doi.org/10.5281/zenodo.5035960)


### Annotation resources included in CPSR - 0.6.2

* [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html) - Variant Effect Predictor v104 (GENCODE v38/v19 as the gene reference dataset), includes [gnomAD r2.1](http://gnomad.broadinstitute.org/), [dbSNP build 154](http://www.ncbi.nlm.nih.gov/SNP/), [1000 Genomes Project - phase3](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)
* [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - Database of variants with clinical significance (June 2021)
* [CIViC](http://civic.genome.wustl.edu) - clinical interpretations of variants in cancer (June 15th 2021)
* [Cancer Hotspots](http://cancerhotspots.org) - Resource for statistically significant mutations in cancer (v2 - 2017)
* [dBNSFP](https://sites.google.com/site/jpopgen/dbNSFP) - Database of non-synonymous functional predictions (v4.2, March 2021)
* [UniProt/SwissProt KnowledgeBase](http://www.uniprot.org) - Resource on protein sequence and functional information (2021_03, June 2021)
* [Pfam](http://pfam.xfam.org) - Database of protein families and domains (v34.0, March 2021)
* [CancerMine](https://zenodo.org/record/4925789#.YNT8K5MzY7Q) - Literature-derived database of tumor suppressor genes/proto-oncogenes (v36, June 2021)
* [Genomics England PanelApp](https://panelapp.genomicsengland.co.uk) - panels as of June 25th 2021
* [NHGRI-EBI GWAS catalog](https://www.ebi.ac.uk/gwas/) - GWAS catalog for cancer phenotypes, June 8th 2021)

### CPSR documentation

[![Documentation Status](https://readthedocs.org/projects/cpsr/badge/?version=latest)](https://cpsr.readthedocs.io/en/latest/?badge=latest)

**IMPORTANT**: If you use CPSR, please cite the following publication:

Sigve Nakken, Vladislav Saveliev, Oliver Hofmann, Pål Møller, Ola Myklebost, and Eivind Hovig. __Cancer Predisposition Sequencing Reporter: a flexible variant report engine for high-throughput germline screening in cancer__ (2021). _Int J Cancer_. doi:[10.1002/ijc.33749](https://doi.org/10.1002/ijc.33749)


### Getting started

#### STEP 0: Python

An installation of Python (version _3.6_) is required to run CPSR. Check that Python is installed by typing `python --version` in your terminal window.

**IMPORTANT NOTE**: STEP 1 & 2 below outline installation guidelines for running CPSR with Docker. If you want to install and run CPSR without the use of Docker (i.e. through Conda), follow [these instructions](https://github.com/sigven/cpsr/tree/master/conda_pkg/README.md)

#### STEP 1: Install PCGR (version 0.9.2)

Make sure you have a working installation of PCGR (**version 0.9.2**) and the accompanying data bundle(s) (walk through [steps 1-2](https://github.com/sigven/pcgr#getting-started)).

#### STEP 2: Download the latest release

Download the [0.6.2 release](https://github.com/sigven/cpsr/releases/tag/v0.6.2) of *cpsr* (run script)


#### STEP 3: Run example

	usage: cpsr.py -h [options]  --input_vcf <INPUT_VCF> --pcgr_dir <PCGR_DIR> --output_dir <OUTPUT_DIR> --genome_assembly  <GENOME_ASSEMBLY> --sample_id <SAMPLE_ID>

	Cancer Predisposition Sequencing Reporter - report of clinically significant cancer-predisposing germline variants

	Required arguments:
	--input_vcf INPUT_VCF
				    VCF input file with germline query variants (SNVs/InDels).
	--pcgr_dir PCGR_DIR   Directory that contains the PCGR data bundle directory, e.g. ~/pcgr-0.9.2
	--output_dir OUTPUT_DIR
				    Output directory
	--genome_assembly {grch37,grch38}
				    Genome assembly build: grch37 or grch38
	--sample_id SAMPLE_ID
				    Sample identifier - prefix for output files

	Panel options:
	--panel_id VIRTUAL_PANEL_ID
				    Comma-separated string with identifier(s) of predefined virtual cancer predisposition gene panels,
					choose any combination of the following identifiers:
				    0 = CPSR exploratory cancer predisposition panel
					(n = 335, Genomics England PanelApp / TCGA Germline Study / Cancer Gene Census / Other)
				    1 = Adult solid tumours cancer susceptibility (Genomics England PanelApp)
				    2 = Adult solid tumours for rare disease (Genomics England PanelApp)
				    3 = Bladder cancer pertinent cancer susceptibility (Genomics England PanelApp)
				    4 = Brain cancer pertinent cancer susceptibility (Genomics England PanelApp)
				    5 = Breast cancer pertinent cancer susceptibility (Genomics England PanelApp)
				    6 = Childhood solid tumours cancer susceptibility (Genomics England PanelApp)
				    7 = Colorectal cancer pertinent cancer susceptibility (Genomics England PanelApp)
				    8 = Endometrial cancer pertinent cancer susceptibility (Genomics England PanelApp)
				    9 = Familial Tumours Syndromes of the central & peripheral Nervous system (Genomics England PanelApp)
				    10 = Familial breast cancer (Genomics England PanelApp)
				    11 = Familial melanoma (Genomics England PanelApp)
				    12 = Familial prostate cancer (Genomics England PanelApp)
				    13 = Familial rhabdomyosarcoma (Genomics England PanelApp)
				    14 = GI tract tumours (Genomics England PanelApp)
				    15 = Genodermatoses with malignancies (Genomics England PanelApp)
				    16 = Haematological malignancies cancer susceptibility (Genomics England PanelApp)
				    17 = Haematological malignancies for rare disease (Genomics England PanelApp)
				    18 = Head and neck cancer pertinent cancer susceptibility (Genomics England PanelApp)
				    19 = Inherited MMR deficiency (Lynch syndrome) - Genomics England PanelApp
				    20 = Inherited non-medullary thyroid cancer (Genomics England PanelApp)
				    21 = Inherited ovarian cancer (without breast cancer) (Genomics England PanelApp)
				    22 = Inherited pancreatic cancer (Genomics England PanelApp)
				    23 = Inherited polyposis (Genomics England PanelApp)
				    24 = Inherited predisposition to acute myeloid leukaemia (AML) - Genomics England PanelApp
				    25 = Inherited predisposition to GIST (Genomics England PanelApp)
				    26 = Inherited renal cancer (Genomics England PanelApp)
				    27 = Inherited phaeochromocytoma and paraganglioma (Genomics England PanelApp)
				    28 = Melanoma pertinent cancer susceptibility (Genomics England PanelApp)
				    29 = Multiple endocrine tumours (Genomics England PanelApp)
				    30 = Multiple monogenic benign skin tumours (Genomics England PanelApp)
				    31 = Neuroendocrine cancer pertinent cancer susceptibility (Genomics England PanelApp)
				    32 = Neurofibromatosis Type 1 (Genomics England PanelApp)
				    33 = Ovarian cancer pertinent cancer susceptibility (Genomics England PanelApp)
				    34 = Parathyroid Cancer (Genomics England PanelApp)
				    35 = Prostate cancer pertinent cancer susceptibility (Genomics England PanelApp)
				    36 = Renal cancer pertinent cancer susceptibility (Genomics England PanelApp)
				    37 = Rhabdoid tumour predisposition (Genomics England PanelApp)
				    38 = Sarcoma cancer susceptibility (Genomics England PanelApp)
				    39 = Sarcoma susceptibility (Genomics England PanelApp)
				    40 = Thyroid cancer pertinent cancer susceptibility (Genomics England PanelApp)
				    41 = Tumour predisposition - childhood onset (Genomics England PanelApp)
				    42 = Upper gastrointestinal cancer pertinent cancer susceptibility (Genomics England PanelApp)

	--custom_list CUSTOM_LIST
				    Provide custom list of genes from virtual panel 0 (single-column txt file with Ensembl gene identifiers),
					alternative to predefined panels provided with --panel_id)
	--custom_list_name CUSTOM_LIST_NAME
				    Set name for custom made panel/list (single word - no whitespace), will be displayed in the report
	--diagnostic_grade_only
				    For panel_id's 1-42 (Genomics England PanelApp) - consider genes with a GREEN status only, default: False

	VEP options:
	--vep_n_forks VEP_N_FORKS
				    Number of forks (option '--fork' in VEP), default: 4
	--vep_buffer_size VEP_BUFFER_SIZE
				    Variant buffer size (variants read into memory simultaneously, option '--buffer_size' in VEP)
				    - set lower to reduce memory usage, default: 500
	--vep_pick_order VEP_PICK_ORDER
				    Comma-separated string of ordered transcript properties for primary variant pick
					( option '--pick_order' in VEP), default: canonical,appris,biotype,ccds,rank,tsl,length,mane
	--vep_no_intergenic   Skip intergenic variants during processing (option '--no_intergenic' in VEP), default: False

	vcfanno options:
	--vcfanno_n_proc VCFANNO_N_PROC
				    Number of vcfanno processes (option '-p' in vcfanno), default: 4

	Other options:
	--force_overwrite     By default, the script will fail with an error if any output file already exists.
					You can force the overwrite of existing result files by using this flag, default: False
	--version             show program's version number and exit
	--basic               Run functional variant annotation on VCF through VEP/vcfanno, omit Tier assignment/report generation (STEP 4), default: False
	--no_vcf_validate     Skip validation of input VCF with Ensembl's vcf-validator, default: False
	--docker_uid DOCKER_USER_ID
				    Docker user ID. Default is the host system user ID. If you are experiencing permission errors,
					try setting this up to root (`--docker_uid root`), default: None
	--no_docker           Run the CPSR workflow in a non-Docker mode, default: False
	--preserved_info_tags PRESERVED_INFO_TAGS
				    Comma-separated string of VCF INFO tags from query VCF that should be kept in CPSR output TSV
	--report_theme {default,cerulean,journal,flatly,readable,spacelab,united,cosmo,lumen,paper,sandstone,simplex,yeti}
				    Visual report theme (rmarkdown), default: default
	--report_nonfloating_toc
				    Do not float the table of contents (TOC) in output HTML report, default: False
	--report_table_display {full,light}
				    Set the level of detail/comprehensiveness in interactive datables of HTML report, very comprehensive (option 'full') or slim/focused ('light')
	--ignore_noncoding    Do not list non-coding variants in HTML report, default: False
	--secondary_findings  Include variants found in ACMG-recommended list for secondary findings (v3.0), default: False
	--gwas_findings       Report overlap with low to moderate cancer risk variants (tag SNPs) identified from genome-wide association studies, default: False
	--gwas_p_value GWAS_P_VALUE
				    Required p-value for variants listed as hits from genome-wide association studies, default: 5e-06
	--pop_gnomad {afr,amr,eas,sas,asj,nfe,fin,global}
				    Population source in gnomAD used for variant frequency assessment (ACMG classification), default: nfe
	--maf_upper_threshold MAF_UPPER_THRESHOLD
				    Upper MAF limit (gnomAD global population frequency) for variants to be included in the report, default: 0.9
	--classify_all        Provide CPSR variant classifications (TIER 1-5) also for variants with exising ClinVar classifications in output TSV, default: False
	--clinvar_ignore_noncancer
				    Ignore (exclude from report) ClinVar-classified variants reported only for phenotypes/conditions NOT related to cancer, default: False
	--debug            Print full docker commands to log, default: False


The *cpsr* software bundle contains an example VCF file.

Report generation with the example VCF, using the [Adult solid tumours cancer susceptibility](https://panelapp.genomicsengland.co.uk/panels/245/) virtual gene panel, can be performed through the following command:

	python ~/cpsr-0.6.2/cpsr.py
	 --input_vcf ~/cpsr-0.6.2/example.vcf.gz
	 --pcgr_dir ~/pcgr-0.9.2
	 --output_dir ~/cpsr-0.6.2
	 --genome_assembly grch37
	 --panel_id 1
	 --sample_id example
	 --secondary_findings
	 --classify_all
	 --maf_upper_threshold 0.2
	 --no_vcf_validate

Note that the example command also refers to the PCGR directory (*pcgr-0.9.2*), which contains the data bundle that are necessary for both *PCGR* and *CPSR*.

This command will run the Docker-based *cpsr* workflow and produce the following output files in the _output_ folder:

  1. __example.cpsr.grch37.vcf.gz (.tbi)__ - Bgzipped VCF file with relevant annotations appended by CPSR
  2. __example.cpsr.grch37.pass.vcf.gz (.tbi)__ - Bgzipped VCF file with relevant annotations appended by CPSR (PASS variants only)
  3. __example.cpsr_config.rds__ - CPSR configuration object (RDS format), mostly for debugging purposes
  4. __example.cpsr.grch37.pass.tsv.gz__ - Compressed TSV file (generated with [vcf2tsv](https://github.com/sigven/vcf2tsv)) of VCF content with relevant annotations appended by CPSR
  5. __example.cpsr.grch37.html__ - Interactive HTML report with clinically relevant variants in cancer predisposition genes organized into tiers
  6. __example.cpsr.grch37.json.gz__ - Compressed JSON dump of HTML report content
  7. __example.cpsr.grch37.snvs_indels.tiers.tsv__ - TSV file with key annotations of tier-structured SNVs/InDels

### Related work

* [CharGer - Characterization of Germline variants](https://github.com/ding-lab/CharGer)
* [PathoMan - Pathogenicity of Mutation Analyzer (Beta)](https://pathoman.mskcc.org/)
* [SherLoc - Variant classification](https://www.invitae.com/en/variant-classification/)

### Contact

sigven AT ifi.uio.no
