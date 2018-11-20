## About

###  What is the Cancer Predisposition Sequencing Preporter (CPSR)?

The *Cancer Predisposition Sequencing Reporter (CPSR)* is a computational workflow that **interprets germline variants** identified from next-generation sequencing **in the context of cancer predisposition**. The workflow is integrated with the framework that underlies the [Personal Cancer Genome Reporter (PCGR)](https://github.com/sigven), utilizing the Docker environment for encapsulation of code and software dependencies. While *PCGR* is intended for reporting and analysis of somatic variants detected in a tumor, *CPSR* is intended for reporting and ranking of germline variants in protein-coding genes that are implicated in cancer predisposition and inherited cancer syndromes.


*CPSR* accepts a query file with raw germline variant calls encoded in the [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format (i.e. analyzing SNVs/InDels). The software performs extensive variant annotation and produces an interactive HTML report, in which the user can investigate three main sets of variants identified in the query set:

1. Germline variants in a selected set of [configurable cancer predisposition genes](predisposition.md), that are **previously reported** as pathogenic or likely pathogenic in ClinVar (with no conflicting interpretations)

2. **Unclassified variants** constitute the set of germline variants within the configurable cancer predisposition gene list that are either:
	* Registered as *variant of uncertain significance (VUS)* in [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), or
	* *Is a novel protein-coding variant* (i.e. not reported in [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), and not found in [gnomAD](http://gnomad.broadinstitute.org/) or [1000 Genomes Project](http://www.internationalgenome.org/) user-defined population datasets), or
	* *Is a rare protein-coding variant* (e.g. minor allele frequency (MAF) < 0.001 in user-defined [gnomAD](http://gnomad.broadinstitute.org/) or [1000 Genomes Project](http://www.internationalgenome.org/) population datasets)
		* *The upper MAF threshold (e.g. 0.001) for listing of unclassified variants can be configured by the user*


3. Variants overlapping with previously identified hits in genome-wide association studies (GWAS) of cancer phenotypes (i.e. low to moderate risk conferring alleles), using [NHGRI-EBI Catalog of published genome-wide association studies](https://www.ebi.ac.uk/gwas/) as the underlying source.

The (**classified** and **unclassified**) variant sets can be interactively explored and ranked further through different types of filters (associated phenotypes, genes, variant consequences, population MAF etc.). Importantly, the unclassified variants are assigned and ranked according to a *pathogenicity score*, which is based on the aggregation of scores according to previously established [ACMG criteria](https://www.ncbi.nlm.nih.gov/pubmed/25741868) and also cancer-specific criteria, as outlined and specified in several previous studies ([Huang et al., *Cell*, 2018](https://www.ncbi.nlm.nih.gov/pubmed/29625052); [Maxwell et al., *Am J Hum Genet.*, 2016](https://www.ncbi.nlm.nih.gov/pubmed/27153395); [Amendola et al., *Am J Hum Genet.*,  2016](https://www.ncbi.nlm.nih.gov/pubmed/27181684)).

#### Cancer predisposition genes

We have compiled a comprehensive list of genes that are implicated in cancer predisposition and cancer syndromes. Three different sources were combined:
* A list of 152 genes that were curated and established within TCGA’s pan-cancer study ([Huang et al., *Cell*, 2018](https://www.ncbi.nlm.nih.gov/pubmed/29625052))
* A list of 107 protein-coding genes that has been manually curated in COSMIC’s [Cancer Gene Census v86](https://cancer.sanger.ac.uk/census)
* A list of 148 protein-coding genes established by experts within the Norwegian Cancer Genomics Consortium (http://cancergenomics.no)

The combination of the three sources resulted in a non-redundant set of 209 protein-coding genes of relevance for predisposition to tumor development. We want to make it explicit that this list of 209 genes is by no means regarded as an international consensus, but should rather be subject to continuous update by the international community that carry expertise on genetic risk factors for cancer.

The final list can be downloaded as a tab-separated text file [here](https://raw.githubusercontent.com/sigven/cpsr/master/predisposition_genes_20181112.tsv).


The Cancer Predisposition Sequencing Reporter has been developed by scientists affiliated with the [Norwegian Cancer Genomics Consortium](http://cancergenomics.no), at the [Institute for Cancer Research/Oslo University Hospital](http://radium.no).

### Example report
* [Cancer predisposition sequencing report](http://folk.uio.no/sigven/example.cpsr.grch37.html)

### Docker-based technology

The PCGR workflow is developed using the [Docker technology](https://www.docker.com/what-docker). The software is thus packaged into isolated containers, in which the installation of all software libraries/tools and required dependencies have been taken care of. In addition to the bundled software, in the form of a Docker image, the workflow only needs to be attached with an [annotation data bundle for precision oncology](annotation_resources.html).

![](docker-logo50.png)

### Contact

sigven@ifi.uio.no
