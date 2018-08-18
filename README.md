## Cancer Predisposition Sequencing Report (CPSR)

### Overview

The *Cancer Predisposition Sequencing Report (CPSR)* is a computational workflow that **interprets germline variants** identified from next-generation sequencing **in the context of cancer predisposition**. The workflow is integrated with the framework that underlies the [Personal Cancer Genome Reporter (PCGR)](https://github.com/sigven), utilizing the Docker environment for encapsulation of code and software dependencies. While *PCGR* is intended for reporting and analysis of somatic variants detected in a tumor, *CPSR* is intended for reporting and ranking of germline variants in protein-coding genes that are implicated in cancer predisposition and inherited cancer syndromes.

*CPSR* accepts a query file encoded in the [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format (i.e. analyzing SNVs and InDels). The software performs extensive variant annotation and produces an interactive HTML report, in which the user can investigate two main sets of variants (in a selected set of configurable cancer predisposition genes):

1. Germline variants that are **previously reported** as pathogenic/likely pathogenic/uncertain significance in ClinVar (with no conflicting interpretations). This set is organized into three tiers:
	* Tier 1 - Pathogenic variants
	* Tier 2 - Likely pathogenic variants
	* Tier 3 - Variants of uncertain significance (VUS)


2. **Unclassified** (i.e. **not reported** in [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)) protein-coding germline variants that are either:
	* *Novel* (i.e. not reported in [gnomAD](http://gnomad.broadinstitute.org/) or [1000 Genomes Project](http://www.internationalgenome.org/) global population datasets), or
	* *Rare* (i.e. minor allele frequency (MAF) < 0.001 in the [gnomAD](http://gnomad.broadinstitute.org/) or [1000 Genomes Project](http://www.internationalgenome.org/) global population datasets)
		* *The upper MAF threshold (e.g. 0.001) for listing of unclassified variants can be configured by the user*

Both variant sets (**previously reported/classified** and **unclassified**) can be interactively explored and ranked further through different types of filters (associated phenotypes, genes, variant consequences, population MAF etc.)

### Example report

* [Cancer predisposition sequencing report](http://folk.uio.no/sigven/example.cpsr.grch37.html)

### Annotation resources included in _cpsr - 0.1.0_

* [VEP v93](http://www.ensembl.org/info/docs/tools/vep/index.html) - Variant Effect Predictor release 93 (GENCODE v19/v28 as the gene reference dataset)
* [dBNSFP v3.5](https://sites.google.com/site/jpopgen/dbNSFP) - Database of non-synonymous functional predictions (August 2017)
* [gnomAD r2](http://gnomad.broadinstitute.org/) - Germline variant frequencies exome-wide (October 2017)
* [dbSNP b150](http://www.ncbi.nlm.nih.gov/SNP/) - Database of short genetic variants (February 2017)
* [1000 Genomes Project - phase3](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) - Germline variant frequencies genome-wide (May 2013)
* [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - Database of clinically related variants (August 2018)
* [DisGeNET](http://www.disgenet.org) - Database of gene-disease associations (v5.0, May 2017)
* [UniProt/SwissProt KnowledgeBase 2018_07](http://www.uniprot.org) - Resource on protein sequence and functional information (July 2018)
* [Pfam v31](http://pfam.xfam.org) - Database of protein families and domains (March 2017)
* [TSGene v2.0](http://bioinfo.mc.vanderbilt.edu/TSGene/) - Tumor suppressor/oncogene database (November 2015)

### News

*COMING SOON: 0.1.0 release (accompanied by PCGR 0.6.3 release)*

### Getting started

#### STEP 0: Install PCGR

Make sure you have a working installation of the latest PCGR release (0.6.3) and latest data bundle (walk through [steps 0-2](https://github.com/sigven/pcgr#getting-started)).

#### STEP 1: Download the latest release

Download the [latest release](https://github.com/sigven/releases/) of *cpsr* (run script and configuration file)

#### STEP 2: Configuration

A few elements of the workflow can be figured using the *cpsr* configuration file, encoded in [TOML](https://github.com/toml-lang/toml) (an easy to read file format).

The initial step of the workflow performs [VCF validation](https://github.com/EBIvariation/vcf-validator) on the input VCF file. This procedure is very strict, and often causes the workflow to return an error due to various violations of the VCF specification. If the user trusts that the most critical parts of the input VCF is properly encoded,  a setting in the configuration file (`vcf_validation = false`) can be used to turn off VCF validation.

An exhaustive, pre-defined list of 212 cancer predisposition genes can also be configured.

#### STEP 3: Run example

Run the workflow with **cpsr.py**, which takes the following arguments and options:

	usage: cpsr.py [-h] [--input_vcf INPUT_VCF] [--force_overwrite]
					 [--version] [--basic]
					 data_dir output_dir {grch37,grch38}
					 configuration_file sample_id

	Cancer Predisposition Sequencing Report (CPSR)

	positional arguments:
	data_dir         Directory that contains the PCGR data bundle
				    directory, e.g. ~/pcgr-0.6.0
	output_dir            Output directory
	{grch37,grch38}       Genome assembly build: grch37 or grch38
	configuration_file    Configuration file (TOML format)
	sample_id             Sample identifier - prefix for output files

	optional arguments:
	-h, --help            show this help message and exit
	--input_vcf INPUT_VCF
				    VCF input file with somatic query variants
				    (SNVs/InDels). (default: None)
	--force_overwrite     By default, the script will fail with an error if any
				    output file already exists. You can force the
				    overwrite of existing result files by using this flag
				    (default: False)
	--version             show program's version number and exit
	--basic               Run functional variant annotation on VCF through
				    VEP/vcfanno, omit report generation (STEP 4) (default:
				    False)



The *cpsr* software bundle contains an example VCF file. It also contains a configuration file (*cpsr.toml*). **NOTE: The example file contains a significant number of known pathogenic and likely pathogenic variants. This is for demonstration purposes only, and NOT a realistic situation in a regular sample.**

Analysis of the example VCF can be performed by the following command:

`python ~/cpsr-0.1.0/cpsr.py --input_vcf ~/cpsr-0.1.0/example.vcf.gz`
` ~/pcgr-0.6.3 ~/cpsr-0.1.0 grch37 ~/cpsr-0.1.0/cpsr.toml example`

Note that the example command also refers to the PCGR directory (*pcgr-0.6.3*), which contains the data bundle that are necessary for both *PCGR* and *CPSR*.

This command will run the Docker-based *cpsr* workflow and produce the following output files in the _cpsr_ folder:

  1. __example.cpsr.grch37.pass.vcf.gz (.tbi)__ - Bgzipped VCF file with functional/clinical annotations
  2. __example.cpsr.grch37.pass.tsv.gz__ - Compressed TSV file (generated with [vcf2tsv](https://github.com/sigven/vcf2tsv)) with functional/clinical annotations
  3. __example.cpsr.grch37.html__ - Interactive HTML report with clinically relevant variants in cancer predisposition genes
  4. __example.cpsr.grch37.json.gz__ - Compressed JSON dump of HTML report content


### Contact

sigven@ifi.uio.no
