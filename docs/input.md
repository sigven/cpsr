## Input

The CPSR workflow accepts a single input file:

  * An unannotated, single-sample **VCF file** (>= v4.2) with germline calls (SNVs/InDels)

### VCF

* We __strongly__ recommend that the input VCF is compressed and indexed using [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html)
* If the input VCF contains multi-allelic sites, these will be subject to [decomposition](http://genome.sph.umich.edu/wiki/Vt#Decompose)
* Variants used for reporting should be designated as 'PASS' in the VCF FILTER column

__IMPORTANT NOTE __: CPSR generates a number of VCF INFO annotation tags that is appended to the query VCF. We will therefore encourage the users to submit query VCF files that have not been subject to annotations by other means, but rather a VCF file that comes directly from variant calling. If not, there are likely to be INFO tags in the query VCF file that coincide with those produced by CPSR.


### CPSR configuration file

The cancer predisposition sequencing report can be flexibly configured in a TOML-formatted configuration file. The default TOML configuration file, are shown below:

	# CPSR configuration options (TOML).

	[popgen]
	## choose population source in gnomAD (non-cancer subset for use as control), defaults to the global set

	## For gnomaAD, this can by any of the following values (three-letter codes):
	## "afr" - African/American (12,020 individuals (7,652 WES / 4,368 WGS))
	## "amr" - Admixed American (17,210 individuals (16,791 WES / 419 WGS))
	## "eas" - East Asian (9,435 individuals (8,624 WES / 811 WGS))
	## "sas" - Sout Asian (15,391 individuals (15,391 WES / 0 WGS))
	## "nfe" - Non-Finnish European (63,369 individuals (55,860 WES / 7,509 WGS))
	## "fin" - Finnish (12,897 individuals (11,150 WES / 1,747 WGS))
	## "global" - All populations (138,632 individuals (123,136 WES / 15,496 WGS))
	pop_gnomad = "global"

	[visual]
	# Choose visual theme of report, any of: "default", "cerulean", "journal", "flatly", "readable", "spacelab", "united", "cosmo", "lumen", "paper", "sandstone", "simplex", or "yeti" (https://bootswatch.com/)
	report_theme = "default"

	[custom_tags]
	## list VCF info tags that should be present in JSON output
	## tags should be comma separated, i.e. custom_tags = "GATK_FILTER,VARSCAN_FILTER"
	custom_tags = ""

	[custom_panel]
	## edit metadata (name, version, url) of relevance for custom-provided gene panel (provided by option --custom_panel)
	name = "Custom gene panel"
	version = "1.0"
	url = ""

	[gwas]
	## Required p-value for reporting of GWAS hits
	p_value_min = 5e-6


	[other]
	n_vcfanno_proc = 4
	n_vep_forks = 4
	vep_skip_intergenic = false
	## choice of how VEP chooses the primary transcript pr. gene
	## https://www.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options
	vep_pick_order = "canonical,appris,tsl,biotype,ccds,rank,length,mane"
