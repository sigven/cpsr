## Getting started


### STEP 0: Python

An installation of Python (version _3.6_) is required to run CPSR. Check that Python is installed by typing `python --version` in your terminal window. In addition, a [Python library](https://github.com/uiri/toml) for parsing configuration files encoded with [TOML](https://github.com/toml-lang/toml) is needed. To install, simply run the following command:

	pip install toml

### STEP 1: Installation of Docker

1. [Install the Docker engine](https://docs.docker.com/engine/installation/) on your preferred platform
	- installing [Docker on Linux](https://docs.docker.com/engine/installation/linux/)
	- installing [Docker on Mac OS](https://docs.docker.com/engine/installation/mac/)
	- NOTE: We have not yet been able to perform enough testing on the Windows platform, and we have received feedback that particular versions of Docker/Windows do not work with CPSR (an example being [mounting of data volumes](https://github.com/docker/toolbox/issues/607))
2. Test that Docker is running, e.g. by typing `docker ps` or `docker images` in the terminal window
3. Adjust the computing resources dedicated to the Docker, i.e.:
	- Memory: minimum 5GB
	- CPUs: minimum 4
	- [How to - Mac OS X](https://docs.docker.com/docker-for-mac/#advanced)

### STEP 2: Download run script/data bundle, and pull Docker image

1. Download and unpack the [CPSR pre-release (0.3.0)](https://github.com/sigven/cpsr/releases/tag/v0.3.0)
2. Pull the PCGR Docker image (*dev version*): `docker pull sigven/pcgr:dev`
3. Download and unpack the latest PCGR data bundles
	* [grch37 data bundle - 20181119](https://drive.google.com/open?id=1OL5C994HDaeadASz7KzMhPoXfdSiyhNy) (approx 9Gb)
	* [grch38 data bundle - 20181119](https://drive.google.com/open?id=1CZNc87E0K5AK2RDSNU57FqLp0H1skpUh) (approx 14Gb)
	* *Unpacking*: `gzip -dc pcgr.databundle.grch37.YYYYMMDD.tgz | tar xvf -`

 	A _data/_ folder should now have been produced

### STEP 3: Input preprocessing

The CPSR workflow accepts one type of input file

* An unannotated, single-sample VCF file (>= v4.2) with germline DNA variants (SNVs/InDels)

* We __strongly__ recommend that the input VCF is compressed and indexed using [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html)
* If the input VCF contains multi-allelic sites, these will be subject to [decomposition](http://genome.sph.umich.edu/wiki/Vt#Decompose)
* Variants used for reporting should be designated as 'PASS' in the VCF FILTER column (non-PASS variants are ignored in the report)

### STEP 4: Configuration

A few elements of the workflow can be figured using the *cpsr* configuration file, encoded in [TOML](https://github.com/toml-lang/toml) (an easy to read file format).

The initial step of the workflow performs [VCF validation](https://github.com/EBIvariation/vcf-validator) on the input VCF file. This procedure is very strict, and often causes the workflow to return an error due to various violations of the VCF specification. If the user trusts that the most critical parts of the input VCF is properly encoded,  a setting in the configuration file (`vcf_validation = false`) can be used to turn off VCF validation.

An exhaustive, predefined list of 209 cancer predisposition/syndrome genes can also be configured.

See section on [Input](input.md) for more details wrt. configuration.

### STEP 5: Run example

Run the workflow with **cpsr.py**, which takes the following arguments and options:

	usage: cpsr.py [-h] [--input_vcf INPUT_VCF] [--force_overwrite] [--version]
			[--basic] [--docker-uid DOCKER_USER_ID] [--no-docker]
			pcgr_base_dir output_dir {grch37,grch38} configuration_file
			sample_id

	Cancer Predisposition Sequencing Reporter (CPSR) - report of cancer-predisposing
	germline variants

	positional arguments:
	pcgr_base_dir         Directory that contains the PCGR data bundle
				    directory, e.g. ~/pcgr-dev
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
	--docker-uid DOCKER_USER_ID
				    Docker user ID. Default is the host system user ID. If
				    you are experiencing permission errors, try setting
				    this up to root (`--docker-uid root`) (default: None)
	--no-docker           Run the CPSR workflow in a non-Docker mode (see
				    install_no_docker/ folder for instructions (default:
				    False)


The *cpsr* software bundle contains an example VCF file. It also contains a configuration file (*cpsr.toml*).

Analysis of the example VCF can be performed by the following command:

`python ~/cpsr-0.3.0/cpsr.py --input_vcf ~/cpsr-0.3.0/example.vcf.gz`
` ~/pcgr-dev ~/cpsr-0.3.0 grch37 ~/cpsr-0.3.0/cpsr.toml example`

Note that the example command also refers to the PCGR directory (*pcgr-dev*), which contains the data bundle that are necessary for both *PCGR* and *CPSR*.

This command will run the Docker-based *cpsr* workflow and produce the following output files in the _cpsr_ folder:

  1. __example.cpsr.grch37.pass.vcf.gz (.tbi)__ - Bgzipped VCF file with functional/clinical annotations
  2. __example.cpsr.grch37.pass.tsv.gz__ - Compressed TSV file (generated with [vcf2tsv](https://github.com/sigven/vcf2tsv)) with functional/clinical annotations
  3. __example.cpsr.grch37.html__ - Interactive HTML report with clinically relevant variants in cancer predisposition genes organized into tiers
  4. __example.cpsr.grch37.json.gz__ - Compressed JSON dump of HTML report content
  5. __example.cpsr.snvs_indels.tiers.grch37.tsv__ - TSV file with most important annotations of tier-structured SNVs/InDels
