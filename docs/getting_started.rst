Getting started
---------------

STEP 0: Python
~~~~~~~~~~~~~~

An installation of Python (version *3.6*) is required to run CPSR. Check
that Python is installed by typing ``python --version`` in your terminal
window. In addition, a `Python library <https://github.com/uiri/toml>`__
for parsing configuration files encoded with
`TOML <https://github.com/toml-lang/toml>`__ is needed. To install,
simply run the following command:

::

   pip install toml

STEP 1: Installation of Docker
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. `Install the Docker
   engine <https://docs.docker.com/engine/installation/>`__ on your
   preferred platform

   -  installing `Docker on
      Linux <https://docs.docker.com/engine/installation/linux/>`__
   -  installing `Docker on Mac
      OS <https://docs.docker.com/engine/installation/mac/>`__
   -  NOTE: We have not yet been able to perform enough testing on the
      Windows platform, and we have received feedback that particular
      versions of Docker/Windows do not work with CPSR (an example being
      `mounting of data
      volumes <https://github.com/docker/toolbox/issues/607>`__)

2. Test that Docker is running, e.g. by typing ``docker ps`` or
   ``docker images`` in the terminal window
3. Adjust the computing resources dedicated to the Docker, i.e.:

   -  Memory: minimum 5GB
   -  CPUs: minimum 4
   -  `How to - Mac OS
      X <https://docs.docker.com/docker-for-mac/#advanced>`__

STEP 2: Download run script/data bundle, and pull Docker image
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Download and unpack the latest `CPSR release
   (0.5.2) <https://github.com/sigven/cpsr/releases/tag/v0.5.2>`__
2. Pull the latest PCGR Docker image (*0.8.4*):
   ``docker pull sigven/pcgr:0.8.4``
3. Download and unpack the latest PCGR data bundles

   -  `grch37 data bundle -
      20191116 <https://drive.google.com/file/d/1TdYagetk-l__aYBsaZJHJvYFStDnIEcq>`__
      (approx 16Gb)
   -  `grch38 data bundle -
      20191116 <https://drive.google.com/file/d/1wpVqlgY5jBKkQaTAxzgf0rgxKxOEJgj->`__
      (approx 17Gb)
   -  *Unpacking*:
      ``gzip -dc pcgr.databundle.grch37.YYYYMMDD.tgz | tar xvf -``

   A *data/* folder should now have been produced

STEP 3: Input preprocessing
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The CPSR workflow accepts one type of input file

-  An unannotated, single-sample VCF file (>= v4.2) with germline DNA
   variants (SNVs/InDels)

-  We **strongly** recommend that the input VCF is compressed and
   indexed using `bgzip <http://www.htslib.org/doc/tabix.html>`__ and
   `tabix <http://www.htslib.org/doc/tabix.html>`__
-  If the input VCF contains multi-allelic sites, these will be subject
   to `decomposition <http://genome.sph.umich.edu/wiki/Vt#Decompose>`__

-  Variants used for reporting should be designated as ‘PASS’ in the VCF
   FILTER column (non-PASS variants are simply ignored in the report)

STEP 4: Configuration
~~~~~~~~~~~~~~~~~~~~~

A few elements of the workflow can be figured using the *cpsr*
configuration file, encoded in
`TOML <https://github.com/toml-lang/toml>`__. The following can be
configured:

-  Choice of gnomAD control population
-  Upper MAF limit for variants considered for inclusion in the report
-  Inclusion of GWAS hits
-  Inclusion of secondary findings
-  Inclusion of ACMG classifications for ClinVar variants
-  VEP/\ *vcfanno* options

See section on `Input <input.html>`__ for more details wrt. default
configuration.

STEP 5: Run example
~~~~~~~~~~~~~~~~~~~

::

   Cancer Predisposition Sequencing Reporter (CPSR) - report of cancer-predisposing germline variants

   positional arguments:
   query_vcf             VCF input file with germline query variants (SNVs/InDels).
   pcgr_base_dir         Directory that contains the PCGR data bundle directory, e.g. ~/pcgr-0.8.4
   output_dir            Output directory
   {grch37,grch38}       Genome assembly build: grch37 or grch38
   configuration_file    Configuration file (TOML format)
   sample_id             Sample identifier - prefix for output files

   optional arguments:
   -h, --help            show this help message and exit
   --force_overwrite     By default, the script will fail with an error if any output file already exists.
                   You can force the overwrite of existing result files by using this flag
   --version             show program's version number and exit
   --basic               Run functional variant annotation on VCF through VEP/vcfanno, omit report generation (STEP 4)
   --panel_id VIRTUAL_PANEL_ID
                   Identifier for choice of predefined virtual cancer predisposition gene panels,
                   choose any between the following identifiers:
                   0 = CPSR exploratory cancer predisposition panel (n = 213, TCGA + Cancer Gene Census + NCGC)
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
                   19 = Inherited non-medullary thyroid cancer (Genomics England PanelApp)
                   20 = Inherited ovarian cancer (without breast cancer) (Genomics England PanelApp)
                   21 = Inherited pancreatic cancer (Genomics England PanelApp)
                   22 = Inherited renal cancer (Genomics England PanelApp)
                   23 = Inherited phaeochromocytoma and paraganglioma (Genomics England PanelApp)
                   24 = Melanoma pertinent cancer susceptibility (Genomics England PanelApp)
                   25 = Multiple endocrine tumours (Genomics England PanelApp)
                   26 = Multiple monogenic benign skin tumours (Genomics England PanelApp)
                   27 = Neuroendocrine cancer pertinent cancer susceptibility (Genomics England PanelApp)
                   28 = Neurofibromatosis Type 1 (Genomics England PanelApp)
                   29 = Ovarian cancer pertinent cancer susceptibility (Genomics England PanelApp)
                   30 = Parathyroid Cancer (Genomics England PanelApp)
                   31 = Prostate cancer pertinent cancer susceptibility (Genomics England PanelApp)
                   32 = Renal cancer pertinent cancer susceptibility (Genomics England PanelApp)
                   33 = Rhabdoid tumour predisposition (Genomics England PanelApp)
                   34 = Sarcoma cancer susceptibility (Genomics England PanelApp)
                   35 = Sarcoma susceptibility (Genomics England PanelApp)
                   36 = Thyroid cancer pertinent cancer susceptibility (Genomics England PanelApp)
                   37 = Tumour predisposition - childhood onset (Genomics England PanelApp)
                   38 = Upper gastrointestinal cancer pertinent cancer susceptibility (Genomics England PanelApp)

   --custom_panel TARGET_BED
                   Define custom screening panel through a four-column BED file (chrom,start,stop,genesymbol). Alternative to predefined panels provided with --panel_id
   --no_vcf_validate     Skip validation of input VCF with Ensembl's vcf-validator
   --diagnostic_grade_only
                   For panel_id's 1-38 (Genomics England PanelApp) - consider genes with a GREEN status only
   --docker-uid DOCKER_USER_ID
                   Docker user ID. Default is the host system user ID. If you are experiencing permission errors,
                   try setting this up to root (`--docker-uid root`)
   --no-docker           Run the CPSR workflow in a non-Docker mode
   --debug               Print full docker commands to log

The *cpsr* software bundle contains an example VCF file. It also comes
with a basic configuration file (*cpsr.toml*).

Report generation with the example VCF, using the `Adult solid tumours
cancer
susceptibility <https://panelapp.genomicsengland.co.uk/panels/245/>`__
virtual gene panel, can be performed through the following command:

``python ~/cpsr-0.5.2/cpsr.py ~/cpsr-0.5.2/example.vcf.gz``
``~/pcgr-0.8.4 ~/cpsr-0.5.2 grch37 --panel_id 1 ~/cpsr-0.5.2/cpsr.toml example``

Note that the example command also refers to the PCGR directory
(*pcgr-0.8.4*), which contains the data bundle that are necessary for
both *PCGR* and *CPSR*.

The command above will run the Docker-based *cpsr* workflow and produce
the following output files in the *cpsr* folder:

1. **example.cpsr.grch37.pass.vcf.gz (.tbi)** - Bgzipped VCF file with
   functional/clinical annotations
2. **example.cpsr.grch37.pass.tsv.gz** - Compressed TSV file (generated
   with `vcf2tsv <https://github.com/sigven/vcf2tsv>`__) with
   functional/clinical annotations
3. **example.cpsr.grch37.html** - Interactive HTML report with variants
   in cancer predisposition genes classified into five levels of
   pathogenicity
4. **example.cpsr.grch37.json.gz** - Compressed JSON dump of HTML report
   content
5. **example.cpsr.snvs_indels.tiers.grch37.tsv** - TSV file with most
   important annotations of tier-structured SNVs/InDels
