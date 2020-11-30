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

**IMPORTANT NOTE**: STEP 1 & 2 below outline installation guidelines for
running CPSR with Docker. If you want to install and run CPSR without
the use of Docker (i.e. through Conda), follow `these
instructions <https://github.com/sigven/cpsr/tree/master/conda_pkg/README.md>`__

STEP 1: Install PCGR (version 0.9.1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Make sure you have a working installation of PCGR (**version 0.9.1**)
and the accompanying data bundle(s) (walk through `steps
1-2 <https://github.com/sigven/pcgr#getting-started>`__).

STEP 2: Download the latest release
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download the `0.6.1
release <https://github.com/sigven/cpsr/releases/tag/v0.6.1>`__ of
*cpsr* (run script and configuration file)

STEP 3: Configuration
~~~~~~~~~~~~~~~~~~~~~

A few elements of the workflow can be figured using the *cpsr*
configuration file, encoded in
`TOML <https://github.com/toml-lang/toml>`__. The following can be
configured:

-  Choice of gnomAD control population
-  VEP/*vcfanno* options
-  Metadata regarding custom panel

See section on
`Input <https://cpsr.readthedocs.io/en/latest/input.html>`__ for more
details wrt. default configuration.

STEP 4: Run example
~~~~~~~~~~~~~~~~~~~

::

   usage: cpsr.py -h [options]  --query_vcf INPUT_VCF --pcgr_dir PCGR_DIR --output_dir OUTPUT_DIR --genome_assembly GENOME_ASSEMBLY --conf CONFIG_FILE --sample_id SAMPLE_ID

   Cancer Predisposition Sequencing Reporter - report of clinically significant cancer-predisposing germline variants

   Required arguments:
   --query_vcf QUERY_VCF
                   VCF input file with germline query variants (SNVs/InDels).
   --pcgr_dir PCGR_DIR   Directory that contains the PCGR data bundle directory, e.g. ~/pcgr-0.9.1
   --output_dir OUTPUT_DIR
                   Output directory
   --genome_assembly {grch37,grch38}
                   Genome assembly build: grch37 or grch38
   --conf CONFIGURATION_FILE
                   Configuration file in TOML format
   --sample_id SAMPLE_ID
                   Sample identifier - prefix for output files

   Optional arguments:
   --force_overwrite     By default, the script will fail with an error if any output file already exists.
                   You can force the overwrite of existing result files by using this flag, default: False
   --version             show program's version number and exit
   --basic               Run functional variant annotation on VCF through VEP/vcfanno, omit Tier assignment/report generation (STEP 4), default: False
   --panel_id VIRTUAL_PANEL_ID
               Identifier for choice of predefined virtual cancer predisposition gene panels,
               choose any between the following identifiers:
               0 = CPSR exploratory cancer predisposition panel (n = 216, TCGA + Cancer Gene Census + NCGC + Other)
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

   --custom_list CUSTOM_LIST   Provide custom list of genes from virtual panel 0 (single-column txt file with gene symbols), alternative to predefined panels provided with --panel_id)
   --no_vcf_validate     Skip validation of input VCF with Ensembl's vcf-validator, default: False
   --diagnostic_grade_only     For panel_id's 1-42 (Genomics England PanelApp) - consider genes with a GREEN status only, default: False
   --docker-uid DOCKER_USER_ID
       Docker user ID. Default is the host system user ID. If you are experiencing permission errors, try setting this up to root (`--docker-uid root`), default: None
   --no-docker           Run the CPSR workflow in a non-Docker mode, default: False
   --ignore_noncoding    Do not list non-coding variants in HTML report
   --secondary_findings    Include variants found in ACMG-recommended list for secondary/incidental findings (v2.0)
   --gwas_findings       Report overlap with low to moderate cancer risk variants (tag SNPs) identified from genome-wide association studies
   --classify_all        Provide CPSR variant classifications (TIER 1-5) also for variants with exising ClinVar classifications in output TSV
   --maf_upper_threshold MAF_UPPER_THRESHOLD
       Upper MAF limit (gnomAD global population frequency) for variants to be included in the report
   --debug               Print full docker commands to log, default: False

The *cpsr* software bundle contains an example VCF file. It also
contains a configuration file (*cpsr.toml*).

Report generation with the example VCF, using the `Adult solid tumours
cancer
susceptibility <https://panelapp.genomicsengland.co.uk/panels/245/>`__
virtual gene panel, can be performed through the following command:

::

   python ~/cpsr-0.6.1/cpsr.py
    --query_vcf ~/cpsr-0.6.1/example.vcf.gz
    --pcgr_dir ~/pcgr-0.9.1
    --output_dir ~/cpsr-0.6.1
    --genome_assembly grch37
    --panel_id 1
    --conf ~/cpsr-0.6.1/cpsr.toml
    --sample_id example
    --incidental_findings
    --classify_all
    --maf_upper_threshold 0.2
    --no_vcf_validate

Note that the example command also refers to the PCGR directory
(*pcgr-0.9.1*), which contains the data bundle that are necessary for
both *PCGR* and *CPSR*.

This command will run the Docker-based *cpsr* workflow and produce the
following output files in the *cpsr* folder:

1. **example.cpsr.grch37.pass.vcf.gz (.tbi)** - Bgzipped VCF file with
   relevant annotations appended by CPSR
2. **example.cpsr.grch37.pass.tsv.gz** - Compressed TSV file (generated
   with `vcf2tsv <https://github.com/sigven/vcf2tsv>`__) of VCF content
   with relevant annotations appended by CPSR
3. **example.cpsr.grch37.html** - Interactive HTML report with
   clinically relevant variants in cancer predisposition genes organized
   into tiers
4. **example.cpsr.grch37.json.gz** - Compressed JSON dump of HTML report
   content
5. **example.cpsr.grch37.snvs_indels.tiers.tsv** - TSV file with most
   important annotations of tier-structured SNVs/InDels
