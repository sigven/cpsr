# Running

  

The CPSR software comes with a range of options to configure the report
and analysis.

## Key settings

Below, we will outline some key settings that are important for the
generation of variant reports with CPSR.

### Virtual gene panel

The user can flexibly choose the set of cancer predisposition genes for
which input variants should be subject to classification and reporting.
The input VCF with germline variants may thus come from any sequencing
assay (whole-genome, whole-exome or targeted assay), yet the results
shown in the report will always be restricted by the choice of a
[virtual gene
panel](https://sigven.github.io/cpsr/dev/articles/virtual_panels.md):

- `--panel_id <panel_number>`

The user can choose from a range of pre-defined gene panels, selected
from the following list of panel identifiers:

- **0**: [Exploratory panel - all cancer predisposition
  genes](https://sigven.github.io/cpsr/articles/virtual_panels.html)
- **1**: [Adult solid tumours cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/245/)
- **2**: [Adult solid tumours for rare
  disease](https://panelapp.genomicsengland.co.uk/panels/391/)
- **3**: [Bladder cancer pertinent cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/208/)
- **4**: [Brain cancer pertinent cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/166/)
- **5**: [Breast cancer pertinent cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/55/)
- **6**: [Childhood solid tumours cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/259/)
- **7**: [Colorectal cancer pertinent cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/244/)
- **8**: [Endometrial cancer pertinent cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/271/)
- **9**: [Familial Tumours Syndromes of the central & peripheral Nervous
  system](https://panelapp.genomicsengland.co.uk/panels/167/)
- **10**: [Familial breast
  cancer](https://panelapp.genomicsengland.co.uk/panels/158/)
- **11**: [Familial
  melanoma](https://panelapp.genomicsengland.co.uk/panels/522/)
- **12**: [Familial prostate
  cancer](https://panelapp.genomicsengland.co.uk/panels/318/)
- **13**: [Familial
  rhabdomyosarcoma](https://panelapp.genomicsengland.co.uk/panels/290/)
- **14**: [GI tract
  tumours](https://panelapp.genomicsengland.co.uk/panels/254/)
- **15**: [Genodermatoses with
  malignancies](https://panelapp.genomicsengland.co.uk/panels/201/)
- **16**: [Haematological malignancies cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/59/)
- **17**: [Haematological malignancies for rare
  disease](https://panelapp.genomicsengland.co.uk/panels/407/)
- **18**: [Head and neck cancer pertinent cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/115/)
- **19**: [Inherited MMR deficiency (Lynch
  syndrome)](https://panelapp.genomicsengland.co.uk/panels/503/)
- **20**: [Inherited non-medullary thyroid
  cancer](https://panelapp.genomicsengland.co.uk/panels/171/)
- **21**: [Inherited ovarian cancer (without breast
  cancer)](https://panelapp.genomicsengland.co.uk/panels/143/)
- **22**: [Inherited pancreatic
  cancer](https://panelapp.genomicsengland.co.uk/panels/524/)
- **23**: [Inherited
  polyposis](https://panelapp.genomicsengland.co.uk/panels/504/)
- **24**: [Inherited predisposition to acute myeloid leukaemia
  (AML)](https://panelapp.genomicsengland.co.uk/panels/525/)
- **25**: [Inherited susceptibility to acute lymphoblastoid leukaemia
  (ALL)](https://panelapp.genomicsengland.co.uk/panels/1349/)
- **26**: [Inherited predisposition to
  GIST](https://panelapp.genomicsengland.co.uk/panels/523/)
- **27**: [Inherited renal
  cancer](https://panelapp.genomicsengland.co.uk/panels/521/)
- **28**: [Inherited phaeochromocytoma and
  paraganglioma](https://panelapp.genomicsengland.co.uk/panels/97/)
- **29**: [Melanoma pertinent cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/133/)
- **30**: [Multiple endocrine
  tumours](https://panelapp.genomicsengland.co.uk/panels/36/)
- **31**: [Multiple monogenic benign skin
  tumours](https://panelapp.genomicsengland.co.uk/panels/558/)
- **32**: [Neuroendocrine cancer pertinent cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/183/)
- **33**: [Neurofibromatosis Type
  1](https://panelapp.genomicsengland.co.uk/panels/255/)
- **34**: [Ovarian cancer pertinent cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/117/)
- **35**: [Parathyroid
  Cancer](https://panelapp.genomicsengland.co.uk/panels/86/)
- **36**: [Prostate cancer pertinent cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/17/)
- **37**: [Renal cancer pertinent cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/154/)
- **38**: [Sarcoma cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/217/)
- **39**: [Sarcoma
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/734/)
- **40**: [Thyroid cancer pertinent cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/421/)
- **41**: [Tumour predisposition - childhood
  onset](https://panelapp.genomicsengland.co.uk/panels/243/)
- **42**: [Upper gastrointestinal cancer pertinent cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/273/)
- **43**: [DNA repair genes pertinent cancer
  susceptibility](https://panelapp.genomicsengland.co.uk/panels/256/)

  

#### Custom-made virtual gene panels

CPSR allows users to create custom virtual gene panels for reporting.
Any set of genes found in the [CPSR superpanel (panel
0)](https://sigven.github.io/cpsr/dev/articles/virtual_panels.html#panel-0)
can be used to design a custom virtual gene panel. Technically, the
users need to create a simple one-column text (TSV) file with Ensembl
gene identifiers, and provide a name for the custom panel, using the
following command line options:

- `--custom_list <custom_list_tsv>`
- `--custom_list_name <custom_list_name`

### ACMG/AMP variant classification

There is an option to set an upper allele frequency limit (gnomAD global
population) for novel variants to be included in the report — a means to
exclude the potentially large number of common, benign variants in the
input VCF:

- `--max_af_gnomad <maf_threshold>`

------------------------------------------------------------------------

CPSR first computes its own ACMG/AMP rule-based classification for all
coding variants, then checks whether an existing ClinVar record should
override that result. The degree to which ClinVar takes precedence is
controlled by:

- `--clinvar_trust_level <0|1|2|3|4>`

| Level | Behaviour |
|----|----|
| **0** | ClinVar trusted — CPSR only overrides conflicted ClinVar records (default) |
| **1** | Override zero gold star ClinVar records |
| **2** | Override zero- and single gold star ClinVar records |
| **3** | Override low-star and non-cancer-phenotype ClinVar records |
| **4** | CPSR always classifies (ClinVar records never take precedence) |

------------------------------------------------------------------------

By default, CPSR does not report variants in cancer predisposition genes
that are associated only with non-cancer phenotypes in ClinVar. To
include these, use:

- `--clinvar_report_noncancer`

  

### Optional report contents

CPSR allows users to report recommended [incidental
findings](https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/), the
occurrence of important variants with respec to chemotherapy toxicity,
and also the genotypes of reported cancer risk loci from [genome-wide
association studies (GWAS)](https://www.ebi.ac.uk/gwas/):

- `--pgx_findings`
- `--secondary_findings`
- `--gwas_findings`

## All options

A cancer predisposition report is generated by running the **cpsr**
command, which takes the following arguments and options:

``` text
usage: cpsr -h [options]
    --input_vcf <INPUT_VCF>
    --vep_dir <VEP_DIR>
    --refdata_dir <REFDATA_DIR>
    --output_dir <OUTPUT_DIR>
    --genome_assembly <GENOME_ASSEMBLY>
    --sample_id <SAMPLE_ID>

Cancer Predisposition Sequencing Reporter - report of clinically significant cancer-predisposing germline variants

Required arguments:
  --input_vcf INPUT_VCF
                        VCF input file with small germline DNA variants (SNVs/InDels).
  --vep_dir VEP_DIR     Directory of VEP cache, e.g. $HOME/.vep
  --refdata_dir REFDATA_DIR
                        Directory that contains the PCGR/CPSR reference data, e.g. ~/pcgr-data-2.3.0
  --output_dir OUTPUT_DIR
                        Output directory
  --genome_assembly {grch37,grch38}
                        Genome assembly build: grch37 or grch38
  --sample_id SAMPLE_ID
                        Sample identifier - prefix for output files

Panel options:
  --panel_id VIRTUAL_PANEL_ID
                        Comma-separated string with identifier(s) of predefined virtual cancer predisposition gene panels,
                        choose any combination of the following identifiers (GEP = Genomics England PanelApp):
                        0 = CPSR exploratory cancer predisposition panel (PanelApp genes / TCGA's germline study / Cancer Gene Census / Other)
                        1 = Adult solid tumours cancer susceptibility (GEP)
                        2 = Adult solid tumours for rare disease (GEP)
                        3 = Bladder cancer pertinent cancer susceptibility (GEP)
                        4 = Brain cancer pertinent cancer susceptibility (GEP)
                        5 = Breast cancer pertinent cancer susceptibility (GEP)
                        6 = Childhood solid tumours cancer susceptibility (GEP)
                        7 = Colorectal cancer pertinent cancer susceptibility (GEP)
                        8 = Endometrial cancer pertinent cancer susceptibility (GEP)
                        9 = Familial Tumours Syndromes of the central & peripheral Nervous system (GEP)
                        10 = Familial breast cancer (GEP)
                        11 = Familial melanoma (GEP)
                        12 = Familial prostate cancer (GEP)
                        13 = Familial rhabdomyosarcoma (GEP)
                        14 = GI tract tumours (GEP)
                        15 = Genodermatoses with malignancies (GEP)
                        16 = Haematological malignancies cancer susceptibility (GEP)
                        17 = Haematological malignancies for rare disease (GEP)
                        18 = Head and neck cancer pertinent cancer susceptibility (GEP)
                        19 = Inherited MMR deficiency (Lynch syndrome) (GEP)
                        20 = Inherited non-medullary thyroid cancer (GEP)
                        21 = Inherited ovarian cancer (without breast cancer) (GEP)
                        22 = Inherited pancreatic cancer (GEP)
                        23 = Inherited polyposis (GEP)
                        24 = Inherited predisposition to acute myeloid leukaemia (AML) (GEP)
                        25 = Inherited susceptibility to acute lymphoblastoid leukaemia (ALL) (GEP)
                        26 = Inherited predisposition to GIST (GEP)
                        27 = Inherited renal cancer (GEP)
                        28 = Inherited phaeochromocytoma and paraganglioma (GEP)
                        29 = Melanoma pertinent cancer susceptibility (GEP)
                        30 = Multiple endocrine tumours (GEP)
                        31 = Multiple monogenic benign skin tumours (GEP)
                        32 = Neuroendocrine cancer pertinent cancer susceptibility (GEP)
                        33 = Neurofibromatosis Type 1 (GEP)
                        34 = Ovarian cancer pertinent cancer susceptibility (GEP)
                        35 = Parathyroid Cancer (GEP)
                        36 = Prostate cancer pertinent cancer susceptibility (GEP)
                        37 = Renal cancer pertinent cancer susceptibility (GEP)
                        38 = Sarcoma cancer susceptibility (GEP)
                        39 = Sarcoma susceptibility (GEP)
                        40 = Thyroid cancer pertinent cancer susceptibility (GEP)
                        41 = Tumour predisposition - childhood onset (GEP)
                        42 = Upper gastrointestinal cancer pertinent cancer susceptibility (GEP)
                        43 = DNA repair genes pertinent cancer susceptibility (GEP)
  --custom_list CUSTOM_LIST
                        Provide custom list of genes from virtual panel 0 (single-column .txt/.tsv file with Ensembl gene identifiers),
                        alternative to predefined panels provided with --panel_id
  --custom_list_name CUSTOM_LIST_NAME
                        Set name for custom made panel/list (single word - no whitespace), will be displayed in the report
  --diagnostic_grade_only
                        For panel_id's 1-43 (Genomics England PanelApp) - consider genes with a GREEN status only, default: False

Variant classification options:
  --secondary_findings  Include variants found in ACMG-recommended list for secondary findings (v3.3), default: False
  --pgx_findings        Report overlap with variants associated with chemotherapy toxicity (PgX findings, CPIC), default: False
  --gwas_findings       Report overlap with low to moderate cancer risk variants (tag SNPs) identified from genome-wide
                        association studies, default: False
  --max_af_gnomad MAX_AF_GNOMAD
                        Ignore novel variants with a gnomAD maximum allele frequency (AF) across populations greater than
                        this value, default: 0.9
  --clinvar_trust_level {0,1,2,3,4}
                        Level of trust assigned to ClinVar records relative to CPSR's own ACMG/AMP classification:
                        0 = ClinVar trusted (override conflicted records only) (default)
                        1 = Override zero gold star ClinVar records
                        2 = Override zero- and single gold star ClinVar records
                        3 = Override low-star and non-cancer-phenotype ClinVar records
                        4 = CPSR always classifies (ClinVar records never take precedence)
  --clinvar_report_noncancer
                        Report also ClinVar-classified variants attributed to phenotypes/conditions NOT directly related
                        to tumor development, default: False

VEP options:
  --vep_n_forks VEP_N_FORKS
                        Number of forks (option '--fork' in VEP), default: 4
  --vep_buffer_size VEP_BUFFER_SIZE
                        Variant buffer size (variants read into memory simultaneously, option '--buffer_size' in VEP)
                        - set lower to reduce memory usage, default: 500
  --vep_gencode_basic   Consider basic GENCODE transcript set only with Variant Effect Predictor (VEP)
                        (option '--gencode_basic' in VEP).
  --vep_pick_order VEP_PICK_ORDER
                        Comma-separated string of ordered transcript properties for primary variant pick
                        (option '--pick_order' in VEP),
                        default: mane_select,mane_plus_clinical,canonical,biotype,ccds,rank,tsl,appris,length
  --vep_no_intergenic   Skip intergenic variants during processing (option '--no_intergenic' in VEP), default: False

vcfanno options:
  --vcfanno_n_proc VCFANNO_N_PROC
                        Number of vcfanno processes (option '-p' in vcfanno), default: 4

Other options:
  --force_overwrite     By default, the script will fail with an error if any output file already exists.
                        You can force the overwrite of existing result files by using this flag, default: False
  --version             show program's version number and exit
  --no_reporting        Run functional variant annotation on VCF through VEP/vcfanno, omit classification/report
                        generation (STEP 4), default: False
  --no_html             Do not generate HTML report, default: False
  --retained_info_tags RETAINED_INFO_TAGS
                        Comma-separated string of VCF INFO tags from query VCF that should be kept in CPSR output TSV
  --ignore_noncoding    Ignore non-coding (i.e. non protein-altering) variants in report, default: False
  --debug               Print full commands to log
  --pcgrr_conda PCGRR_CONDA
                        pcgrr conda env name, default: pcgrr
```

## Example run

The *cpsr* R package comes with a test VCF file (calls from the GRCh37
human genome assembly) that can be used to test the CPSR workflow.
Please note that these are artificial germline calls, not originating
from an actual patient.

Report generation with the example VCF, using the [Adult solid tumours
cancer
susceptibility](https://panelapp.genomicsengland.co.uk/panels/245/) as
the virtual gene panel, can be performed through the following command:

``` bash
$ (base) conda activate pcgr
$ (pcgr)
cpsr \
     --input_vcf ~/cpsr-2.3.0/inst/examples/example.vcf.gz \
     --vep_dir ~/.vep \
     --refdata_dir ~/pcgr_ref_data \
     --output_dir ~/cpsr-2.3.0/ \
     --genome_assembly grch37 \
     --panel_id 1 \
     --sample_id example \
     --secondary_findings \
     --pgx_findings \
     --clinvar_trust_level 1 \
     --max_af_gnomad 0.2 \
     --force_overwrite
```

Note that the example command refers to the PCGR data bundle directory
(*refdata_dir*), which contains the data necessary for both *PCGR* and
*CPSR*.

This command will produce the following output files in the *output*
folder:

1.  **example.cpsr.grch37.vcf.gz (.tbi)** - Bgzipped VCF file with
    various variant annotations appended by CPSR
2.  **example.cpsr.grch37.pass.vcf.gz (.tbi)** - Bgzipped VCF file with
    various variant annotations appended by CPSR (PASS variants only)
3.  **example.cpsr.grch37.conf.yaml** - CPSR configuration file - output
    from pre-reporting annotation (Python) workflow
4.  **example.cpsr.grch37.pass.tsv.gz** - Compressed TSV file (generated
    with [vcf2tsvpy](https://github.com/sigven/vcf2tsvpy)) of VCF
    content with various annotations appended by CPSR
5.  **example.cpsr.grch37.xlsx** - An Excel workbook that contains
    - *i)* information on virtual gene panel interrogated for variants
    - *ii)* classification of clinical significance for variants
      overlapping with cancer predisposition genes
    - *iii)* secondary findings (if any found)
    - *iv)* match of variants with existing biomarkers (if any found)
    - *v)* overlap with pharmacogenomic variants (if any found)
6.  **example.cpsr.grch37.html** - Interactive HTML report with
    clinically relevant variants in cancer predisposition genes
7.  **example.cpsr.grch37.classification.tsv.gz** - TSV file with key
    annotations of germline SNVs/InDels classified according to clinical
    significance
8.  **example.cpsr.grch37.secondary_findings.tsv.gz** - TSV file with
    key annotations of variants found in ACMG-recommended list for
    secondary findings (v3.3)
9.  **example.cpsr.grch37.pgx_findings.tsv.gz** - TSV file with key
    annotations of variants associated with chemotherapy toxicity (PgX
    findings, CPIC)
