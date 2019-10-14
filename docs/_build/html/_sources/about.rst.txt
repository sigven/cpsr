About
-----

What is the Cancer Predisposition Sequencing Reporter (CPSR)?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The *Cancer Predisposition Sequencing Reporter (CPSR)* is a
computational workflow that **interprets germline variants** identified
from next-generation sequencing **in the context of cancer
predisposition**. The workflow is integrated with the framework that
underlies the `Personal Cancer Genome Reporter
(PCGR) <https://github.com/sigven>`__, utilizing the Docker environment
for encapsulation of code and software dependencies. While *PCGR* is
intended for reporting and analysis of somatic variants detected in a
tumor, *CPSR* is intended for reporting and ranking of germline variants
in protein-coding genes that are implicated in cancer predisposition and
inherited cancer syndromes.

*CPSR* accepts a query file with raw germline variant calls encoded in
the `VCF <https://samtools.github.io/hts-specs/VCFv4.2.pdf>`__ format
(i.e. analyzing SNVs/InDels). Furthermore, through the use several
different *virtual cancer predisposition gene panels* harvested from the
`Genomics England PanelApp <https://panelapp.genomicsengland.co.uk/>`__,
the user can flexibly put a restriction on which genes and findings are
displayed in the cancer predisposition report.

The software performs extensive variant annotation on the selected
geneset and produces an interactive HTML report, in which the user can
investigate four types of variants:

1. **ClinVar variants** - pre-classified variants according to a
   five-level tier scheme (Pathogenic to Benign)
2. **Non-ClinVar variants** - classified by CPSR according to a
   five-level tier scheme (Pathogenic to Benign)
3. **Secondary findings (optional)** - pathogenic ClinVar variants in
   the ACMG recommended list for reporting of incidental findings
4. **GWAS hits (optional)** - variants overlapping with previously
   identified hits in genome-wide association studies (GWAS) of cancer
   phenotypes (i.e. low to moderate risk conferring alleles), using
   `NHGRI-EBI Catalog of published genome-wide association
   studies <https://www.ebi.ac.uk/gwas/>`__ as the underlying source.

The variant sets can be interactively explored and filtered further
through different types of filters (phenotypes, genes, variant
consequences, population MAF etc.). Importantly, the unclassified
non-ClinVar variants are assigned a *pathogenicity level* based on the
aggregation of scores according to previously established `ACMG
criteria <https://www.ncbi.nlm.nih.gov/pubmed/25741868>`__. The ACMG
criteria includes cancer-specific criteria, as outlined and specified in
several previous studies (`Huang et al., Cell,
2018 <https://www.ncbi.nlm.nih.gov/pubmed/29625052>`__; `Nykamp et al.,
Genet Med., 2017 <https://www.ncbi.nlm.nih.gov/pubmed/28492532>`__;
`Maxwell et al., Am J Hum Genet.,
2016 <https://www.ncbi.nlm.nih.gov/pubmed/27153395>`__; `Amendola et
al., Am J Hum Genet.,
2016 <https://www.ncbi.nlm.nih.gov/pubmed/27181684>`__). See also
`Related work <https://github.com/sigven/cpsr#related-work>`__

Cancer predisposition genes
'''''''''''''''''''''''''''

The cancer predisposition report can show variants found in a number of
well-known cancer predisposition genes, and the specific set of genes
can be customized by the user by choosing any of the following **virtual
gene panels (0 - 36)**:

-  **Panel 0 (default)** is a comprehensive gene panel assembled through
   known sources on cancer predisposition:

   -  A list of 152 genes that were curated and established within
      TCGA’s pan-cancer study (`Huang et al., Cell,
      2018 <https://www.ncbi.nlm.nih.gov/pubmed/29625052>`__)
   -  A list of 107 protein-coding genes that has been manually curated
      in COSMIC’s `Cancer Gene Census
      v90 <https://cancer.sanger.ac.uk/census>`__,
   -  A list of 152 protein-coding genes established by experts within
      the Norwegian Cancer Genomics Consortium
      (http://cancergenomics.no)

   The combination of the three sources resulted in a non-redundant set
   of `213 protein-coding
   genes <https://github.com/sigven/cpsr/blob/master/predisposition.md>`__
   of relevance for predisposition to tumor development.

   -  **Panels 1 - 36** are panels for inherited cancer syndromes and
      cancer predisposition assembled within the `Genomics England
      PanelApp <https://panelapp.genomicsengland.co.uk/>`__:

      -  `1 = Adult solid tumours cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/245/>`__
      -  `2 = Adult solid tumours for rare
         disease <https://panelapp.genomicsengland.co.uk/panels/391/>`__
      -  `3 = Bladder cancer pertinent cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/208/>`__
      -  `4 = Brain cancer pertinent cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/166/>`__
      -  `5 = Breast cancer pertinent cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/55/>`__
      -  `6 = Childhood solid tumours cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/259/>`__
      -  `7 = Colorectal cancer pertinent cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/244/>`__
      -  `8 = Endometrial cancer pertinent cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/271/>`__
      -  `9 = Familial Tumours Syndromes of the central & peripheral
         Nervous
         system <https://panelapp.genomicsengland.co.uk/panels/167/>`__
      -  `10 = Familial breast
         cancer <https://panelapp.genomicsengland.co.uk/panels/158/>`__
      -  `11 = Familial
         melanoma <https://panelapp.genomicsengland.co.uk/panels/522/>`__
      -  `12 = Familial prostate
         cancer <https://panelapp.genomicsengland.co.uk/panels/318/>`__
      -  `13 = Familial
         rhabdomyosarcoma <https://panelapp.genomicsengland.co.uk/panels/290/>`__
      -  `14 = GI tract
         tumours <https://panelapp.genomicsengland.co.uk/panels/254/>`__
      -  `15 = Haematological malignancies cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/59/>`__
      -  `16 = Head and neck cancer pertinent cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/115/>`__
      -  `17 = Inherited non-medullary thyroid
         cancer <https://panelapp.genomicsengland.co.uk/panels/171/>`__
      -  `18 = Inherited ovarian cancer (without breast
         cancer) <https://panelapp.genomicsengland.co.uk/panels/143/>`__
      -  `19 = Inherited pancreatic
         cancer <https://panelapp.genomicsengland.co.uk/panels/524/>`__
      -  `20 = Inherited renal
         cancer <https://panelapp.genomicsengland.co.uk/panels/521/>`__
      -  `21 = Inherited phaeochromocytoma and
         paraganglioma <https://panelapp.genomicsengland.co.uk/panels/97/>`__
      -  `22 = Melanoma pertinent cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/133/>`__
      -  `23 = Multiple endocrine
         tumours <https://panelapp.genomicsengland.co.uk/panels/36/>`__
      -  `24 = Multiple monogenic benign skin
         tumours <https://panelapp.genomicsengland.co.uk/panels/558/>`__
      -  `25 = Neuroendocrine cancer pertinent cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/183/>`__
      -  `26 - Neurofibromatosis Type
         1 <https://panelapp.genomicsengland.co.uk/panels/255/>`__
      -  `27 = Ovarian cancer pertinent cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/117/>`__
      -  `28 = Parathyroid
         Cancer <https://panelapp.genomicsengland.co.uk/panels/86/>`__
      -  `29 = Prostate cancer pertinent cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/17/>`__
      -  `30 = Renal cancer pertinent cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/154/>`__
      -  `31 = Rhabdoid tumour
         predisposition <https://panelapp.genomicsengland.co.uk/panels/600/>`__
      -  `32 = Sarcoma cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/217/>`__
      -  `33 = Sarcoma
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/734/>`__
      -  `34 = Thyroid cancer pertinent cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/421/>`__
      -  `35 = Tumour predisposition - childhood
         onset <https://panelapp.genomicsengland.co.uk/panels/243/>`__
      -  `36 = Upper gastrointestinal cancer pertinent cancer
         susceptibility <https://panelapp.genomicsengland.co.uk/panels/273/>`__

Example report
~~~~~~~~~~~~~~

-  `Cancer predisposition sequencing
   report <http://folk.uio.no/sigven/example.cpsr.grch37.html>`__

Docker-based technology
~~~~~~~~~~~~~~~~~~~~~~~

The CPSR workflow is developed using the `Docker
technology <https://www.docker.com/what-docker>`__. The software is thus
packaged into isolated containers, in which the installation of all
software libraries/tools and required dependencies have been taken care
of. In addition to the bundled software, in the form of a Docker image,
the workflow needs to be attached with an `annotation data
bundle <annotation_resources.html>`__.

|image0|

Contact
~~~~~~~

sigven@ifi.uio.no

.. |image0| image:: docker-logo50.png

