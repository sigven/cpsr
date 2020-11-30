CHANGELOG
---------

0.6.1 - November 30th 2020
^^^^^^^^^^^^^^^^^^^^^^^^^^

Added
'''''

-  Increased number of genes in panel 0: All genes in 42 virtual panels
   related to cancer conditions in Genomics England PanelApp now also
   contributes toward panel 0
-  Added option in main script (``--clinvar_ignore_noncancer``) that
   will exclude any query variants (from HTML report and TSV/JSON
   output) that have been reported and classified for non-cancer related
   conditions only (in ClinVar)

   -  this to exclude variants associated with non-cancer related
      phenotypes

-  For the variant biomarker table, the resolution of the reported
   biomarker mapping is highlighted with designated background colors
   for the gene (exact/codon - black vs. exon/gene - orange)

Fixed
'''''

-  Bug in GWAS hits retrieval, `Issue
   #30 <https://github.com/sigven/cpsr/issues/18>`__
-  Custom VCF tags (as specified by user in configuration file) not
   shown in output TSV files

Changed
'''''''

-  Removed DisGeNET annotations from output (associations from Open
   Targets Platform serve same purpose)
-  Renamed report section **Genomic Biomarkers** to **Variant
   Biomarkers**
-  Option ``--incidental_findings`` changed back to
   ``--secondary_findings`` - recommended term to use according to ACMG
-  Removed *MOD (mechanism-of-disease)* from TSV output file

0.6.0rc - September 24th 2020
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Data updates: ClinVar, GWAS catalog, GENCODE, CIViC, CancerMine,
   UniProt KB, dbNSFP, Pfam, KEGG, Open Targets Platform, Genomics
   England PanelApp
-  Software updates: VEP 101

.. _fixed-1:

Fixed
'''''

-  Duplicated entries in incidental findings

.. _changed-1:

Changed
'''''''

-  All arguments to ``cpsr.py`` are now non-positional
-  Arguments to ``cpsr.py`` are divided into two groups: *required* and
   *optional*
-  ``secondary_findings`` is now coined ``incidental_findings``
-  Option **gwas:gwas_hits** in CPSR configuration file is now optional
   argument ``--gwas_findings`` in ``cpsr.py``
-  Option **classification:clinvar_cpsr** in CPSR configuration file is
   now optional argument ``--classify_all`` in ``cpsr.py``
-  Option **maf_imits:maf_gnomad** in CPSR configuration file is now
   optional argument ``--maf_upper_threshold`` in ``cpsr.py``
-  Option **secondary_findings:show_sf** in CPSR configuration file is
   now optional argument ``--incidental_findings`` in ``cpsr.py``
-  Virtual panels is now displayed through HTML (previously static
   ggplot plot)
-  **Settings** section of report is now divived into three:

   -  Sample metadata
   -  Report configuration
   -  Virtual panel

-  Classifications of genes as tumor suppressors/oncogenes are now based
   on a combination of CancerMine citation count and presence in Network
   of Cancer Genes

.. _added-1:

Added
'''''

-  Missing ACMG criterion for classification of silent and intronic
   variants outside of splice regions (*ACMG_BP7*)
-  Missing ACMG criterion for classification of variants in promoter and
   untranslated regions (*ACMG_BP3*)
-  Possibility to create custom virtual panel - any combination of genes
   from panel 0 provided as a single-column text file with argument
   ``--custom_list``
-  Ensured that non-empty datatables pr. tier (**ClinVar** and
   **Non-ClinVar**) are set as the active tab
-  Improved documentation of variant classification in the
   **References** section
-  DOIs available for all references

0.5.2 - November 18th 2019
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _changed-2:

Changed
'''''''

-  Definition of pathogenic range (wrt to variant frequency) takes into
   account population- and position-specific allele numbers - no longer
   defined only by allele counts (i.e. *AC*) but by *AC* and *AN*
-  Moved virtual panel identifier from positional argument to optional
   argument (``--panel_id``) in ``cpsr.py``

.. _added-2:

Added
'''''

-  Ability to analyze custom panels, provided through option
   ``--custom_panel``. Needs to be defined as BED file with four
   columns, i.e. chromosome, start, stop, genesymbol

0.5.1 - October 14th 2019
^^^^^^^^^^^^^^^^^^^^^^^^^

.. _fixed-2:

Fixed
'''''

-  Bug in ``cpsr_validate_input.py``, `GitHub
   Issue <https://github.com/sigven/cpsr/issues/18>`__
-  Bug when there are zero variants with a ‘PASS’ status in VCF -
   omitting report generation

0.5.0 - September 23rd 2019
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _fixed-3:

Fixed
'''''

-  Bug in implementation of ACMG criteria; genes without a known
   loss-of-function mechanism were handled inappropriately
-  Bug in assignment of heterozygous/homozygous states (input VCF)
-  Bug in implementation of ACMG_PS1 - Same amino acid change as
   previously pathogenic variant
-  Improved consequence prioritisation for variants with transcript
   consequences in multiple, distinct cancer predisposition genes
-  Upper MAF threshold (as given by user) only applied for unclassified
   (i.e. non-ClinVar variants)
-  Handling of non-coding variants (synonymous, upstream_variants) in
   the report, no longer excluded

.. _added-3:

Added
'''''

-  Section on *genomic biomarkers*; indicating which variants in the
   query VCF that overlaps with existing germline biomarkers (CIViC)
