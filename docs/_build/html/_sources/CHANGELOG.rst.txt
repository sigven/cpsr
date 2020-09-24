CHANGELOG
---------

0.6.0 - September 23rd 2020
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fixed
'''''

-  Duplicated entries in incidental findings

Changed
'''''''

-  All arguments to ``cpsr.py`` is now organized with a ``--`` (no
   positional arguments)
-  Arguments to ``cpsr.py`` are divided into two groups: *required* and
   *optional*
-  ``secondary_findings`` is now coined ``incidental_findings``
-  Option **gwas:gwas_hits** in CPSR configuration file is now option
   ``--gwas_findings`` in ``cpsr.py``
-  Option **classification:clinvar_cpsr** in CPSR configuration file is
   now option ``--classify_all`` in ``cpsr.py``
-  Option **maf_imits:maf_gnomad** in CPSR configuration file is now
   option ``--maf_upper_threshold`` in ``cpsr.py``
-  Option **secondary_findings:show_sf** in CPSR configuration file is
   now option ``--incidental_findings`` in ``cpsr.py``
-  Virtual panels is now displayed through HTML (previously static
   ggplot plot)
-  **Settings** section of report is now divived into three:

   -  Sample metadata
   -  Report configuration
   -  Virtual panel

Added
'''''

-  Missing ACMG criteria for classification of silent and intronic
   variants outside of splice regions (ACMG_BP7)
-  Missing ACMG criterion for classification of variants in promoter and
   untranslated regions (ACMG_BP3)
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

.. _changed-1:

Changed
'''''''

-  Definition of pathogenic range (wrt to variant frequency) takes into
   account population- and position-specific allele numbers - no longer
   defined only by allele counts (i.e. *AC*) but by *AC* and *AN*
-  Moved virtual panel identifier from positional argument to optional
   argument (``--panel_id``) in ``cpsr.py``

.. _added-1:

Added
'''''

-  Ability to analyze custom panels, provided through option
   ``--custom_panel``. Needs to be defined as BED file with four
   columns, i.e. chromosome, start, stop, genesymbol

0.5.1 - October 14th 2019
^^^^^^^^^^^^^^^^^^^^^^^^^

.. _fixed-1:

Fixed
'''''

-  Bug in ``cpsr_validate_input.py``, `GitHub
   Issue <https://github.com/sigven/cpsr/issues/18>`__
-  Bug when there are zero variants with a ‘PASS’ status in VCF -
   omitting report generation

0.5.0 - September 23rd 2019
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _fixed-2:

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

.. _added-2:

Added
'''''

-  Section on *genomic biomarkers*; indicating which variants in the
   query VCF that overlaps with existing germline biomarkers (CIViC)
