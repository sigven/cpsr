CHANGELOG
---------

0.5.2 - November 18th 2019
^^^^^^^^^^^^^^^^^^^^^^^^^^

Changed
'''''''

-  Definition of pathogenic range (wrt to variant frequency) takes into
   account population- and position-specific allele numbers - no longer
   defined only by allele counts (i.e. *AC*) but by *AC* and *AN*
-  Moved virtual panel identifier from positional argument to optional
   argument (``--panel_id``) in ``cpsr.py``

Added
'''''

-  Ability to analyze custom panels, provided through option
   ``--custom_panel``. Needs to be defined as BED file with four
   columns, i.e. chromosome, start, stop, genesymbol

0.5.1 - October 14th 2019
^^^^^^^^^^^^^^^^^^^^^^^^^

Fixed
'''''

-  Bug in ``cpsr_validate_input.py``, `GitHub
   Issue <https://github.com/sigven/cpsr/issues/18>`__
-  Bug when there are zero variants with a ‘PASS’ status in VCF -
   omitting report generation

0.5.0 - September 23rd 2019
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _fixed-1:

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

.. _added-1:

Added
'''''

-  Section on *genomic biomarkers*; indicating which variants in the
   query VCF that overlaps with existing germline biomarkers (CIViC)
