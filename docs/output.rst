Output
------

Interactive HTML report
~~~~~~~~~~~~~~~~~~~~~~~

An interactive and tier-structured HTML report that lists variants in in
known cancer predisposition genes is provided with the following naming
convention:

<*sample_id*>.\ **cpsr**.<*genome_assembly*>.\ **html**

-  The **sample_id** is provided as input by the user, and reflects a
   unique identifier of the tumor-normal sample pair to be analyzed.

The report is structured in five main sections, described in more detail
below:

1. **Settings & annotation sources**

   -  Lists underlying tools and annotation sources (links and versions)
   -  Lists key configurations provided by user

2. **Introduction**

   -  Describes the structure and organization of the report, as well as
      the predisposition genes that have been included for reporting

3. **Summary of findings**

   -  Lists the top findings in Tier1/2/3

4. **Germline SNVs/InDels**

   -  Comprehensive, interactive variant tables for each tier
   -  Tier 1

      -  Pathogenic variants in ClinVar

   -  Tier 2

      -  Likely pathogenic variants in ClinVar

   -  Tier 3

      -  Variants of uncertain significance - color coded by *PATHRANK*
         - a categorical variable indicating the pathogenicity of each
         variant based on an accumulation of ACMG evidence levels

   -  GWAS

      -  Low-risk variants found in genome-wide association studies of
         cancer phenotypes (NHGRI-EBI Catalog)

5. **References**

   -  Supporting scientific literature (Interpretation/implementation of
      ACMG critera etc.)

Interactive datatables
^^^^^^^^^^^^^^^^^^^^^^

The interactive datatables contain a number of hyperlinked annotations
similar to those defined for the annotated VCF file, including the
following:

-  SYMBOL - Gene symbol (Entrez/NCBI)
-  PROTEIN_CHANGE - Amino acid change (VEP)
-  GENE_NAME - gene name description (Entrez/NCBI)
-  PROTEIN_DOMAIN - PFAM protein domain
-  PROTEIN_FEATURE - UniProt feature overlapping variant site
-  CDS_CHANGE - Coding sequence change
-  CONSEQUENCE - VEP consequence (primary transcript)
-  HGVSc - from VEP
-  HGVSp - from VEP
-  ONCOGENE - Known proto-oncogene
-  TUMOR_SUPPRESSOR - known tumor suppressor gene
-  ONCOSCORE - Literature-derived score for oncogenic potential (gene
   level)
-  PREDICTED_EFFECT - Effect predictions from dbNSFP
-  VEP_ALL_CONSEQUENCE - All VEP consequences (multiple transcripts)
-  DBSNP - dbSNP rsID
-  CLINVAR - ClinVar variant origin and associated phenotypes
-  GENOMIC_CHANGE - Variant ID
-  GENOME_VERSION - Genome assembly

Example report
~~~~~~~~~~~~~~

-  `Cancer predisposition sequencing
   report <http://folk.uio.no/sigven/example.cpsr.grch37.html>`__

The HTML reports have been tested using the following browsers:

-  Safari (10.0.3)
-  Mozilla Firefox (52.0.2)
-  Google Chrome (57.0.2987.110)

JSON
~~~~

A JSON file that stores the HTML report content is provided. This file
will easen the process of extracting particular parts of the report for
further analysis. Presently, there is no detailed schema documented for
the JSON structure.

Output files - germline SNVs/InDels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Variant call format - VCF
^^^^^^^^^^^^^^^^^^^^^^^^^

A VCF file containing annotated, germline calls (single nucleotide
variants and insertion/deletions) is generated with the following naming
convention:

<*sample_id*>.\ **cpsr**.<*genome_assembly*>.\ **vcf.gz (.tbi)**

Here, the **sample_id** is provided as input by the user, and reflects a
unique identifier of the tumor-normal sample pair to be analyzed.
Following common standards, the annotated VCF file is compressed with
`bgzip <http://www.htslib.org/doc/tabix.html>`__ and indexed with
`tabix <http://www.htslib.org/doc/tabix.html>`__. Below follows a
description of all annotations/tags present in the VCF INFO column after
processing with the CPSR annotation pipeline:

*VEP consequence annotations*
'''''''''''''''''''''''''''''

-  CSQ - Complete consequence annotations from VEP. Format:
   Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON\|
   INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids\|
   Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS\|
   SYMBOL_SOURCE|HGNC_ID|CANONICAL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC\|
   RefSeq|DOMAINS|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomAD_AF\|
   gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF\|
   gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO\|
   MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE
-  Consequence - Impact modifier for the consequence type (picked by
   VEP’s –flag_pick_allele option)
-  Gene - Ensembl stable ID of affected gene (picked by VEP’s
   –flag_pick_allele option)
-  Feature_type - Type of feature. Currently one of Transcript,
   RegulatoryFeature, MotifFeature (picked by VEP’s –flag_pick_allele
   option)
-  Feature - Ensembl stable ID of feature (picked by VEP’s
   –flag_pick_allele option)
-  cDNA_position - Relative position of base pair in cDNA sequence
   (picked by VEP’s –flag_pick_allele option)
-  CDS_position - Relative position of base pair in coding sequence
   (picked by VEP’s –flag_pick_allele option)
-  CDS_CHANGE - Coding, transcript-specific sequence annotation (picked
   by VEP’s –flag_pick_allele option)
-  AMINO_ACID_START - Protein position indicating absolute start of
   amino acid altered (fetched from Protein_position)
-  AMINO_ACID_END - Protein position indicating absolute end of amino
   acid altered (fetched from Protein_position)
-  Protein_position - Relative position of amino acid in protein (picked
   by VEP’s –flag_pick_allele option)
-  Amino_acids - Only given if the variant affects the protein-coding
   sequence (picked by VEP’s –flag_pick_allele option)
-  Codons - The alternative codons with the variant base in upper case
   (picked by VEP’s –flag_pick_allele option)
-  IMPACT - Impact modifier for the consequence type (picked by VEP’s
   –flag_pick_allele option)
-  VARIANT_CLASS - Sequence Ontology variant class (picked by VEP’s
   –flag_pick_allele option)
-  SYMBOL - Gene symbol (picked by VEP’s –flag_pick_allele option)
-  SYMBOL_SOURCE - The source of the gene symbol (picked by VEP’s
   –flag_pick_allele option)
-  STRAND - The DNA strand (1 or -1) on which the transcript/feature
   lies (picked by VEP’s –flag_pick_allele option)
-  ENSP - The Ensembl protein identifier of the affected transcript
   (picked by VEP’s –flag_pick_allele option)
-  FLAGS - Transcript quality flags: cds_start_NF: CDS 5’, incomplete
   cds_end_NF: CDS 3’ incomplete (picked by VEP’s –flag_pick_allele
   option)
-  SWISSPROT - Best match UniProtKB/Swiss-Prot accession of protein
   product (picked by VEP’s –flag_pick_allele option)
-  TREMBL - Best match UniProtKB/TrEMBL accession of protein product
   (picked by VEP’s –flag_pick_allele option)
-  UNIPARC - Best match UniParc accession of protein product (picked by
   VEP’s –flag_pick_allele option)
-  HGVSc - The HGVS coding sequence name (picked by VEP’s
   –flag_pick_allele option)
-  HGVSp - The HGVS protein sequence name (picked by VEP’s
   –flag_pick_allele option)
-  HGVSp_short - The HGVS protein sequence name, short version (picked
   by VEP’s –flag_pick_allele option)
-  HGVS_OFFSET - Indicates by how many bases the HGVS notations for this
   variant have been shifted (picked by VEP’s –flag_pick_allele option)
-  MOTIF_NAME - The source and identifier of a transcription factor
   binding profile aligned at this position (picked by VEP’s
   –flag_pick_allele option)
-  MOTIF_POS - The relative position of the variation in the aligned
   TFBP (picked by VEP’s –flag_pick_allele option)
-  HIGH_INF_POS - A flag indicating if the variant falls in a high
   information position of a transcription factor binding profile (TFBP)
   (picked by VEP’s –flag_pick_allele option)
-  MOTIF_SCORE_CHANGE - The difference in motif score of the reference
   and variant sequences for the TFBP (picked by VEP’s –flag_pick_allele
   option)
-  CELL_TYPE - List of cell types and classifications for regulatory
   feature (picked by VEP’s –flag_pick_allele option)
-  CANONICAL - A flag indicating if the transcript is denoted as the
   canonical transcript for this gene (picked by VEP’s –flag_pick_allele
   option)
-  CCDS - The CCDS identifier for this transcript, where applicable
   (picked by VEP’s –flag_pick_allele option)
-  INTRON - The intron number (out of total number) (picked by VEP’s
   –flag_pick_allele option)
-  EXON - The exon number (out of total number) (picked by VEP’s
   –flag_pick_allele option)
-  DOMAINS - The source and identifier of any overlapping protein
   domains (picked by VEP’s –flag_pick_allele option)
-  DISTANCE - Shortest distance from variant to transcript (picked by
   VEP’s –flag_pick_allele option)
-  BIOTYPE - Biotype of transcript or regulatory feature (picked by
   VEP’s –flag_pick_allele option)
-  TSL - Transcript support level (picked by VEP’s –flag_pick_allele
   option)>
-  PUBMED - PubMed ID(s) of publications that cite existing variant -
   VEP
-  PHENO - Indicates if existing variant is associated with a phenotype,
   disease or trait - VEP
-  GENE_PHENO - Indicates if overlapped gene is associated with a
   phenotype, disease or trait - VEP
-  ALLELE_NUM - Allele number from input; 0 is reference, 1 is first
   alternate etc - VEP
-  REFSEQ_MATCH - The RefSeq transcript match status; contains a number
   of flags indicating whether this RefSeq transcript matches the
   underlying reference sequence and/or an Ensembl transcript (picked by
   VEP’s –flag_pick_allele option)
-  PICK - Indicates if this block of consequence data was picked by
   VEP’s –flag_pick_allele option
-  VEP_ALL_CONSEQUENCE - All transcript consequences
   (Consequence:SYMBOL:Feature_type:Feature:BIOTYPE) - VEP

*Gene information*
''''''''''''''''''

-  ENTREZ_ID - `Entrez <http://www.ncbi.nlm.nih.gov/gene>`__ gene
   identifier
-  APPRIS - Principal isoform flags according to the `APPRIS principal
   isoform database <http://appris.bioinfo.cnio.es/#/downloads>`__
-  UNIPROT_ID - `UniProt <http://www.uniprot.org>`__ identifier
-  UNIPROT_ACC - `UniProt <http://www.uniprot.org>`__ accession(s)
-  ENSEMBL_GENE_ID - Ensembl gene identifier for VEP’s picked transcript
   (*ENSGXXXXXXX*)
-  ENSEMBL_TRANSCRIPT_ID - Ensembl transcript identifier for VEP’s
   picked transcript (*ENSTXXXXXX*)
-  REFSEQ_MRNA - Corresponding RefSeq transcript(s) identifier for VEP’s
   picked transcript (*NM_XXXXX*)
-  CORUM_ID - Associated protein complexes (identifiers) from
   `CORUM <http://mips.helmholtz-muenchen.de/corum/>`__
-  DISGENET_CUI - Tumor types associated with gene, as found in
   DisGeNET. Tumor types are listed as unique
   `MedGen <https://www.ncbi.nlm.nih.gov/medgen/>`__ concept IDs
   (*CUIs*)
-  TUMOR_SUPPRESSOR - Gene is predicted as tumor suppressor candidate
   according to
   (`CancerMine <https://zenodo.org/record/1336650#.W9do9WJKiL4>`__)
-  ONCOGENE - Gene is predicted as an oncogene according to
   (`CancerMine <https://zenodo.org/record/1336650#.W9do9WJKiL4>`__)
-  ONCOSCORE - Literature-derived score for cancer gene relevance
   `Bioconductor/OncoScore <http://bioconductor.org/packages/release/bioc/html/OncoScore.html>`__,
   range from 0 (low oncogenic potential) to 1 (high oncogenic
   potential)
-  CANCER_SUSCEPTIBILITY_CUI - MedGen concept unique identifier (CUI)
   for cancer phenotype
-  CANCER_SYNDROME_CUI - MedGen concept unique identifier (CUI) for
   cancer syndrome
-  CANCER_PREDISPOSITION_SOURCE - Data source for susceptibility gene
   (NCGC, CGC, TCGA_PANCAN)
-  CANCER_PREDISPOSITION_MOI - Mechanism of inheritance for
   susceptibility gene (AR/AD)
-  CANCER_PREDISPOSITION_MOD - Mechanism of disease for susceptibility
   gene (Lof/GoF)

*Variant effect and protein-coding information*
'''''''''''''''''''''''''''''''''''''''''''''''

-  MUTATION_HOTSPOT - mutation hotspot codon in
   `cancerhotspots.org <http://cancerhotspots.org/>`__. Format:
   gene_symbol \| codon \| q-value
-  MUTATION_HOTSPOT_TRANSCRIPT - hotspot-associated transcripts (Ensembl
   transcript ID)
-  MUTATION_HOTSPOT_CANCERTYPE - hotspot-associated cancer types (from
   cancerhotspots.org)
-  UNIPROT_FEATURE - Overlapping protein annotations from `UniProt
   KB <http://www.uniprot.org>`__
-  PFAM_DOMAIN - Pfam domain identifier (from VEP)
-  EFFECT_PREDICTIONS - All predictions of effect of variant on protein
   function and pre-mRNA splicing from `database of non-synonymous
   functional predictions - dbNSFP
   v3.5 <https://sites.google.com/site/jpopgen/dbNSFP>`__. Predicted
   effects are provided by different sources/algorithms (separated by
   ‘&’):

   1.  `SIFT <http://provean.jcvi.org/index.php>`__ (Jan 2015)
   2.  `LRT <http://www.genetics.wustl.edu/jflab/lrt_query.html>`__
       (2009)
   3.  `MutationTaster <http://www.mutationtaster.org/>`__ (data release
       Nov 2015)
   4.  `MutationAssessor <http://mutationassessor.org/>`__ (release 3)
   5.  `FATHMM <http://fathmm.biocompute.org.uk>`__ (v2.3)
   6.  `PROVEAN <http://provean.jcvi.org/index.php>`__ (v1.1 Jan 2015)
   7.  `FATHMM_MKL <http://fathmm.biocompute.org.uk/fathmmMKL.htm>`__
   8.  `DBNSFP_CONSENSUS_SVM <https://www.ncbi.nlm.nih.gov/pubmed/25552646>`__
       (Ensembl/consensus prediction, based on support vector machines)
   9.  `DBNSFP_CONSENSUS_LR <https://www.ncbi.nlm.nih.gov/pubmed/25552646>`__
       (Ensembl/consensus prediction, logistic regression based)
   10. `SPLICE_SITE_EFFECT_ADA <http://nar.oxfordjournals.org/content/42/22/13534>`__
       (Ensembl/consensus prediction of splice-altering SNVs, based on
       adaptive boosting)
   11. `SPLICE_SITE_EFFECT_RF <http://nar.oxfordjournals.org/content/42/22/13534>`__
       (Ensembl/consensus prediction of splice-altering SNVs, based on
       random forest)
   12. `M-CAP <http://bejerano.stanford.edu/MCAP>`__
   13. `MutPred <http://mutpred.mutdb.org>`__
   14. `GERP <http://mendel.stanford.edu/SidowLab/downloads/gerp/>`__

-  SIFT_DBNSFP - predicted effect from SIFT (dbNSFP)
-  PROVEAN_DBNSFP - predicted effect from PROVEAN (dbNSFP)
-  MUTATIONTASTER_DBNSFP - predicted effect from MUTATIONTASTER (dbNSFP)
-  MUTATIONASSESSOR_DBNSFP - predicted effect from MUTATIONASSESSOR
   (dbNSFP)
-  M_CAP_DBNSFP - predicted effect from M-CAP (dbNSFP)
-  MUTPRED_DBNSFP - score from MUTPRED (dbNSFP)
-  FATHMM_DBNSFP - predicted effect from FATHMM (dbNSFP)
-  FATHMM_MKL_DBNSFP - predicted effect from FATHMM-mkl (dbNSFP)
-  META_LR_DBNSFP - predicted effect from ensemble prediction (logistic
   regression - dbNSFP)
-  SPLICE_SITE_RF_DBNSFP - predicted effect of splice site disruption,
   using random forest (dbscSNV)
-  SPLICE_SITE_ADA_DBNSFP - predicted effect of splice site disruption,
   using boosting (dbscSNV)

*Variant frequencies/annotations in germline databases*
'''''''''''''''''''''''''''''''''''''''''''''''''''''''

-  AFR_AF_GNOMAD - African/American germline allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  AMR_AF_GNOMAD - American germline allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  GLOBAL_AF_GNOMAD - Adjusted global germline allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  SAS_AF_GNOMAD - South Asian germline allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  EAS_AF_GNOMAD - East Asian germline allele frequency (`Genome
   Aggregation Database release
   21 <http://gnomad.broadinstitute.org/>`__)
-  FIN_AF_GNOMAD - Finnish germline allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  NFE_AF_GNOMAD - Non-Finnish European germline allele frequency
   (`Genome Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  OTH_AF_GNOMAD - Other germline allele frequency (`Genome Aggregation
   Database release 2 <http://gnomad.broadinstitute.org/>`__)
-  ASJ_AF_GNOMAD - Ashkenazi Jewish allele frequency (`Genome
   Aggregation Database release
   2 <http://gnomad.broadinstitute.org/>`__)
-  AFR_AF_1KG - `1000G Project - phase 3 <http://www.1000genomes.org>`__
   germline allele frequency for samples from AFR (African)
-  AMR_AF_1KG - `1000G Project - phase 3 <http://www.1000genomes.org>`__
   germline allele frequency for samples from AMR (Ad Mixed American)
-  EAS_AF_1KG - `1000G Project - phase 3 <http://www.1000genomes.org>`__
   germline allele frequency for samples from EAS (East Asian)
-  EUR_AF_1KG - `1000G Project - phase 3 <http://www.1000genomes.org>`__
   germline allele frequency for samples from EUR (European)
-  SAS_AF_1KG - `1000G Project - phase 3 <http://www.1000genomes.org>`__
   germline allele frequency for samples from SAS (South Asian)
-  GLOBAL_AF_1KG - `1000G Project - phase
   3 <http://www.1000genomes.org>`__ germline allele frequency for all
   1000G project samples (global)
-  DBSNPRSID - `dbSNP <http://www.ncbi.nlm.nih.gov/SNP/>`__ reference
   ID, as provided by VEP

*Clinical associations*
'''''''''''''''''''''''

-  CLINVAR_MSID - `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__
   Measure Set/Variant ID
-  CLINVAR_ALLELE_ID - `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__
   allele ID
-  CLINVAR_PMID - Associated Pubmed IDs for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ - germline
   state-of-origin
-  CLINVAR_HGVSP - Protein variant expression using HGVS nomenclature
-  CLINVAR_PMID_SOMATIC - Associated Pubmed IDs for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ - somatic
   state-of-origin
-  CLINVAR_CLNSIG - Clinical significance for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ - germline
   state-of-origin
-  CLINVAR_CLNSIG_SOMATIC - Clinical significance for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ - somatic
   state-of-origin
-  CLINVAR_MEDGEN_CUI - Associated
   `MedGen <https://www.ncbi.nlm.nih.gov/medgen/>`__ concept identifiers
   (*CUIs*) - germline state-of-origin
-  CLINVAR_MEDGEN_CUI_SOMATIC - Associated
   `MedGen <https://www.ncbi.nlm.nih.gov/medgen/>`__ concept identifiers
   (*CUIs*) - somatic state-of-origin
-  CLINVAR_VARIANT_ORIGIN - Origin of variant (somatic, germline, de
   novo etc.) for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__
-  GWAS_HIT - variant associated with cancer phenotype from genome-wide
   association study (NHGRI-EBI GWAS catalog)

Tab-separated values (TSV)
^^^^^^^^^^^^^^^^^^^^^^^^^^

Annotated List of all SNVs/InDels
'''''''''''''''''''''''''''''''''

We provide a tab-separated values file with most important variant/gene
annotations. The file has the following naming convention:

<*sample_id*>.\ **cpsr**.<*genome_assembly*>.\ **snvs_indels.tiers.tsv**

The SNVs/InDels are organized into different **tiers** (as defined above
for the HTML report)

The following variables are included in the tiered TSV file:

::

   1. GENOMIC_CHANGE - Identifier for variant at the genome (VCF) level, e.g. 1:g.152382569A>G
         Format: (<chrom>:g.<position><ref_allele>><alt_allele>)
   2. VAR_ID - Variant identifier
   3. GENOTYPE - Variant genotype (heterozygous/homozygous)
   4. GENOME_VERSION - Assembly version, e.g. GRCh37
   5. VCF_SAMPLE_ID - Sample identifier
   6. VARIANT_CLASS - Variant type, e.g. SNV/insertion/deletion
   7. CODING_STATUS - coding/noncoding (wrt. protein alteration and canonical splice site disruption)
   8. SYMBOL - Gene symbol
   9. GENE_NAME - Gene description
   10. CCDS - CCDS identifier
   11. ENTREZ_ID - Entrez gene identifier
   12. UNIPROT_ID - UniProt protein identifier
   13. ENSEMBL_GENE_ID - Ensembl gene identifier
   14. ENSEMBL_TRANSCRIPT_ID - Ensembl transcript identifier
   15. REFSEQ_MRNA - RefSeq mRNA identifier
   16. ONCOSCORE - Literature-derived score for cancer gene relevance
   17. ONCOGENE - Gene is predicted as an oncogene according to literature mining (CancerMine)
   18. TUMOR_SUPPRESSOR - Gene is predicted as tumor suppressor according to literature mining (CancerMine)
   19. MOD - Mechanism of disease for cancer predisposition gene (Lof/GoF/NA)
   20. PATH_TRUNCATION_RATE - Rate of pathogenic truncations in predisposition gene
   21. BENIGN_MISSENSE_RATE - Rate of benign missense mutations in predisposition gene
   22. CONSEQUENCE - Variant consequence
   23. PROTEIN_CHANGE - Protein change - one letter abbreviation (HGVSp)
   24. PROTEIN_DOMAIN - Protein domain description from PFAM
   25. HGVSp - The HGVS protein sequence name
   26. HGVSc - The HGVS coding sequence name
   27. CDS_CHANGE - Coding, transcript-specific sequence annotation
   28. MUTATION_HOTSPOT - Cancer mutation hotspot (cancerhotspots.org)
   29. RMSK_HIT - RepeatMasker hit
   30. PROTEIN_FEATURE - Protein feature (active sites etc.) from UniProt KnowledgeBase
   31. EFFECT_PREDICTIONS - Functional effect predictions from multiple algorithms (dbNSFP)
   32. LOSS_OF_FUNCTION - Loss-of-function variant, as predicted from VEP's LofTee plugin
   33. DBSNP - dbSNP identifier
   34. CLINVAR_CLINICAL_SIGNIFICANCE - clinical significance of ClinVar variant
   35. CLINVAR_MSID - measureset identifier of ClinVar variant
   36. CLINVAR_VARIANT_ORIGIN - variant origin (somatic/germline) of ClinVar variant
   37. CLINVAR_CONFLICTED - indicator of conflicting interpretations
   38. CLINVAR_PHENOTYPE - associated phenotype(s) for ClinVar variant
   39. VEP_ALL_CONSEQUENCE - all consequences from VEP
   40. PATHSCORE - Composite variant pathogenicity score
   41. PATHRANK - Categorical ranking of pathogenicity (HIGH,MODERATE,LOW,BENIGN)
   42. PATHDOC - Documentation for PATHSCORE/PATHRANK
   43. PVS1 - ACMG criteria: Null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multiexon deletion) - in a gene where LoF is a known mechanism of disease (Dom MoI)
   44. PSC1 - ACMG criteria: Null variant (nonsense, frameshift, canonical ±1 or 2 splice sites, initiation codon, single or multiexon deletion) - in a gene where LoF is a known mechanism of disease (Rec MoI)
   45. PS1 - ACMG criteria: Same amino acid change as a previously established pathogenic variant regardless of nucleotide change
   46. PMC1 - ACMG criteria: Null variant in a gene where LoF is not a known mechanism of disease
   47. PM1 - ACMG criteria: Variant located in mutation hotspot (cancerhotspots.org)
   48. PM2 - ACMG criteria: Absence or extremely low frequency (MAF < 0.0005) in 1000 Genomes Project, or gnomAD
   49. PM4 - ACMG criteria: Protein length changes due to inframe indels or nonstop variant in non-repetitive regions of genes that harbor variants with a Dom MoI.
   50. PM5 - ACMG criteria: Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before
   51. PP2 - ACMG criteria: Missense variant in a gene that has a relatively low rate of benign missense variation (<20%) and where missense variants are a common mechanism of disease (>50% of pathogenic variants)
   52. PPC1 - ACMG criteria: Protein length changes due to inframe indels or nonstop variant in non-repetitive regions of genes that harbor variants with a Rec MoI.
   53. PP3 - ACMG criteria: Multiple lines of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact, etc. - from dbNSFP)
   54. BP4 - ACMG criteria: Multiple lines of computational evidence support a benign effect on the gene or gene product (conservation, evolutionary, splicing impact, etc. - from dbNSFP).
   55. BMC1 - ACMG criteria: Peptide change is at the same location of a known benign change
   56. BSC1 - ACMG criteria: Peptide change is known to be benign
   57. BA1 - ACMG criteria: Common population frequency (MAF > 0.05) in 1000 Genomes Project, or gnomAD
   58. BP1 - ACMG criteria: Missense variant in a gene for which primarily truncating variants (>90%) are known to cause disease
   59. N_INSILICO_CALLED - Number of algorithms with effect prediction (damaging/tolerated) from dbNSFP
   60. N_INSILICO_DAMAGING - Number of algorithms with damaging prediction from dbNSFP
   61. N_INSILICO_TOLERATED - Number of algorithms with tolerated prediction from dbNSFP
   62. N_INSILICO_SPLICING_NEUTRAL - Number of algorithms with splicing neutral prediction from dbscSNV
   63. N_INSILICO_SPLICING_AFFECTED - Number of algorithms with splicing affected prediction from dbscSNV
   64. GLOBAL_AF_GNOMAD - Global MAF in gnomAD
   65. <CUSTOM_POPULATION_GNOMAD> - Population specific MAF in gnomAD (population configured by user)
   66. GLOBAL_AF_1KG - Global MAF in 1000 Genomes Project
   67. <CUSTOM_POPULAION_1KG> - Population-specific MAF in 1000 Genomes Project (population configured by user)
   68. TIER - CPSR tier level
   69. TIER_DESCRIPTION - Description of CPSR tier

**NOTE**: The user has the possibility to append the TSV file with data
from other tags in the input VCF of interest (i.e. using the
*custom_tags* option in the TOML configuration file)
