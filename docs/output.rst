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

1. **Settings**

   -  Lists key configurations provided by user, including the list of
      genes that constitute the virtual gene panel in the report

2. **Summary of findings**

   -  Summarizes the findings through donut charts

      -  Number of variants in each of the five variant classification
         levels

3. **Variant classification**

   -  For all coding variants in the selected cancer predisposition
      geneset, interactive variant tables are shown for each level
      (**ClinVar** and **non-ClinVar (Other)** variants combined):

      -  Pathogenic
      -  Likely Pathogenic
      -  Variants of Uncertain Significance (VUS)
      -  Likely Benign
      -  Benign

   -  Biomarkers

      -  Reported clinical evidence items from
         `CIViC <https://civicdb.org>`__ that overlap with variants in
         the query set are reported in four distinct tabs (Predictive /
         Prognostic / Diagnostic / Predisposing)

   -  Secondary findings

      -  Pathogenic variants in the `ACMG recommended list of genes for
         report of secondary/incidental
         findings <https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/>`__

   -  GWAS hits

      -  Low-risk variants found in genome-wide association studies of
         cancer phenotypes (NHGRI-EBI Catalog)

4. **Documentation**

   -  Introduction

      -  Short overview of the predisposition report - aims and contents

   -  Annotation resources

      -  Underlying tools, databases and annotation sources (with
         versions)

   -  Variant classification

      -  Overview of how CPSR performs variant classification of
         variants not recorded in ClinVar, listing ACMG criteria and
         associated scores

   -  References

      -  Supporting scientific literature (Interpretation/implementation
         of ACMG critera etc.)

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
-  LOSS_OF_FUNCTION - Predicted loss-of-function variant
-  RMSK_HIT - Overlap with repeats as annotated by RepeatMasker
-  HGVSc - from VEP
-  HGVSp - from VEP
-  NCBI_REFSEQ - Transcript accession ID(s) (NCBI RefSeq)
-  ONCOGENE - Predicted proto-oncogene (CancerMine/NCG)
-  TUMOR_SUPPRESSOR - known tumor suppressor gene (CancerMine/NCG)
-  PREDICTED_EFFECT - Effect predictions (deleterious/benign) from
   dbNSFP
-  VEP_ALL_CSQ - All VEP transcript block consequences
-  DBSNP - dbSNP rsID
-  GENOMIC_CHANGE - Variant ID
-  GENOME_VERSION - Genome assembly

Example report
~~~~~~~~~~~~~~

-  `Cancer predisposition genome report
   report <http://insilico.hpc.uio.no/pcgr/example_reports/cpsr/0.6.2/SAMPLE-001.cpsr.grch37.html>`__

The HTML reports have been tested using the following browsers:

-  Safari (Version 14.1.1 (16611.2.7.1.4))
-  Mozilla Firefox (83.0)
-  Google Chrome (Version 90.0.4430.2127 (Official Build) (x86_64))

JSON
~~~~

A JSON file (gzipped) that stores the HTML report content is provided.
This file will easen the process of extracting particular parts of the
report for further analysis.

The JSON contains two main objects, *metadata* and *content*, where the
former contains information about the settings, data versions, and the
latter contains the various sections of the report.

This file can be used as input to
`PCGR <https://github.com/sigven/pcgr>`__, in order to populate a
somatic genome report with germline findings.

At present, there is no detailed schema documented for the JSON
structure.

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
   Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc
   \|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation\|
   ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID\|
   CANONICAL|MANE|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|RefSeq|DOMAINS\|
   HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|gnomAD_AF|gnomAD_AFR_AF\|
   gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF\|
   gnomAD_OTH_AF|gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO|CHECK_REF\|
   MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|NearestExonJB
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
-  SYMBOL_ENTREZ - Official gene symbol as provided by NCBI’s Entrez
   gene
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
-  NearestExonJB - VEP plugin that finds nearest exon junction for a
   coding sequence variant. Format: Ensembl exon identifier+distanceto
   exon boundary+boundary type(start/end)+exon length
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
-  INTRON_POSITION - Relative position of intron variant to nearest
   exon/intron junction (NearestExonJB VEP plugin)
-  EXON_POSITION - Relative position of exon variant to nearest
   intron/exon junction (NearestExonJB VEP plugin)
-  EXON - The exon number (out of total number) (picked by VEP’s
   –flag_pick_allele option)
-  LAST_EXON - Logical indicator for last exon of transcript (picked by
   VEP’s –flag_pick_allele option)
-  LAST_INTRON - Logical indicator for last intron of transcript (picked
   by VEP’s –flag_pick_allele option)
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
-  VEP_ALL_CSQ - All VEP transcript block consequences
   (Consequence:SYMBOL:Feature_type:Feature:BIOTYPE) - VEP
-  EXONIC_STATUS - Indicates if variant consequence type is ‘exonic’ or
   ‘nonexonic’. We define ‘exonic’ as any variants with the following
   consequences:

   -  stop_gained / stop_lost
   -  start_lost
   -  frameshift_variant
   -  missense_variant
   -  splice_donor_variant
   -  splice_acceptor_variant
   -  inframe_insertion / inframe_deletion
   -  synonymous_variant
   -  protein_altering

-  CODING_STATUS - Indicates if primary variant consequence type is
   ‘coding’ or ‘noncoding’. ‘coding’ variants are here defined as those
   with an ‘exonic’ status, with the exception of synonymous variants
-  NULL_VARIANT - Primary variant consequence type is frameshift or
   stop_gained/stop_lost
-  SPLICE_DONOR_RELEVANT - Logical indicating if variant is located at a
   particular location near the splice donor site (+3A/G, +4A or +5G)
-  REGULATORY_ANNOTATION - Comma-separated list of all variant
   annotations of Feature_type RegulatoryFeature and MotifFeature.
   Format: \|\|\|\|\|\|\|\|

*Gene information*
''''''''''''''''''

-  ENTREZ_ID - `Entrez <http://www.ncbi.nlm.nih.gov/gene>`__ gene
   identifier
-  APPRIS - Principal isoform flags according to the `APPRIS principal
   isoform database <http://appris.bioinfo.cnio.es/#/downloads>`__
-  MANE_SELECT - Indicating if the transcript is the MANE Select or MANE
   Plus Clinical transcript for the gene (picked by VEP’s
   –flag_pick_allele_gene option)
-  UNIPROT_ID - `UniProt <http://www.uniprot.org>`__ identifier
-  UNIPROT_ACC - `UniProt <http://www.uniprot.org>`__ accession(s)
-  ENSEMBL_GENE_ID - Ensembl gene identifier for VEP’s picked transcript
   (*ENSGXXXXXXX*)
-  ENSEMBL_TRANSCRIPT_ID - Ensembl transcript identifier for VEP’s
   picked transcript (*ENSTXXXXXX*)
-  ENSEMBL_PROTEIN_ID - Ensembl corresponding protein identifier for
   VEP’s picked transcript
-  REFSEQ_MRNA - Corresponding RefSeq transcript(s) identifier for VEP’s
   picked transcript (*NM_XXXXX*)
-  CORUM_ID - Associated protein complexes (identifiers) from
   `CORUM <http://mips.helmholtz-muenchen.de/corum/>`__
-  TUMOR_SUPPRESSOR - Indicates whether gene is predicted as a tumor
   suppressor gene, from Network of Cancer Genes (NCG) & the CancerMine
   text-mining resource
-  TUMOR_SUPPRESSOR_EVIDENCE - Underlying evidence for gene being a
   tumor suppressor. Format:
   NCG:<TRUE|FALSE>&CancerMine:<LC|MC|HC>:num_citations
-  ONCOGENE - Indicates whether gene is predicted as an oncogene, from
   Network of Cancer Genes (NCG) & the CancerMine text-mining resource
-  ONCOGENE_EVIDENCE - Underlying evidence for gene being an oncogene.
   Format: NCG:<TRUE|FALSE>&CancerMine:<LC|MC|HC>:num_citations
-  CANCER_SUSCEPTIBILITY_CUI - MedGen concept unique identifier (CUI)
   for cancer phenotype
-  CANCER_SYNDROME_CUI - MedGen concept unique identifier (CUI) for
   cancer syndrome
-  CANCER_PREDISPOSITION_SOURCE - Data source for susceptibility gene
   (panel *0*: NCGC, CGC_91, TCGA_PANCAN, PANEL_APP, OTHER)
-  CANCER_PREDISPOSITION_MOI - Mode of inheritance for susceptibility
   gene (AR/AD)
-  CANCER_PREDISPOSITION_MOD - Mechanism of disease for susceptibility
   gene (Lof/GoF)
-  PROB_EXAC_LOF_INTOLERANT - dbNSFP_gene: the probability of being
   loss-of-function intolerant (intolerant of both heterozygous and
   homozygous lof variants) based on ExAC r0.3 data
-  PROB_EXAC_LOF_INTOLERANT_HOM - dbNSFP_gene: the probability of being
   intolerant of homozygous, but not heterozygous lof variants based on
   ExAC r0.3 data
-  PROB_EXAC_LOF_TOLERANT_NULL - dbNSFP_gene: the probability of being
   tolerant of both heterozygous and homozygous lof variants based on
   ExAC r0.3 data
-  PROB_EXAC_NONTCGA_LOF_INTOLERANT - dbNSFP_gene: the probability of
   being loss-of-function intolerant (intolerant of both heterozygous
   and homozygous lof variants) based on ExAC r0.3 nonTCGA subset
-  PROB_EXAC_NONTCGA_LOF_INTOLERANT_HOM - dbNSFP_gene: the probability
   of being intolerant of homozygous, but not heterozygous lof variants
   based on ExAC r0.3 nonTCGA subset
-  PROB_EXAC_NONTCGA_LOF_TOLERANT_NULL - dbNSFP_gene: the probability of
   being tolerant of both heterozygous and homozygous lof variants based
   on ExAC r0.3 nonTCGA subset
-  PROB_GNOMAD_LOF_INTOLERANT - dbNSFP_gene: the probability of being
   loss-of-function intolerant (intolerant of both heterozygous and
   homozygous lof variants based on gnomAD 2.1 data
-  PROB_GNOMAD_LOF_INTOLERANT_HOM - dbNSFP_gene: the probability of
   being intolerant of homozygous, but not heterozygous lof variants
   based on gnomAD 2.1 data
-  PROB_GNOMAD_LOF_TOLERANT_NULL - dbNSFP_gene: the probability of being
   tolerant of both heterozygous and homozygous lof variants based on
   gnomAD 2.1 data
-  PROB_HAPLOINSUFFICIENCY - dbNSFP_gene: Estimated probability of
   haploinsufficiency of the gene (from
   http://dx.doi.org/10.1371/journal.pgen.1001154)
-  ESSENTIAL_GENE_CRISPR - dbNSFP_gene: Essential (E) or Non-essential
   phenotype-changing (N) based on large scale CRISPR experiments. from
   http://dx.doi.org/10.1126/science.aac7041
-  ESSENTIAL_GENE_CRISPR2 - dbNSFP_gene: Essential (E), context-Specific
   essential (S), or Non-essential phenotype-changing (N) based on large
   scale CRISPR experiments. from
   http://dx.doi.org/10.1016/j.cell.2015.11.015

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
   v4.2 <https://sites.google.com/site/jpopgen/dbNSFP>`__. Predicted
   effects are provided by different sources/algorithms (separated by
   ‘&’), T = Tolerated, N = Neutral, D = Damaging:

   1.  `SIFT <https://sift.bii.a-star.edu.sg/>`__
   2.  `MutationTaster <http://www.mutationtaster.org/>`__ (data release
       Nov 2015)
   3.  `MutationAssessor <http://mutationassessor.org/>`__ (release 3)
   4.  `FATHMM <http://fathmm.biocompute.org.uk>`__ (v2.3)
   5.  `PROVEAN <http://provean.jcvi.org/index.php>`__ (v1.1 Jan 2015)
   6.  `FATHMM_MKL <http://fathmm.biocompute.org.uk/fathmmMKL.htm>`__
   7.  `PRIMATEAI <https://www.nature.com/articles/s41588-018-0167-z>`__
   8.  `DEOGEN2 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5570203/>`__
   9.  `DBNSFP_CONSENSUS_RNN <https://www.biorxiv.org/content/10.1101/2021.04.09.438706v1>`__
       (Ensembl/consensus prediction, based on deep learning)
   10. `SPLICE_SITE_EFFECT_ADA <http://nar.oxfordjournals.org/content/42/22/13534>`__
       (Ensembl/consensus prediction of splice-altering SNVs, based on
       adaptive boosting)
   11. `SPLICE_SITE_EFFECT_RF <http://nar.oxfordjournals.org/content/42/22/13534>`__
       (Ensembl/consensus prediction of splice-altering SNVs, based on
       random forest)
   12. `M-CAP <http://bejerano.stanford.edu/MCAP>`__
   13. `MutPred <http://mutpred.mutdb.org>`__
   14. `GERP <http://mendel.stanford.edu/SidowLab/downloads/gerp/>`__
   15. `BayesDel <https://doi.org/10.1002/humu.23158>`__
   16. `LIST-S2 <https://doi.org/10.1093/nar/gkaa288>`__
   17. `ALoFT <https://www.nature.com/articles/s41467-017-00443-5>`__

-  DBNSFP_BAYESDEL_ADDAF - predicted effect from BayesDel (dbNSFP)

-  DBNSFP_LIST_S2 - predicted effect from LIST-S2 (dbNSFP)

-  DBNSFP_SIFT - predicted effect from SIFT (dbNSFP)

-  DBNSFP_PROVEAN - predicted effect from PROVEAN (dbNSFP)

-  DBNSFP_MUTATIONTASTER - predicted effect from MUTATIONTASTER (dbNSFP)

-  DBNSFP_MUTATIONASSESSOR - predicted effect from MUTATIONASSESSOR
   (dbNSFP)

-  DBNSFP_M_CAP - predicted effect from M-CAP (dbNSFP)

-  DBNSFP_ALOFTPRED - predicted effect from ALoFT (dbNSFP)

-  DBNSFP_MUTPRED - score from MUTPRED (dbNSFP)

-  DBNSFP_FATHMM - predicted effect from FATHMM (dbNSFP)

-  DBNSFP_PRIMATEAI - predicted effect from PRIMATEAI (dbNSFP)

-  DBNSFP_DEOGEN2 - predicted effect from DEOGEN2 (dbNSFP)

-  DBNSFP_GERP - evolutionary constraint measure from GERP (dbNSFP)

-  DBNSFP_FATHMM_MKL - predicted effect from FATHMM-mkl (dbNSFP)

-  DBNSFP_META_RNN - predicted effect from ensemble prediction (deep
   learning - dbNSFP)

-  DBNSFP_SPLICE_SITE_RF - predicted effect of splice site disruption,
   using random forest (dbscSNV)

-  DBNSFP_SPLICE_SITE_ADA - predicted effect of splice site disruption,
   using boosting (dbscSNV)

*Variant frequencies/annotations in germline databases*
'''''''''''''''''''''''''''''''''''''''''''''''''''''''

-  AFR_AF_GNOMAD - African/American germline allele frequency (`Genome
   Aggregation Database release
   2.1 <http://gnomad.broadinstitute.org/>`__)
-  AMR_AF_GNOMAD - American germline allele frequency (`Genome
   Aggregation Database release
   2.1 <http://gnomad.broadinstitute.org/>`__)
-  GLOBAL_AF_GNOMAD - Adjusted global germline allele frequency (`Genome
   Aggregation Database release
   2.1 <http://gnomad.broadinstitute.org/>`__)
-  SAS_AF_GNOMAD - South Asian germline allele frequency (`Genome
   Aggregation Database release
   2.1 <http://gnomad.broadinstitute.org/>`__)
-  EAS_AF_GNOMAD - East Asian germline allele frequency (`Genome
   Aggregation Database release
   2.1 <http://gnomad.broadinstitute.org/>`__)
-  FIN_AF_GNOMAD - Finnish germline allele frequency (`Genome
   Aggregation Database release
   2.1 <http://gnomad.broadinstitute.org/>`__)
-  NFE_AF_GNOMAD - Non-Finnish European germline allele frequency
   (`Genome Aggregation Database release
   2.1 <http://gnomad.broadinstitute.org/>`__)
-  OTH_AF_GNOMAD - Other germline allele frequency (`Genome Aggregation
   Database release 2.1 <http://gnomad.broadinstitute.org/>`__)
-  ASJ_AF_GNOMAD - Ashkenazi Jewish allele frequency (`Genome
   Aggregation Database release
   2.1 <http://gnomad.broadinstitute.org/>`__)
-  NON_CANCER_AF_ASJ - Alternate allele frequency for samples of
   Ashkenazi Jewish ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AF_EAS - Alternate allele frequency for samples of East
   Asian ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AF_AFR - Alternate allele frequency for samples of
   African-American/African ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AF_AMR - Alternate allele frequency for samples of Latino
   ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AF_OTH - Alternate allele frequency for samples of Other
   ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AF_NFE - Alternate allele frequency for samples of
   Non-Finnish European ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AF_FIN - Alternate allele frequency for samples of Finnish
   ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AF_SAS - Alternate allele frequency for samples of South
   Asian ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AF_GLOBAL - Alternate allele frequency in the non_cancer
   subset (`gnomAD 2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AC_ASJ - Alternate allele count for samples of Ashkenazi
   Jewish ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AC_EAS - Alternate allele count for samples of East Asian
   ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AC_AFR - Alternate allele count for samples of
   African-American/African ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AC_AMR - Alternate allele count for samples of Latino
   ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AC_OTH - Alternate allele count for samples of Other
   ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AC_NFE - Alternate allele frequency for samples of
   Non-Finnish European ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AC_FIN - Alternate allele count for samples of Finnish
   ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AC_SAS - Alternate allele count for samples of South Asian
   ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AC_GLOBAL - Alternate allele count in the non_cancer
   subset (`gnomAD 2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AN_ASJ - Total number of alleles in samples of Ashkenazi
   Jewish ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AN_EAS - Total number of alleles in samples of East Asian
   ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AN_AFR - Total number of alleles in samples of
   African-American/African ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AN_AMR - Total number of alleles in samples of Latino
   ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AN_OTH - Total number of alleles in samples of Other
   ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AN_NFE - Total number of alleles in samples of Non-Finnish
   European ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AN_FIN - Total number of alleles in samples of Finnish
   ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AN_SAS - Total number of alleles in samples of South Asian
   ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_AN_GLOBAL - Total number of alleles in the non_cancer
   subset (`gnomAD 2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_NHOMALT_ASJ - Count of homozygous individuals in samples
   of Ashkenazi Jewish ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_NHOMALT_EAS - Count of homozygous individuals in samples
   of East Asian ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_NHOMALT_AFR - Count of homozygous individuals in samples
   of African-American/African ancestry in the non_cancer subset
   (`gnomAD 2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_NHOMALT_AMR - Count of homozygous individuals in samples
   of Latino ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_NHOMALT_OTH - Count of homozygous individuals in samples
   of Other ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_NHOMALT_NFE - Count of homozygous individuals in samples
   of Non-Finnish European ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_NHOMALT_FIN - Count of homozygous individuals in samples
   of Finnish ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_NHOMALT_SAS - Count of homozygous individuals in samples
   of South Asian ancestry in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
-  NON_CANCER_NHOMALT_GLOBAL - Count of homozygous individuals in
   samples in the non_cancer subset (`gnomAD
   2.1.1 <http://gnomad.broadinstitute.org>`__)
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
-  CLINVAR_CONFLICTED - Variant has conflicting interpretations
-  CLINVAR_CLNSIG - Clinical significance for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ - germline
   state-of-origin
-  CLINVAR_CLASSIFICATION - Clean clinical significance on a five-level
   scheme
-  CLINVAR_CLNSIG_SOMATIC - Clinical significance for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ - somatic
   state-of-origin
-  CLINVAR_MEDGEN_CUI - Associated
   `MedGen <https://www.ncbi.nlm.nih.gov/medgen/>`__ concept identifiers
   (*CUIs*) - germline state-of-origin
-  CLINVAR_MEDGEN_CUI_SOMATIC - Associated
   `MedGen <https://www.ncbi.nlm.nih.gov/medgen/>`__ concept identifiers
   (*CUIs*) - somatic state-of-origin
-  CLINVAR_MOLECULAR_EFFECT - Variant effect according to ClinVar
   annotation
-  CLINVAR_VARIANT_ORIGIN - Origin of variant (somatic, germline, de
   novo etc.) for variant in
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__
-  CLINVAR_REVIEW_STATUS_STARS - Rating of the
   `ClinVar <http://www.ncbi.nlm.nih.gov/clinvar>`__ variant (0-4 stars)
   with respect to level of review
-  GWAS_HIT - variant associated with cancer phenotype from genome-wide
   association study (NHGRI-EBI GWAS catalog)
-  OPENTARGETS_DISEASE_ASSOCS - Associations between protein targets and
   disease based on multiple lines of evidence (mutations,affected
   pathways,GWAS, literature etc). Format:
   CUI:EFO_ID:IS_DIRECT:OVERALL_SCORE
-  OPENTARGETS_TRACTABILITY_COMPOUND - Confidence for the existence of a
   modulator (small molecule) that interacts with the target to elicit a
   desired biological effect
-  OPENTARGETS_TRACTABILITY_ANTIBODY - Confidence for the existence of a
   modulator (antibody) that interacts with the target to elicit a
   desired biological effect

Tab-separated values (TSV)
^^^^^^^^^^^^^^^^^^^^^^^^^^

Annotated List of all SNVs/InDels
'''''''''''''''''''''''''''''''''

We provide a tab-separated values file with most important variant/gene
annotations. The file has the following naming convention:

<*sample_id*>.\ **cpsr**.<*genome_assembly*>.\ **snvs_indels.tiers.tsv**

The SNVs/InDels are organized into different **tiers** (as defined above
for the HTML report)

The following variables are included in the tiered TSV file (VCF tags
issued by the user will be appended at the end):

::

   1. GENOMIC_CHANGE - Identifier for variant at the genome (VCF) level, e.g. 1:g.152382569A>G
         Format: (<chrom>:g.<position><ref_allele>><alt_allele>)
   2. VAR_ID - Variant identifier
   3. GENOTYPE - Variant genotype (heterozygous/homozygous)
   4. CPSR_CLASSIFICATION_SOURCE - ClinVar or Other (i.e. not present in ClinVar)
   5. GENOME_VERSION - Assembly version, e.g. GRCh37
   6. VCF_SAMPLE_ID - Sample identifier
   7. VARIANT_CLASS - Variant type, e.g. SNV/insertion/deletion
   8. CODING_STATUS - coding/noncoding (wrt. protein alteration and canonical splice site disruption)
   9. SYMBOL - Gene symbol
   10. GENE_NAME - Gene description
   11. CCDS - CCDS identifier
   12. ENTREZ_ID - Entrez gene identifier
   13. UNIPROT_ID - UniProt protein identifier
   14. ENSEMBL_GENE_ID - Ensembl gene identifier
   15. ENSEMBL_TRANSCRIPT_ID - Ensembl transcript identifier
   16. REFSEQ_MRNA - RefSeq mRNA identifier
   17. ONCOGENE - Gene is predicted as an oncogene according to Network of Cancer Genes (NCG) and CancerMine
   18. TUMOR_SUPPRESSOR - Gene is predicted as a tumor suppressor gene according to Network of Cancer Genes (NCG) and CancerMine
   19. CONSEQUENCE - Variant consequence
   20. VEP_ALL_CSQ - All VEP transcript block consequences
   21. PROTEIN_CHANGE - Protein change - one letter abbreviation (HGVSp)
   22. PROTEIN_DOMAIN - Protein domain (Pfam)
   23. DBSNP - dbSNP identifier (rsid)
   24. HGVSp - The HGVS protein sequence name
   25. HGVSc - The HGVS coding sequence name
   26. LAST_EXON - Last exon in gene
   27. EXON_POSITION - Relative position of exon variant to nearest intron/exon junction (NearestExonJB plugin)
   28. INTRON_POSITION - Relative position of intron variant to nearest intron/exon junction (NearestExonJB plugin)
   29. CDS_CHANGE - Coding, transcript-specific sequence annotation
   30. MUTATION_HOTSPOT - Cancer mutation hotspot (cancerhotspots.org)
   31. RMSK_HIT - RepeatMasker hit
   32. PROTEIN_FEATURE - Protein feature (active sites etc.) from UniProt KnowledgeBase
   33. EFFECT_PREDICTIONS - Functional effect predictions from multiple algorithms (dbNSFP)
   34. LOSS_OF_FUNCTION - Loss-of-function variant, as predicted from VEP's LofTee plugin
   35. CANCER_PHENOTYPE - For variants with a ClinVar classification, indication of cancer-associated disease/phenotype (1) or not (0)
   36. CLINVAR_CLASSIFICATION - clinical significance of ClinVar Variant (CPSR category)
   37. CLINVAR_MSID - measureset identifier of ClinVar variant
   38. CLINVAR_VARIANT_ORIGIN - variant origin (somatic/germline) of ClinVar variant
   39. CLINVAR_CONFLICTED - indicator of conflicting interpretations
   40. CLINVAR_PHENOTYPE - associated phenotype(s) for ClinVar variant
   41. CLINVAR_REVIEW_STATUS_STARS
   42. DBMTS - variant with potential effect on microRNA target sites (dbMTS). Format: <ensembl_transcript_id>|<microrna_identifier>|<target_prediction_algorithms>|<gain_loss_consensus>. _Target prediction algorithms_ indicate support by different algorithms (separated by '&'), 'TS' = TargetScan, M = 'miRanda', 'R' = RNAhybrid. *Gain_loss_consensus* indicate whether the variant was predicted to disrupt a binding site ('L' = Loss), or create a new target site ('G' = gain) by the different algorithms
   43. miRNA_TARGET_HIT - loss, gain, or gain|loss, as given from the hits in DBMTS column
   44. miRNA_TARGET_HIT_PREDICTION - links to miRBase, as given from the hits in the DBMTS column
   45. TF_BINDING_SITE_VARIANT - Indicates whether a variant overlaps a critical/non-critical position of a transcription factor binding site (TFBS) - as provided by VEP's--regulatory option ('Overlap: non-critical motif position' or 'Overlap: critical motif position')
   46. TF_BINDING_SITE_VARIANT_INFO - Comma-separated list of transcription factor binding sites affected by variant. Format per factor: <TRANSCRIPTION_FACTOR>|<MOTIF_NAME>|<MOTIF_POS>|<MOTIF_SCORE_CHANGE>|<HIGH_INF_POS>. *HIGH_INF_POS* indicates whether the variant overlapped a critical motif position ('Y'), or non-critical motif position ('N')
   47. GERP_SCORE - Genomic conservation score (GERP)
   48. N_INSILICO_CALLED - Number of algorithms with effect prediction (damaging/tolerated) from dbNSFP
   49. N_INSILICO_DAMAGING - Number of algorithms with damaging prediction from dbNSFP
   50. N_INSILICO_TOLERATED - Number of algorithms with tolerated prediction from dbNSFP
   51. N_INSILICO_SPLICING_NEUTRAL - Number of algorithms with splicing neutral prediction from dbscSNV
   52. N_INSILICO_SPLICING_AFFECTED - Number of algorithms with splicing affected prediction from dbscSNV
   53. GLOBAL_AF_GNOMAD - Global MAF in gnomAD
   54. <CUSTOM_POPULATION_GNOMAD> - Population specific MAF in gnomAD control (non-cancer, population configured by user)
   55. ACMG_BA1_AD - Very high MAF (> 0.5% in gnomAD non-cancer pop subset) - min AN = 12,000 - Dominant mechanism of disease
   56. ACMG_BS1_1_AD - High MAF (> 0.1% in gnomAD non-cancer pop subset) - min AN = 12,000 - Dominant mechanism of disease
   57. ACMG_BS1_2_AD - Somewhat high MAF (> 0.005% in gnomAD non-cancer pop subset) - Dominant mechanism of disease
   58. ACMG_BA1_AR - Very high MAF (> 1% in gnomAD non-cancer pop subset) - min AN = 12,000 - Recessive mechanism of disease
   59. ACMG_BS1_1_AR - High MAF (> 0.3% in gnomAD non-cancer pop subset) - min AN = 12,000 - Recessive mechanism of disease
   60. ACMG_BS1_2_AR - Somewhat high MAF (> 0.005% in gnomAD non-cancer pop subset) - Recessive mechanism of disease
   61. ACMG_PM2_1 - Allele count within pathogenic range (MAF <= 0.005% in the population-specific non-cancer gnomAD subset)
   62. ACMG_PM2_2 - Alternate allele absent in the population-specific non-cancer gnomAD subset
   63. ACMG_PVS1_1 - Null variant (frameshift/nonsense) - predicted as LoF by LOFTEE - within pathogenic range - LoF established for gene
   64. ACMG_PVS1_2 - Null variant (frameshift/nonsense) - not predicted as LoF by LOFTEE - within pathogenic range - LoF established for gene
   65. ACMG_PVS1_3 - Null variant (frameshift/nonsense) - predicted as LoF by LOFTEE - within pathogenic range - LoF not established for gene
   66. ACMG_PVS1_4 - Null variant (frameshift/nonsense) - not predicted as LoF by LOFTEE -- within pathogenic range - LoF not established for gene
   67. ACMG_PVS1_5 - Start (initiator methionine) lost - within pathogenic range - Lof established for gene
   68. ACMG_PVS1_6 - Start (initiator methionine) lost - within pathogenic range - LoF not established for gene
   69. ACMG_PVS1_7 - Donor/acceptor variant - predicted as LoF by LOFTEE - within pathogenic range - not last intron - LoF established for gene
   70. ACMG_PVS1_8 - Donor/acceptor variant - last intron - within pathogenic range - LoF established for gene
   71. ACMG_PVS1_9 - Donor/acceptor variant - not last intron - within pathogenic range - LoF not established for gene
   72. ACMG_PVS1_10 - Donor variant at located at the +3, +4 or +5 position of the intron -  within the pathogenic range (i.e. <9 alleles in ExAC))
   73. ACMG_PS1 - Same amino acid change as a previously established pathogenic variant (ClinVar) regardless of nucleotide change
   74. ACMG_PP2 - Missense variant in a gene that has a relatively low rate of benign missense variation (<20%) and where missense variants are a common mechanism of disease (>50% P/LP (ClinVar))
   75. ACMG_PM1 - Missense variant in a somatic mutation hotspot as determined by cancerhotspots.org
   76. ACMG_PM4 - Protein length changes due to inframe indels or nonstop variant in non-repetitive regions of genes that harbor variants with a dominant mode of inheritance.
   77. ACMG_PPC1 - Protein length changes due to inframe indels or nonstop variant in non-repetitive regions of genes that harbor variants with a recessive mode of inheritance.
   78. ACMG_PM5 - Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before (ClinVar)
   79. ACMG_PP3 - Multiple lines (>=5) of computational evidence support a deleterious effect on the gene or gene product (conservation, evolutionary, splicing impact) with maximum two contradictory predictions - from dbNSFP
   80. ACMG_BP4 - Multiple lines (>=5) of computational evidence support a benign effect on the gene or gene product (conservation, evolutionary, splicing impact) with maximum two contradictory prediction - from dbNSFP
   81. ACMG_BMC1 - Peptide change is at the same location of a known benign change (ClinVar)
   82. ACMG_BSC1 - Peptide change is reported as benign (ClinVar)
   83. ACMG_BP1 - Missense variant in a gene for which primarily truncating variants are known to cause disease (ClinVar)
   84. ACMG_BP3 - Variants in promoter or untranslated regions
   85. ACMG_BP7 - Silent/intronic variant outside of the splice site consensus
   86. FINAL_CLASSIFICATION - Final variant classification based on the combination of CLINVAR_CLASSIFICTION (for ClinVar-classified variants), and CPSR_CLASSIFICATION (for novel variants)
   87. CPSR_CLASSIFICATION - CPSR tier level (P/LP/VUS/LB/B)
   88. CPSR_PATHOGENICITY_SCORE - Aggregated CPSR pathogenicity score
   89. CPSR_CLASSIFICATION_CODE - Combination of CPSR classification codes assigned to the variant (ACMG)
   90. CPSR_CLASSIFICATION_DOC - Verbal description of CPSR classification codes assignted to the variant (ACMG)

**NOTE**: The user has the possibility to append the TSV file with data
from other INFO tags in the input VCF (i.e. using the
*–preserved_info_tags* option)
