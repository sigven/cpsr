# Output files

## Output

### Interactive HTML report

An interactive and structured [quarto](https://quarto.org/)-generated
HTML report, lists variants in known cancer predisposition genes and is
provided with the following naming convention:

- `<sample_id>.cpsr.<genome_assembly>.html`
  - The **sample_id** is provided as input by the user, and reflects a
    unique identifier of the sample to be analyzed.

The report is structured in multiple sections, described briefly below:

1.  **Settings**
    - Lists key configurations provided by user, including the list of
      genes that constitute the virtual gene panel in the report
2.  **Summary of findings**
    - Summarizes the main findings in the sample through value boxes
3.  **Variant classification**
    - For all coding and non-coding variants in the selected cancer
      predisposition geneset, interactive variant tables are shown, with
      variants highlighted by their clinical significance:
      - Pathogenic
        - Likely Pathogenic
        - Variants of Uncertain Significance (VUS)
        - Likely Benign
        - Benign
4.  **Genomic biomarkers**
    - Reported clinical evidence items from [CIViC](https://civicdb.org)
      that match with variants in the query set are reported
    - Pharmacogenetic findings (DPYD, TPMT, NUDT15)
5.  **Secondary findings**
    - Pathogenic variants in the [ACMG-recommended list of genes for
      report of secondary/incidental
      findings](https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/)
6.  **GWAS hits**
    - Status of relatively common, low-risk variants found in
      genome-wide association studies of cancer phenotypes (NHGRI-EBI
      Catalog)
7.  **Documentation**
    - Introduction
      - Short overview of the CPSR variant report - aims and contents
    - Annotation resources
      - Information on annotation sources utilized by CPSR, including
        versions and licensing requirements
    - Variant classification
      - Overview of how CPSR performs variant annotation and
        classification of variants not recorded in ClinVar, listing ACMG
        criteria and associated scores, calibration of classification
        thresholds etc.
8.  **References**
    - Supporting scientific literature - knowledge resources, guideline
      references etc.

  
  

### Variant call format - VCF

A VCF file containing annotated, germline variant calls (single
nucleotide variants and insertions/deletions) is generated with the
following naming convention:

- `<sample_id>.cpsr.<genome_assembly>.vcf.gz (.tbi)`
  - The **sample_id** is provided as input by the user, and reflects a
    unique identifier of the sample to be analyzed. Following common
    standards, the annotated VCF file is compressed with
    [bgzip](http://www.htslib.org/doc/bgzip.md) and indexed with
    [tabix](http://www.htslib.org/doc/tabix.md). Below follows a
    description of all annotations/tags present in the VCF INFO column
    after processing with the CPSR annotation pipeline:

  

##### *VEP consequence annotations*

| Tag | Description |
|----|----|
| `CSQ` | Complete consequence annotations from VEP. Format (separated by a `|`): `Allele`, `Consequence`, `IMPACT`, `SYMBOL`, `Gene`, `Feature_type`, `Feature`, `BIOTYPE`, `EXON`, `INTRON`, `HGVSc`, `HGVSp`, `cDNA_position`, `CDS_position`, `Protein_position`, `Amino_acids`, `Codons`, `Existing_variation`, `ALLELE_NUM`, `DISTANCE`, `STRAND`, `FLAGS`, `PICK`, `VARIANT_CLASS`, `SYMBOL_SOURCE`, `HGNC_ID`, `CANONICAL`, `MANE_SELECT`, `MANE_PLUS_CLINICAL`, `TSL`, `APPRIS`, `CCDS`, `ENSP`, `SWISSPROT`, `TREMBL`, `UNIPARC`, `UNIPROT_ISOFORM`, `RefSeq`, `DOMAINS`, `HGVS_OFFSET`, `gnomADe_AF`, `gnomADe_AFR_AF`, `gnomADe_AMR_AF`, `gnomADe_ASJ_AF`, `gnomADe_EAS_AF`, `gnomADe_FIN_AF`, `gnomADe_NFE_AF`, `gnomADe_OTH_AF`, `gnomADe_SAS_AF`, `CLIN_SIG`, `SOMATIC`, `PHENO`, `CHECK_REF`, `MOTIF_NAME`, `MOTIF_POS`, `HIGH_INF_POS`, `MOTIF_SCORE_CHANGE`, `TRANSCRIPTION_FACTORS`, `NearestExonJB` |
| `Consequence` | Impact modifier for the consequence type (picked by VEP’s `--flag_pick_allele` option) |
| `Gene` | Ensembl stable ID of affected gene (picked by VEP’s `--flag_pick_allele` option) |
| `Feature_type` | Type of feature. Currently one of `Transcript`, `RegulatoryFeature`, `MotifFeature` (picked by VEP’s `--flag_pick_allele` option) |
| `Feature` | Ensembl stable ID of feature (picked by VEP’s `--flag_pick_allele` option) |
| `cDNA_position` | Relative position of base pair in cDNA sequence (picked by VEP’s `--flag_pick_allele` option) |
| `CDS_position` | Relative position of base pair in coding sequence (picked by VEP’s `--flag_pick_allele` option) |
| `CDS_RELATIVE_POSITION` | Ratio of variant coding position to length of coding sequence |
| `CDS_CHANGE` | Coding, transcript-specific sequence annotation (picked by VEP’s `--flag_pick_allele` option) |
| `ALTERATION` | HGVSp/HGVSc identifier |
| `AMINO_ACID_START` | Protein position indicating absolute start of amino acid altered (fetched from `Protein_position`) |
| `AMINO_ACID_END` | Protein position indicating absolute end of amino acid altered (fetched from `Protein_position`) |
| `Protein_position` | Relative position of amino acid in protein (picked by VEP’s `--flag_pick_allele` option) |
| `Amino_acids` | Only given if the variant affects the protein-coding sequence (picked by VEP’s `--flag_pick_allele` option) |
| `GRANTHAM_DISTANCE` | Grantham distance between the reference and variant amino acids |
| `Codons` | The alternative codons with the variant base in upper case (picked by VEP’s `--flag_pick_allele` option) |
| `IMPACT` | Impact modifier for the consequence type (picked by VEP’s `--flag_pick_allele` option) |
| `VARIANT_CLASS` | Sequence Ontology variant class (picked by VEP’s `--flag_pick_allele` option) |
| `SYMBOL` | Gene symbol (picked by VEP’s `--flag_pick_allele` option) |
| `SYMBOL_SOURCE` | The source of the gene symbol (picked by VEP’s `--flag_pick_allele` option) |
| `STRAND` | The DNA strand (1 or -1) on which the transcript/feature lies (picked by VEP’s `--flag_pick_allele` option) |
| `ENSP` | The Ensembl protein identifier of the affected transcript (picked by VEP’s `--flag_pick_allele` option) |
| `FLAGS` | Transcript quality flags: `cds_start_NF`: CDS 5’, incomplete `cds_end_NF`: CDS 3’ incomplete (picked by VEP’s `--flag_pick_allele` option) |
| `SWISSPROT` | Best match UniProtKB/Swiss-Prot accession of protein product (picked by VEP’s `--flag_pick_allele` option) |
| `TREMBL` | Best match UniProtKB/TrEMBL accession of protein product (picked by VEP’s `--flag_pick_allele` option) |
| `UNIPARC` | Best match UniParc accession of protein product (picked by VEP’s `--flag_pick_allele` option) |
| `UNIPROT_ISOFORM` | Best match UniProtKB isoform accession of protein product (picked by VEP’s `--flag_pick_allele` option) |
| `HGVSc` | The HGVS coding sequence name (picked by VEP’s `--flag_pick_allele` option) |
| `HGVSc_RefSeq` | The HGVSc coding sequence name using RefSeq transcript identifiers (MANE select) - picked by VEP’s `--flag_pick_allele` option) |
| `HGVSp` | The HGVS protein sequence name (picked by VEP’s `--flag_pick_allele` option) |
| `HGVSp_short` | The HGVS protein sequence name, short version (picked by VEP’s `--flag_pick_allele` option) |
| `HGVS_OFFSET` | Indicates by how many bases the HGVS notations for this variant have been shifted (picked by VEP’s `--flag_pick_allele` option) |
| `NearestExonJB` | VEP plugin that finds nearest exon junction for a coding sequence variant. Format: `Ensembl exon identifier+distanceto exon boundary+boundary type(start/end)+exon length` |
| `MOTIF_NAME` | The source and identifier of a transcription factor binding profile aligned at this position (picked by VEP’s `--flag_pick_allele` option) |
| `MOTIF_POS` | The relative position of the variation in the aligned TFBP (picked by VEP’s `--flag_pick_allele` option) |
| `HIGH_INF_POS` | A flag indicating if the variant falls in a high information position of a transcription factor binding profile (TFBP) (picked by VEP’s `--flag_pick_allele` option) |
| `MOTIF_SCORE_CHANGE` | The difference in motif score of the reference and variant sequences for the TFBP (picked by VEP’s `--flag_pick_allele` option) |
| `CELL_TYPE` | List of cell types and classifications for regulatory feature (picked by VEP’s `--flag_pick_allele` option) |
| `CANONICAL` | A flag indicating if the transcript is denoted as the canonical transcript for this gene (picked by VEP’s `--flag_pick_allele` option) |
| `CCDS` | The CCDS identifier for this transcript, where applicable (picked by VEP’s `--flag_pick_allele` option) |
| `INTRON` | The intron number (out of total number) (picked by VEP’s `--flag_pick_allele` option) |
| `INTRON_POSITION` | Relative position of intron variant to nearest exon/intron junction (NearestExonJB VEP plugin) |
| `EXON_POSITION` | Relative position of exon variant to nearest intron/exon junction (NearestExonJB VEP plugin) |
| `EXON` | The exon number (out of total number) (picked by VEP’s `--flag_pick_allele` option) |
| `EXON_AFFECTED` | The exon affected by the variant (picked by VEP’s `--flag_pick_allele` option) |
| `LAST_EXON` | Logical indicator for last exon of transcript (picked by VEP’s `--flag_pick_allele` option) |
| `LAST_INTRON` | Logical indicator for last intron of transcript (picked by VEP’s `--flag_pick_allele` option) |
| `INTRON_POSITION` | Relative position of intron variant to nearest exon/intron junction (NearestExonJB VEP plugin) |
| `EXON_POSITION` | Relative position of exon variant to nearest intron/exon junction (NearestExonJB VEP plugin) |
| `DISTANCE` | Shortest distance from variant to transcript (picked by VEP’s `--flag_pick_allele` option) |
| `BIOTYPE` | Biotype of transcript or regulatory feature (picked by VEP’s `--flag_pick_allele` option) |
| `TSL` | Transcript support level (picked by VEP’s `--flag_pick_allele` option)\> |
| `PUBMED` | PubMed ID(s) of publications that cite existing variant - VEP |
| `PHENO` | Indicates if existing variant is associated with a phenotype, disease or trait - VEP |
| `GENE_PHENO` | Indicates if overlapped gene is associated with a phenotype, disease or trait - VEP |
| `ALLELE_NUM` | Allele number from input; 0 is reference, 1 is first alternate etc - VEP |
| `REFSEQ_MATCH` | The RefSeq transcript match status; contains a number of flags indicating whether this RefSeq transcript matches the underlying reference sequence and/or an Ensembl transcript (picked by VEP’s `--flag_pick_allele` option) |
| `PICK` | Indicates if this block of consequence data was picked by VEP’s `--flag_pick_allele` option |
| `VEP_ALL_CSQ` | All VEP transcript block consequences (`Consequence:SYMBOL:Feature_type:Feature:BIOTYPE`) - VEP |
| `EXONIC_STATUS` | Indicates if variant consequence type is ‘exonic’ or ‘nonexonic’. We define ‘exonic’ as any variants with the following consequence types: `stop_gained / stop_lost`, `start_lost`, `frameshift_variant`, `missense_variant`, `splice_donor_variant`, `splice_acceptor_variant`, `inframe_insertion / inframe_deletion`, `synonymous_variant`, `protein_altering` |
| `CODING_STATUS` | Indicates if primary variant consequence type is ‘coding’ or ‘noncoding’. ‘coding’ variants are here defined as those consequence types with an ‘exonic’ status, with the exception of synonymous variants. All other consequence types are considered ‘noncoding’ |
| `NULL_VARIANT` | Primary variant consequence type is `frameshift` or `stop_gained` |
| `LOSS_OF_FUNCTION` | Loss-of-function variant - primary variant consequence being either `stop_gained / stop_lost`, `start_lost`, `frameshift_variant`, `splice_donor_variant`, or `splice_acceptor_variant` |
| `LOF_FILTER` | Loss-of-function filter - exceptions to putative LOF variants - GC to GT at splice donor sites or truncations within the last 5% of coding sequence |
| `SPLICE_DONOR_RELEVANT` | Logical indicating if variant is located at a particular location near the splice donor site (`+3A/G`, `+4A` or `+5G`) |
| `BIOMARKER_MATCH` | Variant matches with germline biomarker evidence in CIViC/CGI. Format: `<db_source>|<db_variant_id>|<db_evidence_id>:<tumor_site>:<clinical_significance>:<evidence_level>:<evidence_type><germline_somatic>|<matching_type>`. Multiple evidence items are separated by ‘&’. Example: civic\|174\|EID445:Colon/Rectum:Sensitivity/Response:D:Predictive:Germline&EID446:Colon/Rectum:Sensitivity/Response:D:Predictive:Germline\|by_gene_mut. Matching type can be any of `by_genomic_coord`, `by_hgvsp_principal`, `by_hgvsc_principal`, `by_hgvsp_nonprincipal`, `by_hgvsc_nonprincipal`, `by_codon_principal`, `by_exon_mut_principal`, `by_gene_mut_lof`, `by_gene_mut` |
| `REGULATORY_ANNOTATION` | Comma-separated list of all variant annotations of `Feature_type`, `RegulatoryFeature`, and `MotifFeature`. Format (separated by a `|`): `<Consequence>`, `<Feature_type>`, `<Feature>`, `<BIOTYPE>`, `<MOTIF_NAME>`, `<MOTIF_POS>`, `<HIGH_INF_POS>`, `<MOTIF_SCORE_CHANGE>`, `<TRANSCRIPTION_FACTORS>` |

  

##### *Gene information*

| Tag | Description |
|----|----|
| `ENTREZGENE` | [Entrez](http://www.ncbi.nlm.nih.gov/gene) gene identifier |
| `APPRIS` | Principal isoform flags according to the [APPRIS principal isoform database](http://appris.bioinfo.cnio.es/#/downloads) |
| `MANE_SELECT` | Indicating if the transcript is the MANE Select for the gene (picked by VEP’s `--flag_pick_allele_gene` option) |
| `MANE_PLUS_CLINICAL` | Indicating if the transcript is MANE Plus Clinical, as required for clinical variant reporting (picked by VEP’s `--flag_pick_allele_gene` option) |
| `UNIPROT_ID` | [UniProt](http://www.uniprot.org) identifier |
| `UNIPROT_ACC` | [UniProt](http://www.uniprot.org) accession(s) |
| `ENSEMBL_GENE_ID` | Ensembl gene identifier for VEP’s picked transcript (*ENSGXXXXXXX*) |
| `ENSEMBL_TRANSCRIPT_ID` | Ensembl transcript identifier for VEP’s picked transcript (*ENSTXXXXXX*) |
| `ENSEMBL_PROTEIN_ID` | Ensembl corresponding protein identifier for VEP’s picked transcript |
| `REFSEQ_TRANSCRIPT_ID` | Corresponding RefSeq transcript(s) identifier for VEP’s picked transcript (*NM_XXXXX*) |
| `REFSEQ_PROTEIN_ID` | RefSeq protein/peptide identifier for VEP’s picked transcript (*NP_XXXXXX*) |
| `MANE_SELECT2` | MANE select transcript identifer: one high-quality representative transcript per protein-coding gene that is well-supported by experimental data and represents the biology of the gene - provided through BioMart |
| `MANE_PLUS_CLINICAL2` | transcripts chosen to supplement MANE Select when needed for clinical variant reporting - provided through BioMart |
| `GENCODE_TAG` | tag for GENCODE transcript (basic etc) |
| `GENCODE_TRANSCRIPT_TYPE` | type of transcript (protein-coding etc.) |
| `TSG` | Indicates whether gene is predicted as a tumor suppressor gene, from Network of Cancer Genes (NCG) & the CancerMine text-mining resource |
| `TSG_SUPPORT` | Underlying evidence for gene being a tumor suppressor. Format: `NCG&CancerMine:num_citations"` |
| `ONCOGENE` | Indicates whether gene is predicted as an oncogene, from Network of Cancer Genes (NCG) & the CancerMine text-mining resource |
| `ONCOGENE_SUPPORT` | Underlying evidence for gene being an oncogene. Format: `NCG&CancerMine:num_citations"` |
| `CPG_SOURCE` | Cancer predisposition gene source (panel *0*: `TCGA_PANCAN_2018`, `CANVAR_UK`, `PANEL_APP`, `CURATED_OTHER`) |
| `NCG_DRIVER` | Cancer driver gene prediction by Network of Cancer Genes (NCG) |
| `INTOGEN_DRIVER` | Indicates whether gene is predicted as cancer driver from IntOGen’s cancer driver prediction algorithm |
| `PROB_EXAC_LOF_INTOLERANT` | `dbNSFP_gene`: the probability of being loss-of-function intolerant (intolerant of both heterozygous and homozygous lof variants) based on ExAC r0.3 data |
| `PROB_EXAC_LOF_INTOLERANT_HOM` | `dbNSFP_gene`: the probability of being intolerant of homozygous, but not heterozygous lof variants based on ExAC r0.3 data |
| `PROB_EXAC_LOF_TOLERANT_NULL` | `dbNSFP_gene`: the probability of being tolerant of both heterozygous and homozygous lof variants based on ExAC r0.3 data |
| `PROB_EXAC_NONTCGA_LOF_INTOLERANT` | `dbNSFP_gene`: the probability of being loss-of-function intolerant (intolerant of both heterozygous and homozygous lof variants) based on ExAC r0.3 nonTCGA subset |
| `PROB_EXAC_NONTCGA_LOF_INTOLERANT_HOM` | `dbNSFP_gene`: the probability of being intolerant of homozygous, but not heterozygous lof variants based on ExAC r0.3 nonTCGA subset |
| `PROB_EXAC_NONTCGA_LOF_TOLERANT_NULL` | `dbNSFP_gene`: the probability of being tolerant of both heterozygous and homozygous lof variants based on ExAC r0.3 nonTCGA subset |
| `PROB_GNOMAD_LOF_INTOLERANT` | `dbNSFP_gene`: the probability of being loss-of-function intolerant (intolerant of both heterozygous and homozygous lof variants based on gnomAD 2.1 data |
| `PROB_GNOMAD_LOF_INTOLERANT_HOM` | `dbNSFP_gene`: the probability of being intolerant of homozygous, but not heterozygous lof variants based on gnomAD 2.1 data |
| `PROB_GNOMAD_LOF_TOLERANT_NULL` | `dbNSFP_gene`: the probability of being tolerant of both heterozygous and homozygous lof variants based on gnomAD 2.1 data |
| `PROB_HAPLOINSUFFICIENCY` | `dbNSFP_gene`: Estimated probability of haploinsufficiency of the gene (from <http://dx.doi.org/10.1371/journal.pgen.1001154>) |
| `ESSENTIAL_GENE_CRISPR` | `dbNSFP_gene`: Essential (`E`) or Non-essential phenotype-changing (`N`) based on large scale CRISPR experiments (from <http://dx.doi.org/10.1126/science.aac7041>) |
| `ESSENTIAL_GENE_CRISPR2` | `dbNSFP_gene`: Essential (`E`), context-Specific essential (`S`), or Non-essential phenotype-changing (`N`) based on large scale CRISPR experiments (from <http://dx.doi.org/10.1016/j.cell.2015.11.015>) |

  

##### *Variant effect and protein-coding information*

| Tag | Description |
|----|----|
| `MUTATION_HOTSPOT` | mutation hotspot codon in [cancerhotspots.org](http://cancerhotspots.org/). Format: `GeneSymbol|Entrez_ID|CodonRefAA|Alt_AA|Q-value` |
| `MUTATION_HOTSPOT_MATCH` | Type of hotspot match (by_hgvsp_principal, by_hgvsc_principal, by_hgvsp_nonprincipal, by_hgvsc_nonprincipal, by_codon_principal, by_codon_nonprincipal) |
| `MUTATION_HOTSPOT_CANCERTYPE` | hotspot-associated cancer types (from cancerhotspots.org) |
| `PFAM_DOMAIN` | Pfam domain identifier (from VEP) |
| `SPLICE_EFFECT` | Effect of splicing, from MutSpliceDB and/or MaxEntScan. Format: |
| MES |  |
| `EFFECT_PREDICTIONS` | Insilico predictions variant effect on protein function and pre-mRNA splicing from [database of non-synonymous functional predictions - dbNSFP v5.0](https://www.dbnsfp.org/). Predicted effects are provided by different sources/algorithms (separated by `&`), `T` = Tolerated, `N` = Neutral, `D` = Damaging |
| `DBNSFP_BAYESDEL_ADDAF` | predicted effect from BayesDel (dbNSFP) |
| `DBNSFP_LIST_S2` | predicted effect from LIST-S2 (dbNSFP) |
| `DBNSFP_SIFT` | predicted effect from SIFT (dbNSFP) |
| `DBNSFP_POLYPHEN2_HVAR` | predicted effect from PolyPhen2 (dbNSFP) |
| `DBNSFP_PROVEAN` | predicted effect from PROVEAN (dbNSFP) |
| `DBNSFP_MUTATIONTASTER` | predicted effect from MUTATIONTASTER (dbNSFP) |
| `DBNSFP_MUTATIONASSESSOR` | predicted effect from MUTATIONASSESSOR (dbNSFP) |
| `DBNSFP_M_CAP` | predicted effect from M-CAP (dbNSFP) |
| `DBNSFP_MUTPRED` | score from MUTPRED (dbNSFP) |
| `DBNSFP_CLINPRED` | predicted effect from ClinPred (dbNSFP) |
| `DBNSFP_FATHMM` | predicted effect from FATHMM-XF (dbNSFP) |
| `DBNSFP_PRIMATEAI` | predicted effect from PRIMATEAI (dbNSFP) |
| `DBNSFP_DEOGEN2` | predicted effect from DEOGEN2 (dbNSFP) |
| `DBNSFP_PHACTBOOST` | predicted effect from PHACTboost (dbNSFP) |
| `DBNSFP_ALPHA_MISSENSE` | predicted effect from AlphaMissense (dbNSFP) |
| `DBNSFP_MUTFORMER` | predicted effect from MutFormer (dbNSFP) |
| `DBNSFP_ESM1B` | predicted effect from ESM1b (dbNSFP) |
| `DBNSFP_GERP` | evolutionary constraint measure from GERP (dbNSFP) |
| `DBNSFP_CADD` | Combined Annotation Dependent Depletion (CADD) score (dbNSFP) |
| `DBNSFP_VEST4` | VEST4 score (dbNSFP) |
| `DBNSFP_FATHMM_XF` | predicted effect from FATHMM-XF (dbNSFP) |
| `DBNSFP_META_RNN` | predicted effect from ensemble prediction (deep learning - dbNSFP) |
| `DBNSFP_SPLICE_SITE_RF` | predicted effect of splice site disruption, using random forest (dbscSNV) |
| `DBNSFP_SPLICE_SITE_ADA` | predicted effect of splice site disruption, using boosting (dbscSNV) |

  

##### *Variant allele frequencies/annotations in germline databases*

| Tag | Description |
|----|----|
| `gnomADe_AFR_AF` | African/American germline allele frequency - exome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADe_AMR_AF` | Latino/Admixed American germline allele frequency - exome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADe_AF` | Adjusted global germline allele frequency - exome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADe_SAS_AF` | South Asian germline allele frequency - exome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADe_EAS_AF` | East Asian germline allele frequency - exome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADe_FIN_AF` | Finnish germline allele frequency - exome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADe_NFE_AF` | Non-Finnish European germline allele frequency - exome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADe_OTH_AF` | Other germline allele frequency - exome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADe_ASJ_AF` | Ashkenazi Jewish allele frequency - exome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADg_AFR_AF` | African/American germline allele frequency - genome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADg_AMR_AF` | Latino/Admixed American germline allele frequency - genome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADg_AF` | Adjusted global germline allele frequency - genome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADg_SAS_AF` | South Asian germline allele frequency - genome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADg_EAS_AF` | East Asian germline allele frequency - genome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADg_FIN_AF` | Finnish germline allele frequency - genome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADg_NFE_AF` | Non-Finnish European germline allele frequency - genome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADg_MID_AF` | Middle Eastern germline allele frequency - genome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADg_OTH_AF` | Other germline allele frequency - genome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADg_ASJ_AF` | Ashkenazi Jewish allele frequency - genome set ([gnomAD release 4.1](http://gnomad.broadinstitute.org/)) |
| `gnomADe_non_cancer_ASJ_AF` | Alternate allele frequency for samples of Ashkenazi Jewish ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_EAS_AF` | Alternate allele frequency for samples of East Asian ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_AFR_AF` | Alternate allele frequency for samples of African-American/African ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_AMR_AF` | Alternate allele frequency for samples of Latino ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_OTH_AF` | Alternate allele frequency for samples of Other ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_NFE_AF` | Alternate allele frequency for samples of Non-Finnish European ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_FIN_AF` | Alternate allele frequency for samples of Finnish ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_SAS_AF` | Alternate allele frequency for samples of South Asian ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_AF` | Alternate allele frequency in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_ASJ_AC` | Alternate allele count for samples of Ashkenazi Jewish ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_EAS_AC` | Alternate allele count for samples of East Asian ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_AFR_AC` | Alternate allele count for samples of African-American/African ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_AMR_AC` | Alternate allele count for samples of Latino ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_OTH_AC` | Alternate allele count for samples of Other ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_NFE_AC` | Alternate allele frequency for samples of Non-Finnish European ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_FIN_AC` | Alternate allele count for samples of Finnish ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_SAS_AC` | Alternate allele count for samples of South Asian ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_AC` | Alternate allele count in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_ASJ_AN` | Total number of alleles in samples of Ashkenazi Jewish ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_EAS_AN` | Total number of alleles in samples of East Asian ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_AFR_AN` | Total number of alleles in samples of African-American/African ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_AMR_AN` | Total number of alleles in samples of Latino ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_OTH_AN` | Total number of alleles in samples of Other ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_NFE_AN` | Total number of alleles in samples of Non-Finnish European ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_FIN_AN` | Total number of alleles in samples of Finnish ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_SAS_AN` | Total number of alleles in samples of South Asian ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_AN` | Total number of alleles in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_ASJ_NHOMALT` | Count of homozygous individuals in samples of Ashkenazi Jewish ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_EAS_NHOMALT` | Count of homozygous individuals in samples of East Asian ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_AFR_NHOMALT` | Count of homozygous individuals in samples of African-American/African ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_AMR_NHOMALT` | Count of homozygous individuals in samples of Latino ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_OTH_NHOMALT` | Count of homozygous individuals in samples of Other ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_NFE_NHOMALT` | Count of homozygous individuals in samples of Non-Finnish European ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_FIN_NHOMALT` | Count of homozygous individuals in samples of Finnish ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_SAS_NHOMALT` | Count of homozygous individuals in samples of South Asian ancestry in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `gnomADe_non_cancer_NHOMALT` | Count of homozygous individuals in samples in the non-cancer subset ([gnomAD 2.1.1](http://gnomad.broadinstitute.org)) |
| `DBSNP_RSID` | [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/) reference ID, as provided by VEP |

  

##### *Clinical associations*

| Tag | Description |
|----|----|
| `CLINVAR_MSID` | [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar) Measure Set/Variant ID |
| `CLINVAR_ALLELE_ID` | [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar) allele ID |
| `CLINVAR_PMID` | Associated Pubmed IDs for variant in [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar) - germline state-of-origin |
| `CLINVAR_HGVSP` | Protein variant expression using HGVS nomenclature - ClinVar |
| `CLINVAR_PMID_SOMATIC` | Associated Pubmed IDs for variant in [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar) - somatic state-of-origin |
| `CLINVAR_CONFLICTED` | ClinVar variant has conflicting interpretations |
| `CLINVAR_CLNSIG` | Clinical significance for variant in [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar) - germline state-of-origin |
| `CLINVAR_CLASSIFICATION` | Clean clinical significance on a five-level scheme - ClinVar |
| `CLINVAR_CONTRIB_CLNS_GERMLINE` | Contributing variant classifications from submissions - germline state-of-origin [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar) |
| `CLINVAR_CLNSIG_SOMATIC` | Clinical significance for variant in [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar) - somatic state-of-origin |
| `CLINVAR_MEDGEN_CUI` | Associated [MedGen](https://www.ncbi.nlm.nih.gov/medgen/) concept identifiers (*CUIs*) - germline state-of-origin |
| `CLINVAR_MEDGEN_CUI_SOMATIC` | Associated [MedGen](https://www.ncbi.nlm.nih.gov/medgen/) concept identifiers (*CUIs*) - somatic state-of-origin |
| `CLINVAR_MOLECULAR_EFFECT` | Variant effect according to ClinVar annotation |
| `CLINVAR_VARIANT_ORIGIN` | Origin of variant (somatic, germline, de novo etc.) for variant in [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar) |
| `CLINVAR_GOLD_STARS` | Rating of the [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar) variant (0-4 stars) with respect to level of review |
| `GWAS_HIT` | variant associated with cancer phenotype from genome-wide association study (NHGRI-EBI GWAS catalog) |

  

##### *Variant/genotype information*

| Tag          | Description                                  |
|--------------|----------------------------------------------|
| `GENOTYPE`   | Variant genotype (*het*/*hom_ref*/*hom_alt*) |
| `DP_CONTROL` | Sequencing depth at variant site (‘DP’)      |

  
  

### Excel workbook - XLSX

We provide an Excel workbook with **five** sheets that lists main
findings and annotations of the predisposition analysis. The file has
the following naming convention:

- `<sample_id>.cpsr.<genome_assembly>.xlsx`

The Excel workbook is populated with the following sheets (the last
three sheets will only be present if any data is found):

- *VIRTUAL_PANEL* - details on the the chosen virtual gene panel
- *CLASSIFICATION* - variant classifications and corresponding gene
  annotations
- *BIOMARKER_EVIDENCE* - matches of variants with genomic biomarkers
- *SECONDARY_FINDINGS* - secondary findings (ACMG recommendations)
- *PHARMACOGENETIC_FINDINGS* - drug toxicity findings (DPYD, TPMT,
  NUDT15)

  
  

### Tab-separated values - TSV

#### *Variant classification*

We provide a compressed tab-separated values file with variant
classifications and the most essential variant/gene annotations. The
file has the following naming convention:

- `<sample_id>.cpsr.<genome_assembly>.classification.tsv.gz`

The SNVs/InDels are classified according to clinical significance
(pathogenicity) (as defined above for the HTML report).

The following variables are included in the tiered TSV file (VCF tags in
the query VCF potentially retained by the user will be appended):

| Variable | Description |
|----|----|
| 1\. `SAMPLE_ID` | Sample identifier |
| 2\. `GENOMIC_CHANGE` | Identifier for variant at the genome (VCF) level, e.g. `1:g.152382569A>G`. Format: `<chrom>:g.<position><ref_allele>><alt_allele>` |
| 3\. `VAR_ID` | Variant identifier - chrom_pos_ref_alt |
| 4\. `GENOME_VERSION` | Assembly version, e.g. grch37/grch38 |
| 5\. `GENOTYPE` | Variant genotype (het/hom_ref/hom_alt) |
| 6\. `DP_CONTROL` | Sequencing depth at variant site (‘DP’) |
| 7\. `ASSERTION_AUTHORITY` | Authority responsible for variant classification - `ClinVar` (ClinVar record meets the configured trust level) or `CPSR` (classified by CPSR’s ACMG/AMP rule-based algorithm) |
| 8\. `ASSERTION_RATIONALE` | Human-readable explanation of why the given assertion authority was chosen, e.g. *“Novel variant (absent from ClinVar)”*, *“Conflicting ClinVar interpretations”*, *“ClinVar review status: 2 gold stars”*, or *“No cancer-related ClinVar phenotype(s)”* |
| 9\. `VARIANT_CLASS` | Variant type, e.g. SNV/insertion/deletion |
| 10\. `CODING_STATUS` | coding/noncoding (wrt. protein alteration and canonical splice site disruption) |
| 11\. `SYMBOL` | Gene symbol |
| 12\. `GENENAME` | Gene description |
| 13\. `CCDS` | CCDS identifier |
| 14\. `ENTREZGENE` | Entrez gene identifier |
| 15\. `UNIPROT_ID` | UniProt protein identifier |
| 16\. `ENSEMBL_GENE_ID` | Ensembl gene identifier |
| 17\. `ENSEMBL_TRANSCRIPT_ID` | Ensembl transcript identifier |
| 18\. `REFSEQ_TRANSCRIPT_ID` | RefSeq mRNA identifier |
| 19\. `MANE_SELECT` | Indicating if the transcript is the MANE Select for the gene (VEP) |
| 20\. `MANE_SELECT2` | MANE Select transcript identifier provided through BioMart |
| 21\. `ONCOGENE` | Gene is predicted as an oncogene according to Network of Cancer Genes (NCG) and/or CancerMine |
| 22\. `TUMOR_SUPPRESSOR` | Gene is predicted as a tumor suppressor gene according to Network of Cancer Genes (NCG) and/or CancerMine |
| 23\. `CPG_MOD` | Gene - cancer predisposition mechanism of disease (e.g. LoF) |
| 24\. `CPG_MOI` | Gene - cancer predisposition mode of inheritance |
| 25\. `CONSEQUENCE` | Variant consequence |
| 26\. `ALTERATION` | Molecular alteration (HGVSp or HGVSc pending on consequence) |
| 27\. `PROTEIN_CHANGE` | Protein change - one letter abbreviation (HGVSp) |
| 28\. `PFAM_DOMAIN` | Protein domain (Pfam identifier) |
| 29\. `PFAM_DOMAIN_NAME` | Protein domain name (Pfam) |
| 30\. `HGVSp` | The HGVS protein sequence name |
| 31\. `GRANTHAM_DISTANCE` | Grantham distance for amino acid change (Grantham score) |
| 32\. `HGVSc` | The HGVS coding sequence name |
| 33\. `HGVSc_RefSeq` | The HGVS coding sequence name (RefSeq - MANE Select) |
| 34\. `CDS_CHANGE` | Coding, transcript-specific sequence annotation |
| 35\. `LAST_EXON` | Last exon in gene |
| 36\. `LAST_INTRON` | Last intron in gene |
| 37\. `EXON` | Exon of variant/total number of exons in transcript (from VEP) |
| 38\. `EXON_AFFECTED` | Transcript exon of variant (from VEP) |
| 39\. `EXON_POSITION` | Relative position of exon variant to nearest intron/exon junction (NearestExonJB plugin) |
| 40\. `INTRON_POSITION` | Relative position of intron variant to nearest intron/exon junction (NearestExonJB plugin) |
| 41\. `NMD` | Nonsense-mediated decay prediction for the variant |
| 42\. `EXON_INTRON_JUNCTION_SPAN` | Indicates whether the variant spans an exon-intron junction boundary |
| 43\. `EXONIC_STATUS` | Indicates if variant consequence type is ‘exonic’ or ‘nonexonic’ |
| 44\. `PROTEIN_RELATIVE_POSITION` | Relative position of the affected amino acid within the protein (ratio of protein position to total protein length) |
| 45\. `MUTATION_HOTSPOT` | Cancer mutation hotspot (cancerhotspots.org) |
| 46\. `RMSK_HIT` | RepeatMasker hit |
| 47\. `EFFECT_PREDICTIONS` | Functional effect predictions from multiple algorithms (dbNSFP) |
| 48\. `MAXENTSCAN` | MaxEntScan splice site scores and effect predictions |
| 49\. `SPLICE_EFFECT` | Splice effect annotations from MutSpliceDB and MaxEntScan (see details above) |
| 50\. `LOSS_OF_FUNCTION` | Loss-of-function variant |
| 51\. `LOF_FILTER` | Loss-of-function filter |
| 52\. `NULL_VARIANT` | Frameshift or stop-gain variant |
| 53\. `DBMTS` | Variant with potential effect on microRNA target sites (dbMTS). Format: `<ensembl_transcript_id>|<microrna_identifier>|<target_prediction_algorithms>|<gain_loss_consensus>`. *Target prediction algorithms* indicate support by different algorithms (separated by ‘&’), `TS` = TargetScan, `M` = miRanda, `R` = RNAhybrid. *Gain_loss_consensus* indicates whether the variant was predicted to disrupt a binding site (`L` = Loss), or create a new target site (`G` = Gain) |
| 54\. `REGULATORY_ANNOTATION` | Overlap of variant with regulatory elements (VEP) |
| 55\. `TF_BINDING_SITE_VARIANT` | Indicates whether a variant overlaps a critical/non-critical position of a transcription factor binding site (TFBS) |
| 56\. `TF_BINDING_SITE_VARIANT_INFO` | Comma-separated list of transcription factor binding sites affected by variant. Format per factor: `<TRANSCRIPTION_FACTOR>|<MOTIF_NAME>|<MOTIF_POS>|<MOTIF_SCORE_CHANGE>|<HIGH_INF_POS>`. *HIGH_INF_POS* indicates whether the variant overlapped a critical motif position (`Y`) or non-critical motif position (`N`) |
| 57\. `VEP_ALL_CSQ` | All VEP transcript block consequences |
| 58\. `GERP_SCORE` | Genomic conservation score (GERP) |
| 59\. `DBSNP_RSID` | dbSNP identifier (rsid) |
| 60\. `CLINVAR_CLASSIFICATION` | Clinical significance of ClinVar-recorded variant (five-level scheme) |
| 61\. `CLINVAR_MSID` | Measureset identifier of ClinVar variant |
| 62\. `CLINVAR_VARIANT_ORIGIN` | Variant origin (somatic/germline) of ClinVar variant |
| 63\. `CLINVAR_CONFLICTED` | Indicator of conflicting interpretations in ClinVar |
| 64\. `CLINVAR_PHENOTYPE` | Associated phenotype(s) for ClinVar variant |
| 65\. `CLINVAR_PHENOTYPE_CANCER` | For variants with a ClinVar classification, indication of cancer-associated disease/phenotype (1) or not (0) |
| 66\. `CLINVAR_GOLD_STARS` | Review confidence rating of the ClinVar variant (0–4 gold stars) |
| 67\. `N_INSILICO_CALLED` | Number of algorithms with effect prediction (damaging/tolerated) from dbNSFP |
| 68\. `N_INSILICO_DAMAGING` | Number of algorithms with damaging prediction from dbNSFP |
| 69\. `N_INSILICO_TOLERATED` | Number of algorithms with tolerated prediction from dbNSFP |
| 70\. `N_INSILICO_SPLICING_NEUTRAL` | Number of algorithms with splicing neutral prediction from dbscSNV |
| 71\. `N_INSILICO_SPLICING_AFFECTED` | Number of algorithms with splicing affected prediction from dbscSNV |
| 72\. `gnomADe_AF` | Global allele frequency in gnomAD controls (exome set, v4.1) |
| 73\. `gnomADg_AF` | Global allele frequency in gnomAD controls (genome set, v4.1) |
| 74\. `CLASSIFICATION` | Final variant classification (P/LP/VUS/LB/B), reflecting the applicable assertion authority: `CLINVAR_CLASSIFICATION` if the ClinVar record meets the configured trust level, otherwise `CPSR_CLASSIFICATION` |
| 75\. `CPSR_CLASSIFICATION` | Variant clinical significance by CPSR’s ACMG/AMP rule-based algorithm (P/LP/VUS/LB/B) |
| 76\. `CPSR_PATHOGENICITY_SCORE` | Aggregated pathogenicity score by CPSR’s algorithm |
| 77\. `ACMG_CODE` | Combination of ACMG/AMP evidence codes assigned to the variant by CPSR |

**NOTE**: The user has the possibility to append the TSV file with data
from other INFO tags in the input VCF (i.e. using the
*–retained_info_tags* option)

#### *Biomarker evidence*

We provide a compressed tab-separated values file with variants
implicated as germline biomarkers. The file has the following naming
convention:

- `<sample_id>.cpsr.<genome_assembly>.biomarker_evidence.tsv.gz`

The following variables are included in the biomarker evidence TSV file:

| Variable | Description |
|----|----|
| 1\. `SAMPLE_ID` | Sample identifier |
| 2\. `GENOMIC_CHANGE` | Identifier for variant at the genome (VCF) level, e.g. `1:g.152382569A>G` |
| 3\. `GENOTYPE` | Variant genotype (het/hom_ref/hom_alt) |
| 4\. `DP_CONTROL` | Sequencing depth at variant site (‘DP’) |
| 5\. `GENOME_VERSION` | Assembly version, e.g. grch37/grch38 |
| 6\. `VARIANT_CLASS` | Variant type, e.g. SNV/insertion/deletion |
| 7\. `SYMBOL` | Gene symbol |
| 8\. `GENENAME` | Gene description |
| 9\. `CONSEQUENCE` | Variant consequence |
| 10\. `PROTEIN_CHANGE` | Protein change - one letter abbreviation (HGVSp) |
| 11\. `ASSERTION_AUTHORITY` | Authority responsible for variant classification - `ClinVar` or `CPSR` |
| 12\. `ASSERTION_RATIONALE` | Human-readable explanation of why the given assertion authority was chosen |
| 13\. `CLASSIFICATION` | Final variant classification (P/LP/VUS/LB/B) |
| 14\. `CPSR_CLASSIFICATION` | Variant clinical significance by CPSR’s ACMG/AMP rule-based algorithm (P/LP/VUS/LB/B) |
| 15\. `CPSR_PATHOGENICITY_SCORE` | Aggregated pathogenicity score by CPSR’s algorithm |
| 16\. `ACMG_CODE` | Combination of ACMG/AMP evidence codes assigned to the variant by CPSR |
| 17\. `CLINVAR_CLASSIFICATION` | Clinical significance of ClinVar-recorded variant (five-level scheme) |
| 18\. `BM_CANCER_TYPE` | Annotated cancer type for biomarker - from CIViC |
| 19\. `BM_DISEASE_ONTOLOGY_ID` | Disease ontology id for cancer type - from CIViC |
| 20\. `BM_PRIMARY_SITE` | Primary tumor type of cancer type - mapped with [phenOncoX](https://github.com/sigven/phenOncoX) |
| 21\. `BM_CLINICAL_SIGNIFICANCE` | Clinical significance of biomarker (drug sensitivity, drug resistance, poor outcome etc.) - from CIViC |
| 22\. `BM_THERAPEUTIC_CONTEXT` | Cancer drugs associated with biomarker (for biomarkers related to drug sensitivity/resistance) - from CIViC |
| 23\. `BM_CITATION` | Reference/source for biomarker - i.e. publication or guidelines - from CIViC |
| 24\. `BM_RATING` | Rating of biomarker - from CIViC |
| 25\. `BM_MOLECULAR_PROFILE` | Associated name of molecular profile - i.e. “BRCA mutation” - from CIViC |
| 26\. `BM_EVIDENCE_TYPE` | Biomarker type - *Predictive*, *Diagnostic*, *Prognostic*, *Predisposing* - from CIViC |
| 27\. `BM_EVIDENCE_LEVEL` | Strength of evidence for the given biomarker - *A* to *D* - from CIViC |
| 28\. `BM_EVIDENCE_DIRECTION` | Direction of biomarker evidence, i.e. *Supports* or *Does Not Support* - from CIViC |
| 29\. `BM_EVIDENCE_DESCRIPTION` | Description of biomarker - from CIViC |
| 30\. `BM_SOURCE_DB` | Biomarker source database - CIViC |
| 31\. `BM_EVIDENCE_ID` | Evidence identifier - from CIViC |
| 32\. `BM_VARIANT_ORIGIN` | Origin of biomarker variant - *germline* |
| 33\. `BM_MATCH` | Match between sample variant and biomarker - *by_genomic_coord*, *by_hgvsp_principal*, *by_gene_mut_lof* etc. |
| 34\. `BM_RESOLUTION` | Highest resolution of mapping between sample variant and biomarker - *genomic*, *hgvsp*, *codon*, *gene* |

#### *Pharmacogenetic findings*

We provide a compressed tab-separated values file with variants
implicated with drug toxicity/dosage effects of cancer chemotherapies.
The file has the following naming convention:

- `<sample_id>.cpsr.<genome_assembly>.pgx_findings.tsv.gz`

The following variables are included in the pharmacogenetic findings TSV
file:

| Variable | Description |
|----|----|
| 1\. `SAMPLE_ID` | Sample identifier |
| 2\. `GENOMIC_CHANGE` | Identifier for variant at the genome (VCF) level, e.g. `1:g.152382569A>G` |
| 3\. `GENOTYPE` | Variant genotype (het/hom_ref/hom_alt) |
| 4\. `DP_CONTROL` | Sequencing depth at variant site (‘DP’) |
| 5\. `GENOME_VERSION` | Assembly version, e.g. grch37/grch38 |
| 6\. `VARIANT_CLASS` | Variant type, e.g. SNV/insertion/deletion |
| 7\. `SYMBOL` | Gene symbol |
| 8\. `GENENAME` | Gene description |
| 9\. `CONSEQUENCE` | Variant consequence |
| 10\. `PROTEIN_CHANGE` | Protein change - one letter abbreviation (HGVSp) |
| 11\. `CLINVAR_CLASSIFICATION` | Clinical significance of ClinVar-recorded variant (five-level scheme) |
| 12\. `CLINVAR_MSID` | Measureset identifier of ClinVar variant |
| 13\. `CLINVAR_VARIANT_ORIGIN` | Variant origin (somatic/germline) of ClinVar variant |
| 14\. `CLINVAR_CONFLICTED` | Indicator of conflicting interpretations in ClinVar |
| 15\. `CLINVAR_PHENOTYPE` | Associated phenotype(s) for ClinVar variant |
| 16\. `CLINVAR_GOLD_STARS` | Review confidence rating of the ClinVar variant (0–4 gold stars) |
| 17\. `ENSEMBL_GENE_ID` | Ensembl gene identifier |
| 18\. `ENSEMBL_TRANSCRIPT_ID` | Ensembl transcript identifier |
| 19\. `REFSEQ_TRANSCRIPT_ID` | RefSeq mRNA identifier |
| 20\. `PFAM_DOMAIN_NAME` | Protein domain name (Pfam) |
| 21\. `HGVSp` | The HGVS protein sequence name |
| 22\. `HGVSc` | The HGVS coding sequence name |
| 23\. `HGVSc_RefSeq` | The HGVS coding sequence name (RefSeq - MANE Select) |
| 24\. `CDS_CHANGE` | Coding, transcript-specific sequence annotation |
| 25\. `CODING_STATUS` | coding/noncoding (wrt. protein alteration and canonical splice site disruption) |
| 26\. `MUTATION_HOTSPOT` | Cancer mutation hotspot (cancerhotspots.org) |
| 27\. `EFFECT_PREDICTIONS` | Functional effect predictions from multiple algorithms (dbNSFP) |
| 28\. `SPLICE_EFFECT` | Splice effect annotations from MutSpliceDB and MaxEntScan |
| 29\. `LOSS_OF_FUNCTION` | Loss-of-function variant |
| 30\. `LOF_FILTER` | Loss-of-function filter |
| 31\. `NULL_VARIANT` | Frameshift or stop-gain variant |
| 32\. `DBSNP_RSID` | dbSNP identifier (rsid) |
| 33\. `gnomADe_AF` | Global allele frequency in gnomAD controls (exome set, v4.1) |
| 34\. `gnomADg_AF` | Global allele frequency in gnomAD controls (genome set, v4.1) |

  
  

### Biomarker annotations

The TSV biomarker evidence output, the interactive HTML report (section
*Genomic biomarkers*), and the Excel workbook (sheet
*BIOMARKER_EVIDENCE*), contain information on matches between potential
pathogenic/likely pathogenic sample variants and reported biomarkers,
the latter referring to clinical evidence items that relate genomic
genomic aberrations to prognosis, diagnosis or sensitivity/resistance to
particular treatments. All biomarker annotations are prefixed with
**BM\_**, and the following is provided per evidence item:

| Variable | Description |
|----|----|
| 1\. `BM_CANCER_TYPE` | Annotated cancer type for biomarker - from CIViC |
| 2\. `BM_DISEASE_ONTOLOGY_ID` | Disease ontology id for cancer type - from CIViC |
| 3\. `BM_PRIMARY_SITE` | Primary tumor type of cancer type - mapped with [phenOncoX](https://github.com/sigven/phenOncoX) |
| 4\. `BM_CLINICAL_SIGNIFICANCE` | Clinical significance of biomarker (drug sensitivity, drug resistance, poor outcome etc.) - from CIViC |
| 5\. `BM_THERAPEUTIC_CONTEXT` | Cancer drugs associated with biomarker (for biomarkers related to drug sensitivity/resistance) - from CIViC |
| 6\. `BM_CITATION` | Reference/source for biomarker - i.e. publication or guidelines - from CIViC |
| 7\. `BM_RATING` | Rating of biomarker - from CIViC |
| 8\. `BM_MOLECULAR_PROFILE` | Associated name of molecular profile - i.e. “BRCA mutation” - from CIViC |
| 9\. `BM_EVIDENCE_TYPE` | Biomarker type - *Predictive*, *Diagnostic*, *Prognostic*, *Predisposing* - from CIViC |
| 10\. `BM_EVIDENCE_LEVEL` | Strength of evidence for the given biomarker - *A* to *D* - from CIViC |
| 11\. `BM_EVIDENCE_DIRECTION` | Direction of biomarker evidence, i.e. *Supports* or *Does Not Support* - from CIViC |
| 12\. `BM_EVIDENCE_DESCRIPTION` | Description of biomarker - from CIViC |
| 13\. `BM_SOURCE_DB` | Biomarker source database - CIViC |
| 14\. `BM_EVIDENCE_ID` | Evidence identifier - from CIViC |
| 15\. `BM_VARIANT_ORIGIN` | Origin of biomarker variant - *germline* |
| 16\. `BM_MATCH` | Match between sample variant and biomarker - *by_genomic_coord*, *by_hgvsp_principal*, *by_gene_mut_lof* etc. |
| 17\. `BM_RESOLUTION` | Highest resolution of mapping between sample variant and biomarker - *genomic*, *hgvsp*, *codon*, *gene* |
