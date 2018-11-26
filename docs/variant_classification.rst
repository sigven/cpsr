Variant classification (ACMG)
-----------------------------

The unclassified, coding variants in the CPSR workflow have all been
assigned a composite pathogenicity score (coined **PATHSCORE** in the
tables below) according to the prioritization scheme outlined in
`CharGer <https://github.com/ding-lab/charger/>`__ (`Huang et al., Cell,
2018 <https://www.ncbi.nlm.nih.gov/pubmed/29625052>`__), as well as
`Maxwell et al., Am J Hum Genet,
2016 <https://www.ncbi.nlm.nih.gov/pubmed/27153395>`__. Specifially, a
cancer-specific adoption of previously established ACMG criteria for
variant classification has been implemented (`Richards et al.,
2015 <https://www.ncbi.nlm.nih.gov/pubmed/25741868>`__, `Amendola et
al., Am J Hum Genet,
2016 <https://www.ncbi.nlm.nih.gov/pubmed/27181684>`__). The following
ACMG evidence items indicate support for pathogenic/benign variants, and
have been implemented as follows (points in parenthesis indicate how
they contribute to overall **PATHSCORE**):

-  *PVS1* (**8**) - null variant (nonsense, frameshift, canonical ±1 or
   2 splice sites - indicated by `VEP’s LofTee
   plugin <https://github.com/konradjk/loftee>`__) in a gene where LoF
   is a known mechanism of disease (*dominant* mode of inheritance
   (MoI))
-  *PS1* (**7**) - Same amino acid change as a previously established
   pathogenic variant
   (`ClinVar <https://www.ncbi.nlm.nih.gov/clinvar/>`__) regardless of
   nucleotide change
-  *PSC1* (**4**) - null variant (nonsense, frameshift, canonical ±1 or
   2 splice sites - indicated by `VEP’s LofTee
   plugin <https://github.com/konradjk/loftee>`__) in a gene where *LoF*
   is a known mechanism of disease (*recessive* MoI)
-  *PMC1* (**2**) - null variant (nonsense, frameshift, canonical ±1 or
   2 splice sites - indicated by `VEP’s LofTee
   plugin <https://github.com/konradjk/loftee>`__) in a gene where *LoF*
   is NOT the known mechanism of disease
-  *PM1* (**2**) - Located in a somatic mutational hotspot
   (`cancerhotspots.org <https://www.cancerhotspots.org>`__)
-  *PM2* (**2**) - Absence/extremely low minor allele frequency (MAF <
   0.0005 in `1000 Genomes
   Project <http://www.internationalgenome.org/>`__/`gnomAD <http://gnomad.broadinstitute.org/>`__
   global population)
-  *PM4* (**2**) - Protein length changes due to inframe indels or
   stoploss variants in non-repetitive regions (as identified by
   `RepeatMasker <http://www.repeatmasker.org/>`__) of known
   susceptibility gene (*dominant* MoI)
-  *PM5* (**2**) - Novel missense change at an amino acid residue where
   a different missense change determined to be pathogenic has been seen
   before (`ClinVar <https://www.ncbi.nlm.nih.gov/clinvar/>`__)
-  *PP2* (**1**) - Missense variant in a gene that has a relatively low
   rate of benign missense variation and where missense variants are a
   common mechanism of disease. Current settings:

   -  Rate of benign missense variants in susceptibility gene is < 20%
   -  Rate of non-truncating pathogenic missense variants in
      susceptibility gene is > 50%

-  *PP3* (**1**) - Multiple lines of computational evidence support a
   deleterious effect on the gene or gene product (conservation,
   evolutionary, splicing impact, etc., from
   `dbNSFP <https://sites.google.com/site/jpopgen/dbNSFP>`__)
-  *PPC1* (**1**) - Protein length changes due to inframe indels or
   stoploss variants in non-repetitive regions (as identified by
   `RepeatMasker <http://www.repeatmasker.org/>`__) of known
   susceptibility gene (*recessive* MoI)
-  *BP4* (**-1**) - Multiple lines of computational evidence support a
   benign effect on the gene or gene product (conservation,
   evolutionary, splicing impact, etc., from
   `dbNSFP <https://sites.google.com/site/jpopgen/dbNSFP>`__)
-  *BMC1* (**-2**) - Peptide change is at the same location (codon) of a
   known benign change
   (`ClinVar <https://www.ncbi.nlm.nih.gov/clinvar/>`__)
-  *BSC1* (**-6**) - Peptide change is known to be benign
   (`ClinVar <https://www.ncbi.nlm.nih.gov/clinvar/>`__)
-  *BA1* (**-8**) - High allele frequency in the general population (MAF
   > 0.05 in `1000 Genomes
   Project <http://www.internationalgenome.org/>`__/`gnomAD <http://gnomad.broadinstitute.org/>`__
   global population)

   -  Exception for `homeostatic iron regulator
      (HFE) <https://www.ncbi.nlm.nih.gov/gene/3077>`__ and
      `SERPINA1 <https://www.ncbi.nlm.nih.gov/gene/5265>`__, requiring
      MAF > 0.25

For *PP3* and *BP4* CPSR requires by default that a minimum of five
distinct algorithms must constitute a majority vote (damaging/benign),
and that maximally one algorithm votes for the minority
(damaging/benign)

The composite **PATHSCORE** is finally assigned one of four different
**PATHRANK** levels:

-  HIGH: *PATHSCORE* > **8**
-  MODERATE: *PATHSCORE* <= **8** AND *PATHSCORE* > **4**
-  LOW: *PATHSCORE* <= **4** AND *PATHSCORE* >= **0**
-  BENIGN: *PATHSCORE* < **0**

The contribution of ACMG evidence items pr. variant can be seen in the
**PATHDOC** variable.
