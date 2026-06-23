# Input files

## VCF file

The CPSR workflow accepts a single input file, which is an unannotated,
single-sample **VCF file** (\>= v4.2) with germline calls (SNVs/InDels)

- We **strongly** recommend that the input VCF is compressed and indexed
  using [bgzip](http://www.htslib.org/doc/bgzip.md) and
  [tabix](http://www.htslib.org/doc/tabix.md)
- If the input VCF contains multi-allelic sites, these will be subject
  to [decomposition](http://genome.sph.umich.edu/wiki/Vt#Decompose).
  Either way, we encourage that users prepare the input VCF *without the
  presence of multi-allelic sites*.
- Variants used for reporting should be designated as ‘PASS’ in the VCF
  FILTER column
- In order for CPSR to report the genotypes of individual variants in
  cancer predisposition genes (**hom/het**), the VCF file needs to
  contain the sample-specific genotype data (i.e. a FORMAT field +
  sample-specific genotype field pr. VCF record)

**IMPORTANT NOTE** : CPSR generates a number of VCF INFO annotation tags
that is appended to the query VCF. We will therefore encourage the users
to submit query VCF files that have not been subject to annotations by
other means, but rather a VCF file that comes directly from variant
calling. If not, there are likely to be INFO tags in the query VCF file
that coincide with those produced by CPSR.
