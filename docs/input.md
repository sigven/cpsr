## Input

The CPSR workflow accepts a single input file:

  * An unannotated, single-sample **VCF file** (>= v4.2) with germline calls (SNVs/InDels)

### VCF

* We __strongly__ recommend that the input VCF is compressed and indexed using [bgzip](http://www.htslib.org/doc/tabix.html) and [tabix](http://www.htslib.org/doc/tabix.html)
* If the input VCF contains multi-allelic sites, these will be subject to [decomposition](http://genome.sph.umich.edu/wiki/Vt#Decompose)
* Variants used for reporting should be designated as 'PASS' in the VCF FILTER column

__IMPORTANT NOTE __: CPSR generates a number of VCF INFO annotation tags that is appended to the query VCF. We will therefore encourage the users to submit query VCF files that have not been subject to annotations by other means, but rather a VCF file that comes directly from variant calling. If not, there are likely to be INFO tags in the query VCF file that coincide with those produced by CPSR.
