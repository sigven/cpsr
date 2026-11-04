# Function that determines whether variants in cancer predisposition genes are linked to cancer-related conditions according to ClinVar

Each ClinVar trait of a variant (trait names in *CLINVAR_TRAITS*, and
UMLS concept names of *CLINVAR_UMLS_CUI*) is given one of three states:
cancer-related (1), not cancer-related (0) or unknown (NA). A trait is
cancer-related if it matches the cancer phenotype regex, maps to an
OncoTree concept, or is a '-related disorder' for a primary target gene
of the virtual panel. A trait is unknown if it is uninformative (e.g.
'not provided') or if a UMLS concept cannot be mapped to a name. The
variant-level state in *CLINVAR_PHENOTYPE_CANCER* is 1 if any trait is
cancer-related, otherwise 0 if any trait is not cancer-related,
otherwise NA (also NA for variants not found in ClinVar).

## Usage

``` r
is_clinvar_cancer_phenotype(cpg_calls, ref_data, primary_target_entrez = NULL)
```

## Arguments

- cpg_calls:

  data frame with variants in predisposition_genes

- ref_data:

  object with PCGR/CPSR reference data

- primary_target_entrez:

  character vector with Entrez gene identifiers of the primary targets
  (virtual panel genes); used to recognize gene-named conditions
  ('-related disorder') as cancer-related
