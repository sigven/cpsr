CPSR is primarily using the non-cancer subset of gnomAD (v3.1.2) for variant allele
frequency estimates in its classification procedure. In the cases where there are
_no reported variation/allele frequencies_ in the non-cancer subset (gnomAD v3.1.2),
CPSR considers frequencies in the general and more comprehensive gnomAD dataset (v4.1.0)
as _a safeguard_, specifically for the application of the __PM2__, __BA1__ and __BS1__
criteria.
