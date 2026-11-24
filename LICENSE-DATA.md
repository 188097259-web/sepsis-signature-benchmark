# Licence for data, results and tables

Copyright (c) 2026 The sepsis-signature-benchmark authors

Unless stated otherwise below, the non-code material in this repository is licensed under the
Creative Commons Attribution 4.0 International licence (CC BY 4.0):
https://creativecommons.org/licenses/by/4.0/legalcode

This covers the authors' own material:

- `results/`: archived outputs of the analysis, including the archived random signatures
- `tables/`: source tables and the numbers cited in the article
- `docs/` and `data/signatures/`: the extraction tables and records compiled by the authors from the original reports
- `data/prior_work/`
- `params/` and `environment/`

Material that is not the authors' own:

- `data/cohorts/` holds processed copies of public Gene Expression Omnibus series (GSE65682, GSE185263, GSE63042,
  GSE236713, GSE272769, GSE95233 and GSE54514). These data remain subject to the terms under which the original
  investigators and GEO made them available. The CC BY 4.0 licence applies only to the processing added by the authors.
  Please cite the original studies when using these data.
- The gene lists and coefficients in `data/signatures/` were published by the authors of the original reports, which are
  identified by PubMed identifier. Please cite the original reports when using a signature.

The scripts in `scripts/` are licensed separately under the MIT licence (see `LICENSE`).
