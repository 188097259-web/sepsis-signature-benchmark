# Sepsis blood transcriptomic signature benchmark

Code, input data, archived random signatures and result tables for a benchmark of 36 published sepsis prognostic gene signatures against random signatures built with a common modelling pipeline, with two further analyses: how much of each signature score is explained by leukocyte composition, and whether the labels assigned by endotyping tools depend on the other samples in the analysis batch.

## What is in this repository

| Folder | Content |
|---|---|
| `scripts/` | All analysis and figure scripts (R). The table under “How to run” gives the purpose of each script |
| `params/` | Parameters of the random-signature pipeline (`pipeline_params.json`), read by the scripts |
| `docs/` | `pipeline_extraction.csv`: modelling steps extracted from the 36 original reports, from which the pipeline parameters were taken |
| `data/cohorts/` | Processed expression matrices and sample sheets of the seven public cohorts |
| `data/signatures/` | Definitions of the published signatures and endotyping tools (genes, coefficients, direction, eligibility per cohort) and the extraction records of the 36 reports |
| `data/prior_work/` | Pooled AUROCs of the published signatures that the authors had computed before the analysis plan was written; input to the precision simulation |
| `data/external/`, `data/cmaise/` | Empty. Places for two inputs that are not distributed here (see below) |
| `results/` | Archived outputs of every script, including all random signatures (`results/part1/random_models.rds`) |
| `tables/` | Source tables of the article (Table 1, Table 2, signature key) and, in `tables/numbers/`, every number cited in the text with its source file |
| `environment/` | `sessionInfo()` of each script at the time of the archived run |

## What is not in this repository

- **CMAISE patient-level data** (counts and clinical table of the Chinese multicentre cohort). These data are under controlled access; see the data availability statement of the article. Only aggregate CMAISE results are archived in `results/cmaise/`. Script `08_cmaise_validation.R` needs the two files `counts.tsv.gz` and `clinical_samples.csv`; set the environment variable `CMAISE_DIR` to the folder that holds them, or place them in `data/cmaise/`.
- **LM22 signature matrix.** Obtain `LM22.txt` from the supplementary material of Newman et al., *Nat Methods* 2015;12:453–457, and either place it at `data/external/LM22.txt` or set the environment variable `LM22_FILE`. The scripts check its MD5 (`f315a17039596ef73f83fbb7466677dd`).
- **Raw GEO downloads.** The processed matrices in `data/cohorts/` were derived from the GEO series listed below.

## Cohorts

| File in `data/cohorts/` | GEO series | Role |
|---|---|---|
| `GSE65682.rds`, `GSE65682_series_matrix_header.txt.gz` | GSE65682 | Training of random signatures. The second file is the header of the GEO series matrix (sample characteristics, no expression values), from which the 28-day outcome and time are read |
| `GSE185263.rds` | GSE185263 | Main evaluation |
| `GSE63042.rds` | GSE63042 | Main evaluation |
| `GSE236713.rds` | GSE236713 | Main evaluation (leukocyte RNA; deposited values are centred gene by gene) |
| `GSE272769_expr.rds`, `GSE272769_sample_status.csv` | GSE272769 | Replication |
| `GSE95233.rds` | GSE95233 | Sensitivity analyses |
| `GSE54514.rds`, `pdata_GSE54514.csv` | GSE54514 | Sensitivity analyses; measured neutrophil fractions for the accuracy check of the composition estimates |

Each `.rds` file of the first six cohorts is a list with `expr` (genes × samples, log2 scale) and `pheno` (sample sheet with group and death). `GSE272769_expr.rds` holds the expression matrix only; its sample sheet is the CSV file. Processing of each cohort is described in the Supplementary Methods of the article.

## Signatures

- `model_table.csv`: one row per model. The 36 signatures evaluated in the article are those with `tier` equal to `A_coef_complete` that do not belong to the ten reports of reference models (the PubMed identifiers are listed in `scripts/_common.R`, `REF_PMIDS`). The columns named after cohorts give the eligibility of each signature in each cohort: a signature is not evaluated in a cohort that its original report used in any role.
- `genes.csv`: genes, coefficients and direction of each model, including the endotyping tools and reference scores used in Part III.
- `cohort_usage.csv`: the cohorts used by each original report and their role.
- `extraction_records/`: one JSON file per report (named by PubMed identifier) with the extracted formula, direction, training cohort and reported performance.
- Nine signatures have an inconsistency in the original report (for example, a coefficient sign that contradicts the text). They were used as printed and are flagged in `tables/TableS1_signature_key.csv`, which also maps the codes S01–S36 used in the article to the reports.

## Requirements

R 4.4.3 with data.table, jsonlite, limma, edgeR, org.Hs.eg.db, AnnotationDbi, glmnet, survival, matrixStats, metafor, ggplot2 and their dependencies. Script `05_batch_composition.R` additionally needs SepstratifieR (1.0), ConsensusTranscriptomicSubtype (0.1.4), sva and batchelor. Exact versions are recorded in `environment/`.

## How to run

Run every script from the repository root, for example:

```bash
Rscript scripts/03_composition.R
```

To run from another directory, set `SEPSIS_BENCH_ROOT` to the repository root.

| Order | Script | Purpose | Approximate time |
|---|---|---|---|
| — | `_common.R` | Paths, data loading, scoring and AUROC functions shared by the other scripts | — |
| before the analysis | `p0_precision_simulation.R` | Precision simulation used to set the margins | minutes |
| 1 | `03_composition.R` | Leukocyte composition estimates and quality gate; no outcomes used (needs LM22) | < 1 min |
| 2 | `05_batch_composition.R` | Part III: batch-composition scenarios and hypothesis H4; no outcomes used | about 1 min |
| 3 | `01_build_random_models.R` | Builds and archives all random signatures in the training cohort | about 5 h on 10 cores |
| 4 | `02_evaluate.R` | Part I: evaluation, bootstrap, hypotheses H1 and H2, sensitivity analyses | about 2 h |
| 5 | `04_composition_adjusted.R` | Part II: hypothesis H3 and composition-adjusted results | about 4 min |
| 6 | `06_holm.R` | Holm adjustment of H2–H4 and final decisions | seconds |
| 7 | `08_cmaise_validation.R` | Local validation in CMAISE (needs the controlled data and LM22) | about 18 min |
| 8 | `10_prespecified_descriptive.R` | Three pre-specified descriptive analyses, computed from archived outputs | about 2 min |
| 9 | `11_per_cohort.R` | Post hoc analysis: each main cohort alone | about 4 min |
| 10 | `07_figures.R`, `07b_fig1_design.R`, `07c_graphical_abstract.R` | Figures and source tables | minutes |
| 11 | `09_manuscript_numbers.R` | Collects the numbers cited in the article into `tables/numbers/` | seconds |

Scripts write to `results/`, `tables/` and `figures/` and overwrite the archived files. Setting `DRY_RUN=1` permutes all outcomes and writes to `results_dryrun/` instead; the output is meaningless and serves only to check that the code runs.

Step 3 can be skipped: `02_evaluate.R` and the later scripts read the archived random signatures in `results/part1/random_models.rds` and verify the file against `results/part1/random_models_checksum.txt`.

The base seed is 20261116. Each attempt of the random-signature pipeline uses a seed determined by the signature and the attempt number, so the archive does not depend on parallel execution.

## Archived random signatures

`results/part1/random_models.rds` is an R list:

| Element | Content |
|---|---|
| `U` | Gene universe (9,471 genes measured in the training, main and replication cohorts) |
| `deg` | Genes differentially expressed between patients and healthy controls in the training cohort (1,478) |
| `R2` | 7,200 random-theme pipeline signatures, 200 per published signature (`matched`). Each has `genes` and `coef` (multivariable Cox weights), the LASSO and stepwise variants (`lasso_genes`, `lasso_coef`, `step_genes`, `step_coef`), the seed, the number of candidate genes and the training sample size |
| `R2b` | 1,800 signatures from the same pipeline without the differential expression step |
| `R1` | 12,999 multivariable Cox models on random genes without selection, about 1,000 per gene count |
| `R3` | 200 runs of the pipeline with the whole gene universe as theme |
| `params` | The pipeline parameters used |

## Results

| Path | Content |
|---|---|
| `results/hypotheses_final.csv` | Final decisions on the four hypotheses, with Holm adjustment |
| `results/part1/` | `part1_estimates.csv` (main analysis, replication and all sensitivity analyses), `H1_H2_result.csv`, `E3_percentiles*.csv` (percentile of each published signature among its matched random signatures), `per_cohort_*.csv` (post hoc), `evaluated_cohorts.csv`, `random_scores_main*.rds` (scores of the random signatures in the main cohorts), `bootstrap_*.rds` |
| `results/part2/` | `H3_result.csv`, `H3_estimates.csv`, `score_vs_composition.csv`, `tnlr_auroc.csv`, `random_models_adjusted.csv` |
| `results/composition/` | Composition estimates of each public cohort and the quality gate in GSE54514 |
| `results/batch/` | `batch_comparisons.csv`, `H4_result.csv`, `batch_outputs_long.csv` (label and score of every sample in every scenario), agreement between tools and the association between class and death by scenario |
| `results/cmaise/` | Aggregate results of the local validation |
| `results/power/` | Outputs of the precision simulation |
