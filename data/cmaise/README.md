# CMAISE data (controlled access; not distributed)

Script `scripts/08_cmaise_validation.R` reads two files from this folder, or from the folder named by the environment
variable `CMAISE_DIR`:

- `counts.tsv.gz`: gene-level RNA-sequencing counts
- `clinical_samples.csv`: sample sheet with day, batch, centre, outcome and SOFA score

The script checks the SHA-256 of both files. Access conditions are given in the data availability statement of the article.
Only aggregate results of the CMAISE analysis are archived in `results/cmaise/`.
