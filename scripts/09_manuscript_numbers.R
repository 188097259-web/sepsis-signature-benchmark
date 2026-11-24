suppressPackageStartupMessages(library(data.table))
PROJ <- Sys.getenv("SEPSIS_BENCH_ROOT", "."); R <- file.path(PROJ, "results")
OUTD <- file.path(PROJ, "tables/numbers"); dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
N <- list(); add <- function(k, v, src) N[[length(N) + 1]] <<- data.table(item = k, value = as.character(v), source = src)
r3 <- function(x) sprintf("%.3f", x); pc <- function(x, d = 1) sprintf(paste0("%.", d, "f"), 100 * x)

rm_ <- readRDS(file.path(R, "part1/random_models.rds"))
ng <- sapply(rm_$R2, function(m) length(m$genes)); nc <- sapply(rm_$R2, function(m) m$n_cand)
add("R2_n", length(rm_$R2), "random_models.rds"); add("R2b_n", length(rm_$R2b), "random_models.rds"); add("R3_n", length(rm_$R3), "random_models.rds"); add("R1_n", length(rm_$R1), "random_models.rds")
add("R2_genes_median", median(ng), "random_models.rds"); add("R2_genes_range", paste(range(ng), collapse = "–"), "random_models.rds")
add("R2_candidates_median", median(nc), "random_models.rds"); add("R2_candidates_IQR", paste(quantile(nc, c(.25, .75)), collapse = "–"), "random_models.rds")
add("universe_genes", length(rm_$U), "random_models.rds"); add("DEG_n", length(rm_$deg), "random_models.rds")
add("R2_unique_genes", length(unique(unlist(lapply(rm_$R2, `[[`, "genes")))), "random_models.rds")
f2 <- fread(file.path(R, "part1/R2_pipeline_failures.csv"))
add("R2_tries_total", sum(f2$tries), "R2_pipeline_failures.csv"); add("R2_success_total", sum(f2$success), "R2_pipeline_failures.csv")
add("R2_failure_rate_overall_pct", pc(1 - sum(f2$success) / sum(f2$tries)), "R2_pipeline_failures.csv")
fr <- 1 - f2$success / f2$tries; add("R2_failure_rate_by_sig_range_pct", paste(pc(min(fr)), pc(max(fr)), sep = "–"), "R2_pipeline_failures.csv")
add("n_train_lt479", sum(f2$n_train < 479), "R2_pipeline_failures.csv")

e3 <- fread(file.path(R, "part1/E3_percentiles.csv")); c3 <- fread(file.path(R, "cmaise/cmaise_E3_percentiles.csv"))
add("E3_main_gt95", sum(e3$percentile > 0.95), "E3_percentiles.csv"); add("E3_main_lt05", sum(e3$percentile < 0.05), "E3_percentiles.csv")
add("E3_main_median_pct", r3(median(e3$percentile)), "E3_percentiles.csv")
add("E3_main_below_random_median", sum(e3$pooled_auc < e3$random_median), "E3_percentiles.csv")
add("pub_auc_range_main", paste(r3(range(e3$pooled_auc)), collapse = "–"), "E3_percentiles.csv")
add("pub_auc_ge070_main", sum(e3$pooled_auc >= 0.70), "E3_percentiles.csv")
pcol <- intersect(c("percentile", "pct"), names(c3))[1]
add("E3_cmaise_gt95", sum(c3[[pcol]] > 0.95, na.rm = TRUE), "cmaise_E3_percentiles.csv")
mm <- merge(e3[, .(model_id, pm = percentile)], c3[, .(model_id, pc_ = get(pcol))], by = "model_id")[!is.na(pc_)]
add("E3_both_gt95", sum(mm$pm > 0.95 & mm$pc_ > 0.95), "E3 merge"); add("E3_spearman_main_cmaise", sprintf("%.2f", cor(mm$pm, mm$pc_, method = "spearman")), "E3 merge")
add("E3_n_cmaise", nrow(mm), "E3 merge")
aucol <- intersect(c("pooled_auc", "auc", "pub"), names(c3))[1]
if (!is.na(aucol)) add("pub_auc_range_cmaise", paste(r3(range(c3[[aucol]], na.rm = TRUE)), collapse = "–"), "cmaise_E3_percentiles.csv")

qg <- fread(file.path(R, "composition/quality_gate_GSE54514.csv")); add("gate_rho_all", sprintf("%.3f", qg$rho_nnls), "quality_gate_GSE54514.csv")
add("gate_rho_marker", sprintf("%.3f", qg$rho_marker), "quality_gate_GSE54514.csv")
sv <- fread(file.path(R, "part2/score_vs_composition.csv"))[cohort %in% c("GSE185263", "GSE63042")]
NCOH <- c(GSE185263 = 345, GSE63042 = 106); adjr2 <- function(r2, n, p = 5) 1 - (1 - r2) * (n - 1) / (n - p - 1)
sv[, n := NCOH[cohort]][, `:=`(lm22 = adjr2(r2_comp, n), marker = adjr2(r2_marker, n))]
r2s <- sv[, .(lm22 = sum(lm22 * n) / sum(n), marker = sum(marker * n) / sum(n), rho = sum(rho_tnlr * n, na.rm = TRUE) / sum(n)), by = model_id]
add("R2_lm22_median", r3(median(r2s$lm22)), "score_vs_composition.csv"); add("R2_lm22_range", paste(r3(range(r2s$lm22)), collapse = "–"), "score_vs_composition.csv")
add("R2_lm22_n_gt05", sum(r2s$lm22 > 0.5), "score_vs_composition.csv"); add("R2_marker_median", r3(median(r2s$marker)), "score_vs_composition.csv")
add("R2_n_signatures", nrow(r2s), "score_vs_composition.csv")
add("rho_tnlr_median", sprintf("%.2f", median(r2s$rho)), "score_vs_composition.csv")
add("rho_tnlr_abs_median", sprintf("%.2f", median(abs(r2s$rho))), "score_vs_composition.csv")
rma <- fread(file.path(R, "part2/random_models_adjusted.csv"))
add("random_raw_median_part2", r3(median(rma$raw)), "random_models_adjusted.csv"); add("random_adj_median_part2", r3(median(rma$adj)), "random_models_adjusted.csv")
dc <- fread(file.path(R, "composition/deconvolution_summary.csv")); for (i in seq_len(nrow(dc))) add(paste0("neutro_median_", dc$cohort[i]), r3(dc$median_neutrophil[i]), "deconvolution_summary.csv")

bc <- fread(file.path(R, "batch/batch_comparisons.csv"))
TOOLS <- c("SRS_davenport", "SRS_extended", "SRS_davenport_k20", "CTS", "SUBSPACE_quadrant_healthyref", "SUBSPACE_quadrant_median")
hm <- bc[tool %in% TOOLS & !is.na(prop_changed) & !grepl("^seed|C4", contrast), .(cohort, tool, contrast, p = prop_changed, kappa)]
seedd <- bc[tool == "CTS" & grepl("^seed", contrast), .(p = median(prop_changed), kappa = median(kappa)), by = cohort][, `:=`(tool = "CTS", contrast = "seed")]
c5 <- rbindlist(lapply(list.files(file.path(R, "batch"), "^C5_subspace", full.names = TRUE), fread))[, .(p = median(change)), by = cohort][, `:=`(tool = "SUBSPACE_quadrant_healthyref", contrast = "C5", kappa = NA_real_)]
hm <- rbind(hm, seedd, c5, use.names = TRUE)
hs <- hm[, .(k = .N, median = median(p), min = min(p), max = max(p), kappa_median = median(kappa, na.rm = TRUE)), by = .(tool, contrast)]
fwrite(hs, file.path(OUTD, "batch_table.csv"))
srs <- bc[tool == "SRS_davenport" & contrast == "C0 vs C1", .(cohort, n, changed, kappa)]
for (i in seq_len(nrow(srs))) add(paste0("SRS7_C1_", srs$cohort[i]), sprintf("%d/%d (%.1f%%), kappa %.2f", srs$changed[i], srs$n[i], 100 * srs$changed[i] / srs$n[i], srs$kappa[i]), "batch_comparisons.csv")
h4 <- fread(file.path(R, "batch/H4_result.csv")); add("H4_I2", sprintf("%.1f", h4$I2), "H4_result.csv")
cts1 <- bc[tool == "CTS" & contrast == "C0 vs C1", .(cohort, n, changed)]
add("CTS_C1_total", sprintf("%d/%d", sum(cts1$changed), sum(cts1$n)), "batch_comparisons.csv")
add("SRS7_C1_total", sprintf("%d/%d", sum(srs$changed), sum(srs$n)), "batch_comparisons.csv")
ld <- fread(file.path(R, "batch/label_distribution.csv"))
hc <- ld[tool %in% c("CTS", "SRS_davenport", "SUBSPACE_quadrant_healthyref") & scenario %in% c("C4_healthy_alone", "C1"), .(n = sum(n)), by = .(tool, scenario, value)][, share := n / sum(n), by = .(tool, scenario)]
for (i in seq_len(nrow(hc))) add(paste("HC", hc$tool[i], hc$scenario[i], hc$value[i], sep = "_"), sprintf("%d (%.1f%%)", hc$n[i], 100 * hc$share[i]), "label_distribution.csv")
hk <- bc[contrast == "C4_healthy_alone vs C1" & tool %in% c("CTS", "SRS_davenport"), .(n = sum(n), changed = sum(changed)), by = tool]
for (i in seq_len(nrow(hk))) add(paste0("HC_relabel_", hk$tool[i]), sprintf("%d/%d (%.1f%%)", hk$changed[i], hk$n[i], 100 * hk$changed[i] / hk$n[i]), "batch_comparisons.csv")
CONT <- c("Mars1_ratio", "Sweeney2018NC_Duke18_raw", "Sweeney2018NC_Stanford12_raw", "SUBSPACE_myeloid_raw", "SUBSPACE_lymphoid_raw", "SRSq_davenport", "SRSq_extended",
          "Sweeney2018NC_Duke18_batchz", "Sweeney2018NC_Stanford12_batchz", "Jiang2022_IRRG9_batchz", "Chen2023_IRG7_batchz", "Feng2024_ROS37_batchz", "Li2026_GMGscore6_batchz")
cs <- bc[contrast == "C0 vs C1" & tool %in% CONT, .(k = .N, shift_median = median(median_abs_shift_sd), shift_max = max(median_abs_shift_sd), rank_shift_median = median(median_abs_rank_pct_shift), spearman_min = min(spearman)), by = tool]
fwrite(cs, file.path(OUTD, "continuous_shift_table.csv"))

cov <- fread(file.path(R, "cmaise/cmaise_signature_coverage.csv")); add("cmaise_cov_rows", nrow(cov), "cmaise_signature_coverage.csv")

e3r <- fread(file.path(R, "part1/E3_percentiles_replication_GSE272769.csv"))
add("E3_rep_gt95", sum(e3r$percentile > 0.95), "E3_percentiles_replication_GSE272769.csv"); add("E3_rep_lt05", sum(e3r$percentile < 0.05), "E3_percentiles_replication_GSE272769.csv")
add("E3_rep_below_random_median", sum(e3r$pooled_auc < e3r$random_median), "E3_percentiles_replication_GSE272769.csv")
add("E3_rep_median_pct", pc(median(e3r$percentile), 0), "E3_percentiles_replication_GSE272769.csv")
m3 <- merge(merge(e3[, .(model_id, pm = percentile)], e3r[, .(model_id, pr = percentile)], by = "model_id"), c3[, .(model_id, pc_ = get(pcol))], by = "model_id", all.x = TRUE)
add("E3_spearman_main_rep", sprintf("%.2f", cor(m3$pm, m3$pr, method = "spearman")), "E3 merge")
add("E3_gt95_in_two_or_more", sum((m3$pm > 0.95) + (m3$pr > 0.95) + ifelse(is.na(m3$pc_), 0, m3$pc_ > 0.95) >= 2), "E3 merge")
add("E3_gt95_distinct_signatures", sum(m3$pm > 0.95 | m3$pr > 0.95 | (!is.na(m3$pc_) & m3$pc_ > 0.95)), "E3 merge")
add("E3_rep_gt95_percentiles_main", paste(pc(sort(m3[pr > 0.95]$pm), 0), collapse = ", "), "E3 merge"); add("E3_rep_gt95_percentiles_cmaise", paste(pc(sort(m3[pr > 0.95]$pc_), 0), collapse = ", "), "E3 merge")
cb <- merge(bc[tool == "CTS" & contrast == "C0 vs C1", .(cohort, p = prop_changed)], seedd[, .(cohort, seed = p)], by = "cohort")[, net := p - seed]
add("CTS_C1_minus_seed_median_pct", pc(median(cb$net)), "batch_comparisons.csv"); add("CTS_C1_minus_seed_range_pct", paste(pc(min(cb$net)), pc(max(cb$net)), sep = "–"), "batch_comparisons.csv")
h3e <- fread(file.path(R, "part2/H3_estimates.csv")); rfm <- h3e[analysis == "main" & stat == "rf"]
add("RF_main", sprintf("%.2f (%.2f–%.2f)", rfm$estimate, rfm$lo, rfm$hi), "H3_estimates.csv")
add("rho_tnlr_IQR", paste(sprintf("%.2f", quantile(r2s$rho, c(.25, .75))), collapse = " to "), "score_vs_composition.csv")
ar <- fread(file.path(R, "batch/ari_pairwise_summary.csv"))
hr <- ar[scenario == "C1" & (tool_a == "SUBSPACE_quadrant_healthyref" | tool_b == "SUBSPACE_quadrant_healthyref")]
add("ARI_C1_healthyref_median_range", paste(sprintf("%.2f", range(hr$ari_median)), collapse = "–"), "ari_pairwise_summary.csv")
la <- fread(file.path(R, "batch/label_death_association.csv"))
add("V_change_C0C1_n_combinations", nrow(la[contrast == "C0 vs C1"]), "label_death_association.csv")
add("V_change_C0C1_n_ci_excluding_0", sum(la[contrast == "C0 vs C1", diff_lo > 0 | diff_hi < 0]), "label_death_association.csv")
add("V_change_C0C2_n_ci_excluding_0", sum(la[contrast == "C0 vs C2", diff_lo > 0 | diff_hi < 0]), "label_death_association.csv")
add("V_change_C0C1_n_range_small_datasets", paste(range(la[contrast == "C0 vs C1" & cohort != "GSE65682"]$n), collapse = "–"), "label_death_association.csv")

wb <- srs[cohort != "GSE236713"]
add("SRS7_C1_wholeblood", sprintf("%d/%d (%.1f%%)", sum(wb$changed), sum(wb$n), 100 * sum(wb$changed) / sum(wb$n)), "batch_comparisons.csv")
add("SRS7_C1_GSE236713_share", sprintf("%d/%d", srs[cohort == "GSE236713"]$changed, sum(srs$changed)), "batch_comparisons.csv")
add("CTS_seed_baseline_range_pct", paste(pc(min(seedd$p)), pc(max(seedd$p)), sep = "–"), "batch_comparisons.csv")
r3g <- sapply(rm_$R3, function(m) paste(sort(m$genes), collapse = "|"))
add("R3_distinct_signatures", length(unique(r3g)), "random_models.rds"); add("R3_genes_range", paste(range(sapply(rm_$R3, function(m) length(m$genes))), collapse = "–"), "random_models.rds")
cs_ <- fread(file.path(R, "cmaise/cmaise_sensitivity.csv"))[analysis == "shock_and_refractory_only" & stat == "Delta"]
add("cmaise_shock_share", sprintf("%d/355 (%.0f%%)", cs_$n, 100 * cs_$n / 355), "cmaise_sensitivity.csv")
pe_ <- fread(file.path(PROJ, "docs/pipeline_extraction.csv"))
add("orig_weights_counts", paste(sprintf("%s %d", c("multivariable_cox", "lasso_coefficients", "stepwise_cox", "logistic", "unweighted"),
    sapply(c("multivariable_cox", "lasso_coefficients", "stepwise_cox", "logistic", "unweighted"), function(k) sum(pe_$final_weights_from == k))), collapse = "; "), "pipeline_extraction.csv")
add("orig_weights_unclear", sum(grepl("^other", pe_$final_weights_from)), "pipeline_extraction.csv")
add("orig_deg_step", sprintf("yes %d; no %d; unclear %d", sum(pe_$deg_step %in% TRUE), sum(pe_$deg_step %in% FALSE), sum(is.na(pe_$deg_step))), "pipeline_extraction.csv")
add("orig_no_named_category", sum(grepl("^none", pe_$theme) | is.na(pe_$theme) | pe_$theme == ""), "pipeline_extraction.csv")
add("pub_auc_lt050_main", sum(e3$pooled_auc < 0.5), "E3_percentiles.csv")
x185 <- as.data.table(readRDS(file.path(PROJ, "data/cohorts/GSE185263.rds"))$pheno)[group == "sepsis"]
add("sofa_n_GSE185263", sprintf("%d/%d", sum(!is.na(suppressWarnings(as.numeric(x185$sofa)))), nrow(x185)), "GSE185263 pheno")
rmS <- readRDS(file.path(R, "part1/random_models.rds")); add("R2_stepwise_n", sprintf("%d/%d", sum(sapply(rmS$R2, function(m) length(m$step_genes)) >= 1), length(rmS$R2)), "random_models.rds"); rm(rmS)
add("cmaise_random_scored", sprintf("%d/%d", as.integer(sum(c3$n_random, na.rm = TRUE)), as.integer(200L * sum(!is.na(c3[[pcol]])))), "cmaise_E3_percentiles.csv")
svr <- fread(file.path(R, "part2/score_vs_composition.csv"))[cohort == "GSE272769"]
add("R2_lm22_median_GSE272769", r3(median(adjr2(svr$r2_comp, 161))), "score_vs_composition.csv"); add("R2_n_signatures_GSE272769", nrow(svr), "score_vs_composition.csv")
add("V_change_C0C1_n_relabelled_combinations", nrow(la[contrast == "C0 vs C1" & changed > 0]), "label_death_association.csv")
add("V_change_C0C1_relabelled_n_ci_excluding_0", sum(la[contrast == "C0 vs C1" & changed > 0, diff_lo > 0 | diff_hi < 0]), "label_death_association.csv")
add("V_change_CTS_n_positive", sprintf("%d/%d", sum(la[contrast == "C0 vs C1" & tool == "CTS"]$diff > 0), nrow(la[contrast == "C0 vs C1" & tool == "CTS"])), "label_death_association.csv")

out <- rbindlist(N); fwrite(out, file.path(OUTD, "manuscript_numbers.csv"))
print(out, nrows = 300); print(hs, nrows = 100); print(cs)
