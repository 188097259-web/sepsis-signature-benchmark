source("scripts/_common.R")
set.seed(SEED)
od <- outdir("part2"); comp_dir <- file.path(OUT, "composition")
B <- if (DRY_RUN) 20 else 2000
COVARS <- c("neutrophil", "mono_mac_dc", "t_cell", "nk", "b_plasma")
MK_COVARS <- c("mk_neutrophil", "mk_mono_mac_dc", "mk_t_cell", "mk_nk", "mk_b_plasma")
M <- load_models()
gate <- fread(file.path(comp_dir, "quality_gate_GSE54514.csv"))
message("quality gate pass: ", gate$pass)

COH <- c(MAIN, REPL)
dat <- lapply(setNames(COH, COH), function(g) {
  x <- load_cohort(g); cp <- fread(file.path(comp_dir, paste0("composition_", g, ".csv")))[match(x$ph$sample, sample)]
  list(expr = x$expr, y = x$ph$death, comp = cp, sofa = if ("sofa" %in% names(x$ph)) suppressWarnings(as.numeric(x$ph$sofa)) else NULL)
})
resid_on <- function(s, C) { ok <- complete.cases(C) & !is.na(s); r <- rep(NA_real_, length(s)); r[ok] <- resid(lm(s[ok] ~ ., data = as.data.frame(C[ok, , drop = FALSE]))); r }

S <- list()
for (id in M$mt$model_id) for (g in COH) if (eligible(M, id, g)) {
  s <- score_published(M, id, dat[[g]]$expr); if (is.null(s)) next
  adj <- resid_on(s, dat[[g]]$comp[, ..COVARS]); adj_mk <- resid_on(s, dat[[g]]$comp[, ..MK_COVARS])
  S[[id]][[g]] <- list(raw = s, adj = adj, adj_mk = adj_mk, adj_s = if (g %in% NNLS_INVALID) adj_mk else adj)
}
desc <- rbindlist(lapply(names(S), function(id) rbindlist(lapply(names(S[[id]]), function(g) { s <- S[[id]][[g]]$raw; C <- dat[[g]]$comp
  data.table(model_id = id, cohort = g, nnls_valid = !(g %in% NNLS_INVALID),
             r2_comp = summary(lm(s ~ ., data = as.data.frame(C[, ..COVARS])))$r.squared,
             r2_marker = summary(lm(s ~ ., data = as.data.frame(C[, ..MK_COVARS])))$r.squared,
             rho_tnlr = if (g %in% NNLS_INVALID) NA_real_ else cor(s, C$tnlr, method = "spearman")) }))))
fwrite(desc, file.path(od, "score_vs_composition.csv"))

r2_once <- function(idx, cohorts = PART2_MAIN, marker_for = character(0)) {
  per <- sapply(names(S), function(id) { gs <- intersect(names(S[[id]]), cohorts); if (!length(gs)) return(NA)
    r2 <- sapply(gs, function(g) { ii <- idx[[g]]; cv <- if (g %in% marker_for) MK_COVARS else COVARS
      summary(lm(S[[id]][[g]]$raw[ii] ~ ., data = as.data.frame(dat[[g]]$comp[ii, ..cv])))$adj.r.squared })
    w <- sapply(gs, function(g) length(idx[[g]])); sum(r2 * w) / sum(w) })
  median(per, na.rm = TRUE)
}
nidx <- lapply(dat[PART2_MAIN], function(d) seq_along(d$y))
r2_0 <- r2_once(nidx); r2_bt <- replicate(B, r2_once(lapply(nidx, function(i) sample(i, replace = TRUE))))
lo_b <- 2 * r2_0 - unname(quantile(r2_bt, 0.975)); hi_b <- 2 * r2_0 - unname(quantile(r2_bt, 0.025))
p3 <- (1 + sum(r2_bt >= 2 * r2_0 - 0.5)) / (B + 1)
saveRDS(list(r2_0 = r2_0, r2_bt = r2_bt), file.path(od, "H3A_bootstrap.rds"))
h3a <- data.table(form = "A_adjR2", median_adj_r2 = r2_0, bias_corrected = 2 * r2_0 - mean(r2_bt), lo = lo_b, hi = hi_b, pct_lo = unname(quantile(r2_bt, 0.025)), pct_hi = unname(quantile(r2_bt, 0.975)), p_one_sided = p3, supported = lo_b > 0.5)
fwrite(h3a, file.path(od, "H3A_R2_result.csv")); print(h3a)
sidx <- lapply(dat[MAIN], function(d) seq_along(d$y))
s0 <- r2_once(sidx, MAIN, NNLS_INVALID); s_bt <- replicate(B, r2_once(lapply(sidx, function(i) sample(i, replace = TRUE)), MAIN, NNLS_INVALID))
fwrite(data.table(analysis = "main3_GSE236713_marker", median_adj_r2 = s0, bias_corrected = 2 * s0 - mean(s_bt),
                  lo = 2 * s0 - unname(quantile(s_bt, 0.975)), hi = 2 * s0 - unname(quantile(s_bt, 0.025))),
       file.path(od, "H3A_sensitivity_GSE236713_marker.csv"))

pooled_set <- function(which, cohorts, idx) {
  sapply(names(S), function(id) { gs <- intersect(names(S[[id]]), cohorts); if (!length(gs)) return(NA)
    r <- sapply(gs, function(g) { ok <- !is.na(S[[id]][[g]][[which]]); ii <- idx[[g]][ok[idx[[g]]]]; delong(S[[id]][[g]][[which]][ii], dat[[g]]$y[ii]) })
    pool_fixed(r[1, ], r[2, ])["auc"] })
}
est <- function(cohorts, idx, adj_field = "adj") {
  raw <- pooled_set("raw", cohorts, idx); adj <- pooled_set(adj_field, cohorts, idx); mk <- pooled_set("adj_mk", cohorts, idx)
  mr <- median(raw, na.rm = TRUE); ma <- median(adj, na.rm = TRUE); mm <- median(mk, na.rm = TRUE)
  c(median_raw = mr, median_adj = ma, rf = (ma - 0.5) / (mr - 0.5), diff = mr - ma, half_excess_margin = (mr - ma) - 0.5 * (mr - 0.5),
    median_adj_marker = mm, rf_marker = (mm - 0.5) / (mr - 0.5))
}
Y <- lapply(dat, `[[`, "y")
res <- list()
for (set in list(list(name = "main", co = PART2_MAIN, af = "adj"), list(name = "replication", co = REPL, af = "adj"),
                 list(name = "main3_GSE236713_marker", co = MAIN, af = "adj_s"))) {
  Yc <- Y[set$co]; e0 <- est(set$co, full_idx(Yc), set$af)
  bt <- replicate(B, est(set$co, strat_idx(Yc), set$af))
  res[[set$name]] <- data.table(analysis = set$name, stat = names(e0), estimate = e0,
                                lo = apply(bt, 1, quantile, 0.025, na.rm = TRUE), hi = apply(bt, 1, quantile, 0.975, na.rm = TRUE))
}
res <- rbindlist(res); fwrite(res, file.path(od, "H3_estimates.csv"))
h3 <- copy(h3a)
h3[, gate_pass := gate$pass]; fwrite(h3, file.path(od, "H3_result.csv")); print(h3)

tn <- rbindlist(lapply(setdiff(COH, NNLS_INVALID), function(g) { a <- delong(dat[[g]]$comp$tnlr, Y[[g]]); data.table(cohort = g, auc = a[1], se = a[2]) }))
tn_pool <- tn[cohort %in% PART2_MAIN, as.list(pool_fixed(auc, se))]
fwrite(rbind(tn[, .(cohort, auc, lo = auc - 1.96 * se, hi = auc + 1.96 * se)], data.table(cohort = "pooled_main", tn_pool)), file.path(od, "tnlr_auroc.csv"))
d <- dat[["GSE185263"]]; ok <- !is.na(d$sofa)
sofa <- data.table(score = c("SOFA", "SOFA_adjusted"), rbind(delong(d$sofa[ok], d$y[ok]), delong(resid_on(d$sofa, d$comp[, ..COVARS])[ok], d$y[ok])))
fwrite(sofa, file.path(od, "sofa_adjusted_GSE185263.csv"))

rf <- file.path(OUT, "part1", "random_scores_main.rds")
if (file.exists(rf)) {
  R <- readRDS(rf)
  RMm <- readRDS(file.path(OUT, "part1", "random_models.rds")); matched_of <- setNames(sapply(RMm$R2, `[[`, "matched"), sapply(RMm$R2, `[[`, "key")); rm(RMm)
  stopifnot(all(names(R) %in% names(matched_of)))
  rr <- rbindlist(lapply(names(R), function(k) { gs <- intersect(intersect(names(R[[k]]), PART2_MAIN), names(S[[matched_of[[k]]]])); if (!length(gs)) return(NULL)
    a0 <- sapply(gs, function(g) delong(R[[k]][[g]], Y[[g]])); a1 <- sapply(gs, function(g) delong(resid_on(R[[k]][[g]], dat[[g]]$comp[, ..COVARS]), Y[[g]]))
    data.table(model = k, raw = pool_fixed(a0[1, ], a0[2, ])["auc"], adj = pool_fixed(a1[1, ], a1[2, ])["auc"]) }))
  fwrite(rr, file.path(od, "random_models_adjusted.csv"))
}
writeLines(capture.output(sessionInfo()), file.path(od, "sessionInfo.txt"))
message("done")
