source("scripts/_common.R")
suppressPackageStartupMessages(library(matrixStats))
set.seed(SEED + 11)
IN <- file.path(PROJ, "results"); od <- outdir("part1")
B <- if (DRY_RUN) 20 else 2000
M <- load_models(); pub_ids <- M$mt$model_id
RSM <- readRDS(file.path(IN, "part1", "random_scores_main_matrix.rds"))
RM <- readRDS(file.path(IN, "part1", "random_models.rds")); r2 <- data.table(rkey = sapply(RM$R2, `[[`, "key"), matched = sapply(RM$R2, `[[`, "matched")); rm(RM)
ev <- fread(file.path(IN, "part1", "evaluated_cohorts.csv"))
stopifnot(identical(names(RSM), MAIN), all(sapply(RSM, function(m) identical(rownames(m), r2$rkey))))
rowAUC <- function(S, y) { n1 <- sum(y == 1); n0 <- sum(y == 0); r <- rowRanks(S, ties.method = "average"); (rowSums(r[, y == 1, drop = FALSE]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
rowDelong <- function(S, y) {
  pos <- which(y == 1); neg <- which(y == 0); n1 <- length(pos); n0 <- length(neg)
  Ra <- rowRanks(S[, c(pos, neg), drop = FALSE], ties.method = "average")
  V10 <- (Ra[, seq_len(n1), drop = FALSE] - rowRanks(S[, pos, drop = FALSE], ties.method = "average")) / n0
  V01 <- 1 - (Ra[, n1 + seq_len(n0), drop = FALSE] - rowRanks(S[, neg, drop = FALSE], ties.method = "average")) / n1
  cbind(auc = rowMeans(V10), se = sqrt(rowVars(V10) / n1 + rowVars(V01) / n0)) }

res <- list(); per_sig <- list(); DL <- list()
for (g in MAIN) {
  x <- load_cohort(g); y <- x$ph$death; stopifnot(ncol(RSM[[g]]) == length(y), !anyNA(y))
  ids <- ev[cohort == g & evaluated == TRUE]$model_id
  P <- t(sapply(ids, function(id) score_published(M, id, x$expr))); rownames(P) <- ids; stopifnot(!anyNA(P))
  keep <- r2$matched %in% ids; Rg <- RSM[[g]][keep, , drop = FALSE]
  st <- function(i) { a <- rowAUC(P[, i, drop = FALSE], y[i]); r <- rowAUC(Rg[, i, drop = FALSE], y[i]); c(pub_median = median(a), E1 = median(r), Delta = median(a) - median(r)) }
  s0 <- st(seq_along(y)); bt <- replicate(B, st(strat_idx(list(y))[[1]]))
  res[[g]] <- data.table(cohort = g, n = length(y), deaths = sum(y), n_published = length(ids), n_random = sum(keep), stat = names(s0), estimate = s0,
                         lo = apply(bt, 1, quantile, 0.025), hi = apply(bt, 1, quantile, 0.975))
  a0 <- rowAUC(P, y); r0 <- rowAUC(Rg, y); mk <- r2$matched[keep]
  per_sig[[g]] <- data.table(cohort = g, model_id = ids, auc = a0, random_median = sapply(ids, function(id) median(r0[mk == id])), percentile = sapply(ids, function(id) mean(r0[mk == id] < a0[[id]])))
  DL[[g]] <- list(pub = rowDelong(P, y), ids = ids, r2 = rowDelong(RSM[[g]], y))
  message(g, ": ", length(ids), " published signatures, ", sum(keep), " random signatures")
}
res <- rbindlist(res); per_sig <- rbindlist(per_sig)
fwrite(res, file.path(od, "per_cohort_estimates.csv")); fwrite(per_sig, file.path(od, "per_cohort_per_signature.csv")); print(res)
sm <- per_sig[, .(n_published = .N, n_below_random_median = sum(auc < random_median), n_gt95 = sum(percentile > 0.95), n_auc_lt050 = sum(auc < 0.5)), by = cohort]
fwrite(sm, file.path(od, "per_cohort_summary.csv")); print(sm)

if (!DRY_RUN) {
  e3 <- fread(file.path(IN, "part1", "E3_percentiles.csv")); p0 <- readRDS(file.path(IN, "part1", "bootstrap_main.rds"))$st0$pooled
  chk <- sapply(pub_ids, function(id) { gs <- MAIN[sapply(MAIN, function(g) id %in% DL[[g]]$ids)]; a <- sapply(gs, function(g) DL[[g]]$pub[match(id, DL[[g]]$ids), "auc"]); s <- sapply(gs, function(g) DL[[g]]$pub[match(id, DL[[g]]$ids), "se"]); pool_fixed(a, s)[["auc"]] })
  d1 <- max(abs(chk - e3$pooled_auc[match(pub_ids, e3$model_id)]))
  chk2 <- sapply(seq_len(nrow(r2)), function(j) { gs <- ev[model_id == r2$matched[j] & cohort %in% MAIN & evaluated == TRUE]$cohort; pool_fixed(sapply(gs, function(g) DL[[g]]$r2[j, "auc"]), sapply(gs, function(g) DL[[g]]$r2[j, "se"]))[["auc"]] })
  d2 <- max(abs(chk2 - p0[r2$rkey]))
  message(sprintf("check: max |pooled from per-cohort − archived pooled|: published %.2e, random %.2e", d1, d2)); stopifnot(d1 < 1e-9, d2 < 1e-9)
  writeLines(sprintf("published %.3e\nrandom %.3e", d1, d2), file.path(od, "per_cohort_check.txt"))
}
writeLines(capture.output(sessionInfo()), file.path(od, "sessionInfo_11.txt"))
message("done")
