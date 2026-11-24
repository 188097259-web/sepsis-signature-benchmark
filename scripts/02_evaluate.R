source("scripts/_common.R")
suppressPackageStartupMessages({ library(matrixStats); library(metafor) })
set.seed(SEED + 2)
od <- outdir("part1")
P <- fromJSON(file.path(PROJ, "params/pipeline_params.json")); DELTA <- P$h1_margin; B <- P$bootstrap_B
if (DRY_RUN) B <- 20
rf <- file.path(od, "random_models.rds")
ck <- readLines(file.path(od, "random_models_checksum.txt"))[1]; stopifnot(grepl(md5_of(rf), ck, fixed = TRUE))
RM <- readRDS(rf); M <- load_models()
VAR <- P$main_weights
pick <- function(r, v) switch(v, multicox = list(genes = r$genes, coef = r$coef), lasso = list(genes = r$lasso_genes, coef = r$lasso_coef),
                             stepwise = list(genes = r$step_genes, coef = r$step_coef))

COH <- c(MAIN, REPL, SENS)
dat <- lapply(setNames(COH, COH), function(g) { x <- load_cohort(g); list(expr = x$expr, y = x$ph$death) })
Y <- lapply(dat, `[[`, "y")

pub_ids <- M$mt$model_id
rows <- list()
for (id in pub_ids) rows[[paste0("PUB_", id)]] <- list(type = "PUB", matched = id, f = function(e, id) score_published(M, id, e), arg = id)
mk_lin <- function(genes, coef) list(genes = genes, coef = coef)
for (r in RM$R2) { w <- pick(r, VAR); rows[[r$key]] <- list(type = "R2", matched = r$matched, lin = mk_lin(w$genes, w$coef)) }
for (v in setdiff(c("lasso", "stepwise", "multicox"), VAR)) for (r in RM$R2) { w <- pick(r, v); if (length(w$genes) >= 1)
  rows[[paste0(r$key, "_", v)]] <- list(type = paste0("R2_", v), matched = r$matched, lin = mk_lin(w$genes, w$coef)) }
if (!is.null(RM$R2b)) for (r in RM$R2b) { w <- pick(r, VAR); rows[[r$key]] <- list(type = "R2_nodeg", matched = r$matched, lin = mk_lin(w$genes, w$coef)) }
for (r in RM$R1) rows[[r$key]] <- list(type = "R1", k = r$k, lin = mk_lin(r$genes, r$coef))
for (r in RM$R3) { w <- pick(r, VAR); rows[[r$key]] <- list(type = "R3", lin = mk_lin(w$genes, w$coef)) }
meta <- data.table(key_ = names(rows), type = sapply(rows, `[[`, "type"),
                   matched = sapply(rows, function(r) if (is.null(r$matched)) NA_character_ else r$matched),
                   k = sapply(rows, function(r) if (is.null(r$k)) NA_integer_ else as.integer(r$k)))
setnames(meta, "key_", "key")
score_mat <- function(g) {
  e <- dat[[g]]$expr
  S <- t(sapply(rows, function(r) { s <- if (r$type == "PUB") r$f(e, r$arg) else lin_score(e, r$lin$genes, r$lin$coef, "center")
    if (is.null(s)) rep(NA_real_, ncol(e)) else s }))
  S
}
SM <- lapply(setNames(COH, COH), score_mat)
message("score matrices: ", nrow(SM[[1]]), " models")
saveRDS(lapply(setNames(MAIN, MAIN), function(g) SM[[g]][meta$type == "R2", , drop = FALSE]), file.path(od, "random_scores_main_matrix.rds"))
R2keys <- meta[type == "R2"]$key
saveRDS(setNames(lapply(R2keys, function(k) setNames(lapply(MAIN, function(g) SM[[g]][k, ]), MAIN)), R2keys), file.path(od, "random_scores_main.rds"))

rowDelong <- function(S, y) {
  pos <- which(y == 1); neg <- which(y == 0); n1 <- length(pos); n0 <- length(neg)
  Ra <- rowRanks(S[, c(pos, neg), drop = FALSE], ties.method = "average")
  V10 <- (Ra[, seq_len(n1), drop = FALSE] - rowRanks(S[, pos, drop = FALSE], ties.method = "average")) / n0
  V01 <- 1 - (Ra[, n1 + seq_len(n0), drop = FALSE] - rowRanks(S[, neg, drop = FALSE], ties.method = "average")) / n1
  cbind(auc = rowMeans(V10), se = sqrt(rowVars(V10) / n1 + rowVars(V01) / n0))
}
elig_of <- function(id, cohorts) cohorts[sapply(cohorts, function(g) eligible(M, id, g))]
EL <- lapply(setNames(pub_ids, pub_ids), elig_of, cohorts = c(MAIN, REPL, SENS))
EV <- lapply(setNames(pub_ids, pub_ids), function(id) EL[[id]][sapply(EL[[id]], function(g) !anyNA(SM[[g]][paste0("PUB_", id), ]))])
ev_tab <- rbindlist(lapply(pub_ids, function(id) data.table(model_id = id, cohort = c(MAIN, REPL, SENS),
  eligible = c(MAIN, REPL, SENS) %in% EL[[id]], evaluated = c(MAIN, REPL, SENS) %in% EV[[id]])))
fwrite(ev_tab, file.path(od, "evaluated_cohorts.csv"))
message("eligible but not scorable (coverage < 80%): ", paste(ev_tab[eligible & !evaluated, paste(model_id, cohort)], collapse = "; "))
k_of <- sapply(pub_ids, function(id) nrow(M$gn[model_id == id][!duplicated(gene)]))
mask_for <- function(set) {
  mk <- matrix(FALSE, nrow(meta), length(set), dimnames = list(meta$key, set))
  for (j in seq_len(nrow(meta))) { ty <- meta$type[j]
    if (ty == "PUB" || grepl("^R2", ty)) mk[j, ] <- set %in% EV[[meta$matched[j]]] else if (ty == "R3") mk[j, ] <- TRUE }
  mk
}
MASKS <- list()
pool_mat <- function(AUC, SE, mk) {
  L <- qlogis(AUC); W <- 1 / (SE / (AUC * (1 - AUC)))^2
  bad <- !mk | !is.finite(L) | !is.finite(W) | is.na(W); W[bad] <- 0; L[bad] <- 0
  sw <- rowSums(W); out <- plogis(rowSums(W * L) / sw); out[sw == 0] <- NA; out
}
r1_rows <- which(meta$type == "R1")
stats_once <- function(idx, set = MAIN) {
  key <- paste(set, collapse = "+"); if (is.null(MASKS[[key]])) MASKS[[key]] <<- mask_for(set); mk <- MASKS[[key]]
  A <- lapply(set, function(g) rowDelong(SM[[g]][, idx[[g]], drop = FALSE], Y[[g]][idx[[g]]]))
  AUC <- sapply(A, function(a) a[, "auc"]); SE <- sapply(A, function(a) a[, "se"]); dim(AUC) <- dim(SE) <- c(nrow(meta), length(set))
  pooled <- setNames(pool_mat(AUC, SE, mk), meta$key)
  combo <- sapply(pub_ids, function(id) paste(k_of[[id]], paste(set %in% EV[[id]], collapse = ""), sep = "|"))
  cm <- sapply(unique(combo), function(cb) { id <- names(combo)[combo == cb][1]; rr <- r1_rows[meta$k[r1_rows] == k_of[[id]]]
    if (!length(rr) || !any(set %in% EV[[id]])) return(NA)
    m1 <- matrix(set %in% EV[[id]], length(rr), length(set), byrow = TRUE)
    median(pool_mat(AUC[rr, , drop = FALSE], SE[rr, , drop = FALSE], m1), na.rm = TRUE) })
  list(pooled = pooled, r1_med = setNames(cm[combo], pub_ids))
}
summ <- function(st, pub_subset = pub_ids) {
  p <- st$pooled; pub <- p[paste0("PUB_", pub_subset)]
  r2 <- p[meta[type == "R2" & matched %in% pub_subset]$key]
  out <- c(pub_median = median(pub, na.rm = TRUE), E1 = median(r2, na.rm = TRUE))
  out["Delta"] <- out["pub_median"] - out["E1"]
  for (v in c("R2_lasso", "R2_stepwise", "R2_multicox", "R2_nodeg")) { kk <- meta[type == v & matched %in% pub_subset]$key
    if (length(kk)) { m <- median(p[kk], na.rm = TRUE); out[paste0("E1_", v)] <- m; out[paste0("Delta_", v)] <- out["pub_median"] - m } }
  out["E1_R1"] <- median(st$r1_med[pub_subset], na.rm = TRUE); out["Delta_R1"] <- out["pub_median"] - out["E1_R1"]
  out["R3_median"] <- median(p[meta[type == "R3"]$key], na.rm = TRUE)
  out
}
SUBSETS <- list(all = pub_ids, trained_GSE65682 = P$trained_on_GSE65682, exclude_flagged = setdiff(pub_ids, FLAGGED))

run_set <- function(set, label) {
  Ys <- Y[set]; st0 <- stats_once(full_idx(Ys), set)
  est <- sapply(SUBSETS, function(s) summ(st0, s))
  bt <- replicate(B, { st <- stats_once(strat_idx(Ys), set); sapply(SUBSETS, function(s) summ(st, s)) }, simplify = "array")
  saveRDS(list(st0 = st0, bt = bt), file.path(od, paste0("bootstrap_", label, ".rds")))
  rbindlist(lapply(colnames(est), function(s) data.table(analysis = label, subset = s, stat = rownames(est), estimate = est[, s],
    lo = apply(bt[, s, , drop = FALSE], 1, quantile, 0.025, na.rm = TRUE), hi = apply(bt[, s, , drop = FALSE], 1, quantile, 0.975, na.rm = TRUE))))
}
res <- rbind(run_set(MAIN, "main"), run_set(REPL, "replication_GSE272769"), run_set(c(MAIN, SENS), "S4_main_plus_sensitivity_cohorts"))
fwrite(res, file.path(od, "part1_estimates.csv"))

m <- res[analysis == "main" & subset == "all"]
H <- data.table(hypothesis = c("H1", "H2"),
                estimate = c(m[stat == "Delta"]$estimate, m[stat == "E1"]$estimate),
                lo = c(m[stat == "Delta"]$lo, m[stat == "E1"]$lo), hi = c(m[stat == "Delta"]$hi, m[stat == "E1"]$hi),
                rule = c(sprintf("upper 95%% CI of Delta < %.2f", DELTA), "lower 95% CI of E1 > 0.50"))
H[, supported := c(hi[1] < DELTA, lo[2] > 0.5)]
H[hypothesis == "H1", note := fifelse(lo > DELTA, "published models better than random by more than margin", fifelse(supported, "non-superiority shown", "inconclusive"))]
rep_d <- res[analysis == "replication_GSE272769" & subset == "all" & stat == "Delta"]
H <- rbind(H, data.table(hypothesis = "replication_consistency", estimate = rep_d$estimate, lo = rep_d$lo, hi = rep_d$hi,
  rule = "same sign as main Delta and |difference| < 0.05", supported = sign(rep_d$estimate) == sign(H$estimate[1]) & abs(rep_d$estimate - H$estimate[1]) < 0.05, note = ""), fill = TRUE)
fwrite(H, file.path(od, "H1_H2_result.csv")); print(H)

st0 <- readRDS(file.path(od, "bootstrap_main.rds"))$st0; p <- st0$pooled
E3 <- rbindlist(lapply(pub_ids, function(id) { a <- p[paste0("PUB_", id)]; r <- p[meta[type == "R2" & matched == id]$key]
  data.table(model_id = id, pooled_auc = a, random_median = median(r, na.rm = TRUE), percentile = mean(r < a, na.rm = TRUE), n_random = sum(!is.na(r))) }))
fwrite(E3, file.path(od, "E3_percentiles.csv"))
message(sprintf("E3: %d of %d published models above the 95th percentile of their matched random models", sum(E3$percentile > 0.95, na.rm = TRUE), nrow(E3)))

A0 <- lapply(setNames(MAIN, MAIN), function(g) { r <- rowDelong(SM[[g]], Y[[g]]); rownames(r) <- meta$key; r })
re_pool <- function(k, gs) { a <- sapply(gs, function(g) A0[[g]][k, "auc"]); s <- sapply(gs, function(g) A0[[g]][k, "se"]); if (length(gs) < 2) return(NA)
  f <- tryCatch(rma(yi = qlogis(a), sei = s / (a * (1 - a)), method = "REML"), error = function(e) NULL); if (is.null(f)) NA else plogis(f$b[1]) }
pubre <- sapply(pub_ids, function(id) re_pool(paste0("PUB_", id), intersect(EV[[id]], MAIN)))
r2m <- meta[type == "R2"]
r2re <- mapply(function(kk, mm) re_pool(kk, intersect(EV[[mm]], MAIN)), r2m$key, r2m$matched)
fwrite(data.table(stat = c("pub_median_RE", "E1_RE", "Delta_RE"), estimate = c(median(pubre, na.rm = TRUE), median(r2re, na.rm = TRUE), median(pubre, na.rm = TRUE) - median(r2re, na.rm = TRUE))),
       file.path(od, "S5_random_effects.csv"))
writeLines(capture.output(sessionInfo()), file.path(od, "sessionInfo_02.txt"))
message("done")
