source("scripts/_common.R")
suppressPackageStartupMessages({ library(edgeR); library(matrixStats); library(parallel) })
od <- outdir("cmaise")
P <- fromJSON(file.path(PROJ, "params/pipeline_params.json")); B <- if (DRY_RUN) 20 else P$bootstrap_B
NCORES <- as.integer(Sys.getenv("NCORES", "8"))
log_lines <- character(0); note <- function(...) { s <- paste0(...); message(s); log_lines <<- c(log_lines, s) }

CM_DIR <- CMAISE_DIR
CM_SHA <- c("counts.tsv.gz" = "9ad8d17b644cbaaa8a5fbead14462a3ea647df80f204006e50454afa1e4127b5",
            "clinical_samples.csv" = "dd178168d304ddeeeef3a535d7b8466db3cb7c198732568b6af70bc06019138d")
sha256 <- function(f) sub(" .*$", "", system2("shasum", c("-a", "256", shQuote(f)), stdout = TRUE))
for (f in names(CM_SHA)) stopifnot(sha256(file.path(CM_DIR, f)) == CM_SHA[[f]])
rf <- file.path(PROJ, "results/part1/random_models.rds")
stopifnot(md5_of(rf) == "77f8f686a062b8fb10bb07385694ee81")

cl <- fread(file.path(CM_DIR, "clinical_samples.csv"))
d1 <- cl[day == "d1"]; stopifnot(!anyDuplicated(d1$patient))
d1[, death := as.integer(mort_flg == "Y")][, death28 := as.integer(mort_flg == "Y" & Hospital_days <= 28)]
cnt <- fread(file.path(CM_DIR, "counts.tsv.gz"))
ens <- sub("\\..*$", "", cnt$gene_id); cm <- as.matrix(cnt[, d1$SampleName, with = FALSE])
mito <- grepl("^MT-", cnt$symbol) | cnt$symbol %in% c("MTND2P28", "MTCO1P12", "MTATP6P1"); rm(cnt)
mt_cnt <- colSums(cm[mito, , drop = FALSE])
note("mitochondrial rows removed: ", sum(mito), "; median mito counts per sample A = ", median(mt_cnt[d1$batch == "A"]), ", B = ", median(mt_cnt[d1$batch == "B"]))
cm <- cm[!mito, ]; ens <- ens[!mito]
sym <- suppressMessages(mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
sym <- official(unname(sym)); keep <- !is.na(sym)
cm <- rowsum(cm[keep, ], sym[keep])
note("mapped genes: ", nrow(cm))
prep <- function(m) { y <- DGEList(m); y <- calcNormFactors(y, method = "TMM"); kg <- rowMeans(cpm(y) < 1) <= 0.9
  y <- y[kg, , keep.lib.sizes = FALSE]; y <- calcNormFactors(y); cpm(y, log = TRUE, prior.count = 1) }
BAT <- c("A", "B")
dat <- lapply(setNames(BAT, BAT), function(b) { s <- d1[batch == b]; list(expr = prep(cm[, s$SampleName]), ph = s) })
rm(cm); invisible(gc())
if (DRY_RUN) for (b in BAT) { ph <- copy(dat[[b]]$ph); o <- sample(nrow(ph)); ph[, c("death", "death28") := .(death[o], death28[o])]; dat[[b]]$ph <- ph }
Y <- lapply(dat, function(d) d$ph$death)
for (b in BAT) note("batch ", b, ": ", ncol(dat[[b]]$expr), " patients, ", sum(Y[[b]]), " deaths, ", nrow(dat[[b]]$expr), " genes after filter, ",
                    uniqueN(dat[[b]]$ph$centre), " centres")

RM <- readRDS(rf); M <- load_models(); VAR <- P$main_weights
pick <- function(r, v) switch(v, multicox = list(genes = r$genes, coef = r$coef), lasso = list(genes = r$lasso_genes, coef = r$lasso_coef),
                             stepwise = list(genes = r$step_genes, coef = r$step_coef))
pub_ids <- M$mt$model_id; rows <- list()
for (id in pub_ids) rows[[paste0("PUB_", id)]] <- list(type = "PUB", matched = id)
for (r in RM$R2) { w <- pick(r, VAR); rows[[r$key]] <- list(type = "R2", matched = r$matched, genes = w$genes, coef = w$coef) }
for (v in setdiff(c("lasso", "stepwise", "multicox"), VAR)) for (r in RM$R2) { w <- pick(r, v); if (length(w$genes) >= 1)
  rows[[paste0(r$key, "_", v)]] <- list(type = paste0("R2_", v), matched = r$matched, genes = w$genes, coef = w$coef) }
for (r in RM$R2b) { w <- pick(r, VAR); rows[[r$key]] <- list(type = "R2_nodeg", matched = r$matched, genes = w$genes, coef = w$coef) }
for (r in RM$R1) rows[[r$key]] <- list(type = "R1", k = r$k, genes = r$genes, coef = r$coef)
for (r in RM$R3) { w <- pick(r, VAR); rows[[r$key]] <- list(type = "R3", genes = w$genes, coef = w$coef) }
meta <- data.table(key_ = names(rows), type = sapply(rows, `[[`, "type"),
                   matched = sapply(rows, function(r) if (is.null(r$matched)) NA_character_ else r$matched),
                   k = sapply(rows, function(r) if (is.null(r$k)) NA_integer_ else as.integer(r$k)))
setnames(meta, "key_", "key")
SM <- lapply(dat, function(d) { e <- d$expr
  t(sapply(rows, function(r) { s <- if (r$type == "PUB") score_published(M, r$matched, e) else lin_score(e, r$genes, r$coef, "center")
    if (is.null(s)) rep(NA_real_, ncol(e)) else s })) })
note("score matrices: ", nrow(SM[[1]]), " models; published signatures scorable in both batches: ",
     sum(sapply(pub_ids, function(id) all(sapply(SM, function(S) !anyNA(S[paste0("PUB_", id), ]))))), " of 36")
cov_tab <- rbindlist(lapply(pub_ids, function(id) { g <- unique(M$gn[model_id == id]$gene)
  data.table(model_id = id, n_genes = length(g), coverage_A = mean(g %in% rownames(dat$A$expr)), coverage_B = mean(g %in% rownames(dat$B$expr))) }))
fwrite(cov_tab, file.path(od, "cmaise_signature_coverage.csv"))
for (i in which(pmin(cov_tab$coverage_A, cov_tab$coverage_B) < 0.8)) note("coverage < 80%: ", cov_tab$model_id[i], " A = ", round(cov_tab$coverage_A[i], 3), ", B = ", round(cov_tab$coverage_B[i], 3))

rowDelong <- function(S, y) {
  pos <- which(y == 1); neg <- which(y == 0); n1 <- length(pos); n0 <- length(neg)
  Ra <- rowRanks(S[, c(pos, neg), drop = FALSE], ties.method = "average")
  V10 <- (Ra[, seq_len(n1), drop = FALSE] - rowRanks(S[, pos, drop = FALSE], ties.method = "average")) / n0
  V01 <- 1 - (Ra[, n1 + seq_len(n0), drop = FALSE] - rowRanks(S[, neg, drop = FALSE], ties.method = "average")) / n1
  cbind(auc = rowMeans(V10), se = sqrt(rowVars(V10) / n1 + rowVars(V01) / n0))
}
pool_mat <- function(AUC, SE, mk = matrix(TRUE, nrow(AUC), ncol(AUC))) { L <- qlogis(AUC); W <- 1 / (SE / (AUC * (1 - AUC)))^2
  bad <- !mk | !is.finite(L) | !is.finite(W) | is.na(W); W[bad] <- 0; L[bad] <- 0
  sw <- rowSums(W); out <- plogis(rowSums(W * L) / sw); out[sw == 0] <- NA; out }
k_of <- sapply(pub_ids, function(id) nrow(M$gn[model_id == id][!duplicated(gene)]))
r1_rows <- which(meta$type == "R1")
SCB <- sapply(BAT, function(b) sapply(pub_ids, function(id) !anyNA(SM[[b]][paste0("PUB_", id), ])))
MK <- matrix(TRUE, nrow(meta), length(BAT), dimnames = list(meta$key, BAT))
hasm <- !is.na(meta$matched); MK[hasm, ] <- SCB[meta$matched[hasm], , drop = FALSE]
note("published signatures not scorable in batch A / B: ", paste(pub_ids[!SCB[, "A"]], collapse = ", "), " / ", paste(pub_ids[!SCB[, "B"]], collapse = ", "),
     "; matched random rows excluded from at least one batch: ", sum(hasm & !apply(MK, 1, all)))
auc_se <- function(idx, Yl, rsel) {
  A <- lapply(names(idx), function(b) rowDelong(SM[[b]][rsel, idx[[b]], drop = FALSE], Yl[[b]][idx[[b]]]))
  AUC <- sapply(A, function(a) a[, "auc"]); SE <- sapply(A, function(a) a[, "se"]); dim(AUC) <- dim(SE) <- c(length(rsel), length(idx))
  list(AUC = AUC, SE = SE)
}
pooled_once <- function(idx, Yl, rsel = seq_len(nrow(meta))) {
  a <- auc_se(idx, Yl, rsel); setNames(pool_mat(a$AUC, a$SE, MK[rsel, names(idx), drop = FALSE]), meta$key[rsel])
}
stats_once <- function(idx) {
  a <- auc_se(idx, Y, seq_len(nrow(meta))); bs <- names(idx)
  p <- setNames(pool_mat(a$AUC, a$SE, MK[, bs, drop = FALSE]), meta$key)
  combo <- sapply(pub_ids, function(id) paste(k_of[[id]], paste(SCB[id, bs], collapse = ""), sep = "|"))
  cm <- sapply(unique(combo), function(cb) { id <- names(combo)[combo == cb][1]; rr <- r1_rows[meta$k[r1_rows] == k_of[[id]]]
    if (!length(rr) || !any(SCB[id, bs])) return(NA)
    m1 <- matrix(SCB[id, bs], length(rr), length(bs), byrow = TRUE)
    median(pool_mat(a$AUC[rr, , drop = FALSE], a$SE[rr, , drop = FALSE], m1), na.rm = TRUE) })
  list(pooled = p, r1_med = setNames(cm[combo], pub_ids))
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
boot_idx <- function(Yl, nb) { set.seed(SEED + 100); replicate(nb, strat_idx(Yl), simplify = FALSE) }

st0 <- stats_once(full_idx(Y)); est <- sapply(SUBSETS, function(s) summ(st0, s))
IDX <- boot_idx(Y, B)
btl <- mclapply(IDX, function(ix) sapply(SUBSETS, function(s) summ(stats_once(ix), s)), mc.cores = NCORES)
bt <- simplify2array(btl); saveRDS(list(st0 = st0, bt = bt), file.path(od, "bootstrap_cmaise.rds"))
q <- function(a, p) apply(a, c(1, 2), quantile, p, na.rm = TRUE)
lo <- q(bt, 0.025); hi <- q(bt, 0.975)
res <- rbindlist(lapply(colnames(est), function(s) data.table(subset = s, stat = rownames(est), estimate = est[, s], lo = lo[, s], hi = hi[, s])))
fwrite(res, file.path(od, "cmaise_estimates.csv"))
m <- res[subset == "all"]; dC <- m[stat == "Delta"]; e1 <- m[stat == "E1"]
mainbt <- readRDS(file.path(PROJ, "results/part1/bootstrap_main.rds"))
dmain0 <- summ(mainbt$st0)["Delta"]; dmain_bt <- mainbt$bt["Delta", "all", ]
dd_bt <- bt["Delta", "all", ] - dmain_bt[seq_len(dim(bt)[3])]
L <- data.table(item = c("L1", "L2", "difference_vs_main", "consistency_replication_rule"),
                estimate = c(dC$estimate, e1$estimate, dC$estimate - dmain0, NA),
                lo = c(dC$lo, e1$lo, quantile(dd_bt, 0.025, na.rm = TRUE), NA), hi = c(dC$hi, e1$hi, quantile(dd_bt, 0.975, na.rm = TRUE), NA),
                rule = c("upper 95% CI of Delta < 0.05", "lower 95% CI of E1 > 0.50", "Delta_CMAISE - Delta_main (paired bootstrap)",
                         "same sign as main Delta and |difference| < 0.05"))
L[, met := c(dC$hi < P$h1_margin, e1$lo > 0.5, NA, sign(dC$estimate) == sign(dmain0) & abs(dC$estimate - dmain0) < 0.05)]
fwrite(L, file.path(od, "cmaise_L_result.csv")); print(L)
E3 <- rbindlist(lapply(pub_ids, function(id) { a <- st0$pooled[paste0("PUB_", id)]; r <- st0$pooled[meta[type == "R2" & matched == id]$key]
  data.table(model_id = id, pooled_auc = a, random_median = median(r, na.rm = TRUE), percentile = mean(r < a, na.rm = TRUE), n_random = sum(!is.na(r))) }))
fwrite(E3, file.path(od, "cmaise_E3_percentiles.csv"))
note(sprintf("E3: %d of %d published signatures above the 95th percentile of their matched random signatures", sum(E3$percentile > 0.95, na.rm = TRUE), nrow(E3)))
pa <- lapply(BAT, function(b) { a <- rowDelong(SM[[b]], Y[[b]]); d <- data.table(key_ = meta$key, type = meta$type, batch = b, auc = a[, "auc"], se = a[, "se"]); setnames(d, "key_", "model"); d })
fwrite(rbindlist(pa)[type %in% c("PUB", "R3")], file.path(od, "cmaise_auc_by_batch_pub_R3.csv"))

rs <- which(meta$type %in% c("PUB", "R2"))
summ2 <- function(p) { pub <- median(p[paste0("PUB_", pub_ids)], na.rm = TRUE); r2 <- median(p[meta[type == "R2"]$key], na.rm = TRUE)
  c(pub_median = pub, E1 = r2, Delta = pub - r2) }
ctr_auc <- function(idx, Yl) {
  num <- 0; den <- 0
  for (b in names(idx)) { ii <- idx[[b]]; y <- Yl[[b]][ii]; ct <- dat[[b]]$ph$centre[ii]; S <- SM[[b]][rs, ii, drop = FALSE]
    mk <- MK[rs, b]
    for (cc in unique(ct)) { j <- which(ct == cc); n1 <- sum(y[j] == 1); n0 <- sum(y[j] == 0); if (!n1 || !n0) next
      R <- rowRanks(S[, j, drop = FALSE], ties.method = "average")
      u <- rowSums(R[, y[j] == 1, drop = FALSE]) - n1 * (n1 + 1) / 2
      num <- num + ifelse(mk, u, 0); den <- den + n1 * n0 * mk } }
  out <- num / den; out[den == 0] <- NA; setNames(out, meta$key[rs])
}
run_sens <- function(label, Yl, sel = NULL, fun = "pooled") {
  bs <- if (is.null(sel)) names(Yl) else names(sel)
  base <- lapply(setNames(bs, bs), function(b) if (is.null(sel)) seq_along(Yl[[b]]) else sel[[b]])
  sub_idx <- function(ix) lapply(setNames(bs, bs), function(b) base[[b]][ix[[b]]])
  Ys <- lapply(setNames(bs, bs), function(b) Yl[[b]][base[[b]]])
  f <- function(ix) { full <- sub_idx(ix); p <- if (fun == "centre") ctr_auc(full, Yl) else pooled_once(full, Yl, rs); summ2(p) }
  e0 <- f(full_idx(Ys)); ib <- boot_idx(Ys, B)
  b2 <- simplify2array(mclapply(ib, f, mc.cores = NCORES))
  data.table(analysis = label, stat = names(e0), estimate = e0, lo = apply(b2, 1, quantile, 0.025, na.rm = TRUE), hi = apply(b2, 1, quantile, 0.975, na.rm = TRUE),
             n = sum(lengths(base)), deaths = sum(sapply(bs, function(b) sum(Yl[[b]][base[[b]]]))))
}
Y28 <- lapply(dat, function(d) d$ph$death28)
shock <- lapply(dat, function(d) which(d$ph$group %in% c("shock", "refractory")))
sens <- rbind(run_sens("centre_stratified", Y, fun = "centre"),
              run_sens("shock_and_refractory_only", Y, sel = shock),
              run_sens("death_within_28_days", Y28),
              run_sens("batch_A_only", Y["A"]),
              run_sens("batch_B_only", Y["B"]))
fwrite(sens, file.path(od, "cmaise_sensitivity.csv")); print(sens)

sofa_b <- rbindlist(lapply(BAT, function(b) { a <- delong(dat[[b]]$ph$SOFA, Y[[b]]); data.table(batch = b, auc = a[1], se = a[2]) }))
sofa_pool <- pool_fixed(sofa_b$auc, sofa_b$se)
ally <- unlist(Y, use.names = FALSE); allb <- factor(rep(BAT, sapply(Y, length))); alls <- unlist(lapply(dat, function(d) d$ph$SOFA), use.names = FALSE)
auc_lp <- function(f) delong(predict(f), ally)[1]
base_auc <- auc_lp(glm(ally ~ allb + alls, family = binomial))
inc <- sapply(meta$key[rs], function(k) { sc <- unlist(lapply(BAT, function(b) SM[[b]][k, ]), use.names = FALSE); if (anyNA(sc) || !all(MK[k, ])) return(NA)
  auc_lp(glm(ally ~ allb + alls + sc, family = binomial)) - base_auc })
sofa_out <- data.table(item = c("SOFA_pooled_auc", "SOFA_pooled_lo", "SOFA_pooled_hi", "base_auc_batch_SOFA_insample",
                                "median_increment_published", "median_increment_random", "share_published_increment_gt_0.01", "share_random_increment_gt_0.01"),
                       value = c(sofa_pool, base_auc, median(inc[grepl("^PUB_", names(inc))], na.rm = TRUE), median(inc[!grepl("^PUB_", names(inc))], na.rm = TRUE),
                                 mean(inc[grepl("^PUB_", names(inc))] > 0.01, na.rm = TRUE), mean(inc[!grepl("^PUB_", names(inc))] > 0.01, na.rm = TRUE)))
fwrite(sofa_out, file.path(od, "cmaise_sofa.csv")); fwrite(sofa_b, file.path(od, "cmaise_sofa_by_batch.csv")); print(sofa_out)

stopifnot(md5_of(LM22_FILE) == LM22_MD5)
lm22 <- as.matrix(fread(LM22_FILE), rownames = 1); rownames(lm22) <- official(rownames(lm22)); lm22 <- lm22[!duplicated(rownames(lm22)), ]
GROUPS <- list(neutrophil = "Neutrophils", eosinophil = "Eosinophils",
  mono_mac_dc = c("Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2", "Dendritic cells resting", "Dendritic cells activated"),
  t_cell = c("T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting", "T cells CD4 memory activated",
             "T cells follicular helper", "T cells regulatory (Tregs)", "T cells gamma delta"),
  nk = c("NK cells resting", "NK cells activated"), b_plasma = c("B cells naive", "B cells memory", "Plasma cells"),
  mast = c("Mast cells resting", "Mast cells activated"))
COVARS <- c("neutrophil", "mono_mac_dc", "t_cell", "nk", "b_plasma")
nnls1 <- function(A, b) { f <- function(x) { r <- A %*% x - b; sum(r^2) }; g <- function(x) as.numeric(2 * crossprod(A, A %*% x - b))
  optim(rep(1 / ncol(A), ncol(A)), f, g, method = "L-BFGS-B", lower = 0, control = list(maxit = 1000, factr = 1e5))$par }
deconv <- function(expr) { gg <- intersect(rownames(lm22), rownames(expr)); A <- lm22[gg, ]; Yx <- 2^expr[gg, , drop = FALSE]
  A <- A / max(A); Yx <- Yx / max(Yx)
  Pm <- t(apply(Yx, 2, function(b) { x <- nnls1(A, b); if (sum(x) > 0) x / sum(x) else x })); colnames(Pm) <- colnames(lm22)
  list(groups = sapply(GROUPS, function(cl) rowSums(Pm[, cl, drop = FALSE])), n_genes = length(gg)) }
mk_sets <- fread(file.path(PROJ, "results/composition/marker_sets_sensitivity.csv"))
tnlr <- function(G) log((G[, "neutrophil"] + 1e-3) / (G[, "t_cell"] + G[, "nk"] + G[, "b_plasma"] + 1e-3))
COMP <- lapply(dat, function(d) { dc <- deconv(d$expr); G <- dc$groups
  mk <- sapply(split(mk_sets$gene, mk_sets$group), function(gs) { gs <- intersect(gs, rownames(d$expr)); colMeans(d$expr[gs, , drop = FALSE]) })
  list(G = G, mk = mk, tnlr = tnlr(G), n_genes = dc$n_genes) })
for (b in BAT) note("batch ", b, " LM22 genes used: ", COMP[[b]]$n_genes, "; median neutrophil ", round(median(COMP[[b]]$G[, "neutrophil"]), 3))
fwrite(rbindlist(lapply(BAT, function(b) data.table(batch = b, sample = colnames(dat[[b]]$expr), COMP[[b]]$G, tnlr = COMP[[b]]$tnlr))), file.path(od, "cmaise_composition.csv"))
pubk <- paste0("PUB_", pub_ids)
r2_med <- function(idx, cv = "lm22") median(sapply(pubk, function(k) {
  r2 <- sapply(BAT, function(b) { ii <- idx[[b]]; s <- SM[[b]][k, ii]; if (anyNA(s)) return(NA_real_)
    X <- if (cv == "lm22") COMP[[b]]$G[ii, COVARS] else COMP[[b]]$mk[ii, c("neutrophil", "mono_mac_dc", "t_cell", "nk", "b_plasma")]
    summary(lm(s ~ X))$adj.r.squared })
  w <- sapply(idx, length) * !is.na(r2); if (sum(w) == 0) NA_real_ else sum(r2 * w, na.rm = TRUE) / sum(w) }), na.rm = TRUE)
nid <- lapply(Y, seq_along); set.seed(SEED + 101); ridx <- replicate(B, lapply(nid, function(i) sample(i, replace = TRUE)), simplify = FALSE)
r2_0 <- r2_med(nid); r2_bt <- unlist(mclapply(ridx, r2_med, mc.cores = NCORES))
mk_0 <- r2_med(nid, "marker")
resid_mat <- function(S, X) { X <- cbind(1, X); S - t(X %*% solve(crossprod(X), crossprod(X, t(S)))) }
SMadj <- lapply(setNames(BAT, BAT), function(b) { S <- SM[[b]][rs, , drop = FALSE]; ok <- rowSums(is.na(S)) == 0; out <- S; out[ok, ] <- resid_mat(S[ok, , drop = FALSE], COMP[[b]]$G[, COVARS]); out })
padj <- { A <- lapply(BAT, function(b) rowDelong(SMadj[[b]], Y[[b]])); AUC <- sapply(A, function(a) a[, "auc"]); SE <- sapply(A, function(a) a[, "se"]); setNames(pool_mat(AUC, SE, MK[rs, BAT, drop = FALSE]), meta$key[rs]) }
praw <- st0$pooled[meta$key[rs]]
tn <- rbindlist(lapply(BAT, function(b) { a <- delong(COMP[[b]]$tnlr, Y[[b]]); data.table(batch = b, auc = a[1], se = a[2]) })); tn_pool <- pool_fixed(tn$auc, tn$se)
p2 <- data.table(item = c("median_adj_r2_lm22", "median_adj_r2_lm22_basic_lo", "median_adj_r2_lm22_basic_hi", "median_adj_r2_marker",
                          "published_median_raw", "published_median_adj", "random_median_raw", "random_median_adj",
                          "tnlr_pooled_auc", "tnlr_pooled_lo", "tnlr_pooled_hi"),
                 value = c(r2_0, 2 * r2_0 - quantile(r2_bt, 0.975), 2 * r2_0 - quantile(r2_bt, 0.025), mk_0,
                           median(praw[grepl("^PUB_", names(praw))], na.rm = TRUE), median(padj[grepl("^PUB_", names(padj))], na.rm = TRUE),
                           median(praw[!grepl("^PUB_", names(praw))], na.rm = TRUE), median(padj[!grepl("^PUB_", names(padj))], na.rm = TRUE), tn_pool))
fwrite(p2, file.path(od, "cmaise_part2_descriptive.csv")); print(p2)

writeLines(log_lines, file.path(od, "cmaise_run_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(od, "sessionInfo_08.txt"))
message("done")
