source("scripts/_common.R")
suppressPackageStartupMessages({ library(AnnotationDbi); library(metafor); library(SepstratifieR); library(ConsensusTranscriptomicSubtype) })
set.seed(SEED)
od <- outdir("batch"); comp_dir <- file.path(OUT, "composition")
DSETS <- c("GSE185263", "GSE236713", "GSE95233", "GSE54514", "GSE65682")
CTS_EXCLUDE <- "GSE65682"
M <- load_models(); gn <- M$gn
bench_gn <- fread(file.path(DATA, "signatures/genes.csv"), colClasses = "character"); bench_gn[, gene := official(gene_symbol)]
genes_of <- function(id) bench_gn[model_id == id][!duplicated(gene)]

gmean_rows <- function(m) exp(colMeans(log(m)))
SUB_SHIFT <- 0
subspace_raw <- function(expr) {
  my <- genes_of("Moore2025_myeloid_dysregulation"); ly <- genes_of("Moore2025_lymphoid_dysregulation"); m <- expr + SUB_SHIFT
  up <- intersect(my[direction == "risk"]$gene, rownames(m)); dn <- intersect(my[direction == "protective"]$gene, rownames(m))
  lg <- intersect(ly$gene, rownames(m))
  cbind(myeloid = gmean_rows(m[up, , drop = FALSE]) - gmean_rows(m[dn, , drop = FALSE]), lymphoid = -gmean_rows(m[lg, , drop = FALSE]))
}
quadrant <- function(zm, zl, thr) ifelse(zm >= thr & zl >= thr, "system", ifelse(zm >= thr, "myeloid", ifelse(zl >= thr, "lymphoid", "balanced")))
subspace_healthy <- function(raw, is_ctrl) {
  z <- sweep(sweep(raw, 2, colMeans(raw[is_ctrl, , drop = FALSE])), 2, apply(raw[is_ctrl, , drop = FALSE], 2, sd), "/")
  quadrant(z[, "myeloid"], z[, "lymphoid"], 1.65)
}
subspace_median <- function(raw, is_pt) {
  med <- apply(raw[is_pt, , drop = FALSE], 2, median)
  quadrant(raw[, "myeloid"] - med["myeloid"], raw[, "lymphoid"] - med["lymphoid"], 0)
}
sweeney_raw <- function(expr, id) {
  g <- genes_of(id); lin <- 2^expr; up <- intersect(g[direction == "risk"]$gene, rownames(lin)); dn <- intersect(g[direction == "protective"]$gene, rownames(lin))
  s <- gmean_rows(lin[up, , drop = FALSE]) - gmean_rows(lin[dn, , drop = FALSE]); if (id == "Sweeney2018NC_Stanford12") s <- gmean_rows(lin[up, , drop = FALSE]) - gmean_rows(lin[dn, , drop = FALSE]) * length(up) / length(dn); s
}
mars1 <- function(expr) if (all(c("BPGM", "TAP2") %in% rownames(expr))) expr["BPGM", ] - expr["TAP2", ] else NULL
srs <- function(expr, set = "davenport", k = NULL) {
  need <- colnames(if (set == "davenport") SepstratifieR::reference_set_davenport else SepstratifieR::reference_set_extended)
  sym <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = need, columns = "SYMBOL", keytype = "ENSEMBL")); sym$SYMBOL <- official(sym$SYMBOL)
  hit <- sapply(need, function(e) { c0 <- intersect(sym$SYMBOL[sym$ENSEMBL == e], rownames(expr)); if (length(c0)) c0[1] else NA_character_ })
  if (anyNA(hit)) return(NULL)
  dat <- as.data.frame(t(expr[hit, , drop = FALSE])); colnames(dat) <- need
  if (is.null(k)) k <- max(5, round(0.25 * nrow(dat)))
  p <- suppressMessages(stratifyPatients(dat, gene_set = set, k = k, verbose = FALSE))
  data.table(sample = colnames(expr), SRS = as.character(p@SRS), SRSq = p@SRSq, mnn_outlier = p@mNN_outlier)
}
cts <- function(expr, seed) {
  g <- genes_of("Scicluna2025_CTS18_RF")$gene; pres <- intersect(g, rownames(expr)); if (length(pres) < 0.8 * length(g)) return(NULL)
  en <- mapIds(org.Hs.eg.db, keys = pres, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
  m <- expr[pres, , drop = FALSE]; rownames(m) <- unname(en); set.seed(seed); pdf(NULL); on.exit(dev.off())
  r <- suppressMessages(suppressWarnings(run_subtype_classifier(new_expr_data = m)))
  as.character(r$predictions$CTS[match(colnames(m), r$predictions$Sample)])
}
zmodel <- function(expr, id) score_published(M, id, expr)

run_batch <- function(expr, is_ctrl, is_pt, tag, cts_ok, cts_seeds = 1) {
  out <- list(); add <- function(tool, value, type) out[[length(out) + 1]] <<- data.table(sample = colnames(expr), tool = tool, value = as.character(value), type = type)
  raw <- subspace_raw(expr)
  add("SUBSPACE_myeloid_raw", raw[, "myeloid"], "continuous"); add("SUBSPACE_lymphoid_raw", raw[, "lymphoid"], "continuous")
  if (sum(is_ctrl) >= 3) add("SUBSPACE_quadrant_healthyref", subspace_healthy(raw, is_ctrl), "label")
  if (sum(is_pt) >= 3) add("SUBSPACE_quadrant_median", subspace_median(raw, is_pt), "label")
  for (id in c("Sweeney2018NC_Duke18", "Sweeney2018NC_Stanford12")) { s <- sweeney_raw(expr, id)
    add(paste0(id, "_raw"), s, "continuous"); add(paste0(id, "_batchz"), as.numeric(scale(s)), "continuous") }
  s <- mars1(expr); if (!is.null(s)) add("Mars1_ratio", s, "continuous")
  for (id in Z_MODELS) { s <- zmodel(expr, id); if (!is.null(s)) add(paste0(id, "_batchz"), s, "continuous") }
  for (set in c("davenport", "extended")) { r <- tryCatch(srs(expr, set), error = function(e) { message("SRS ", tag, ": ", conditionMessage(e)); NULL })
    if (!is.null(r)) { add(paste0("SRS_", set), r$SRS, "label"); add(paste0("SRSq_", set), r$SRSq, "continuous"); add(paste0("SRS_", set, "_mnn_outlier"), r$mnn_outlier, "flag") } }
  r <- tryCatch(srs(expr, "davenport", k = 20), error = function(e) NULL); if (!is.null(r)) add("SRS_davenport_k20", r$SRS, "label")
  if (cts_ok) for (i in seq_len(cts_seeds)) { r <- tryCatch(cts(expr, SEED + i - 1), error = function(e) { message("CTS ", tag, ": ", conditionMessage(e)); NULL })
    if (!is.null(r)) add(if (i == 1) "CTS" else paste0("CTS_seed", i), r, "label") }
  res <- rbindlist(out); res[, batch := tag]; res
}

long_file <- file.path(od, "batch_outputs_long.csv")
REUSE <- identical(Sys.getenv("REUSE_BATCH"), "1") && file.exists(long_file)
all_res <- list()
if (!REUSE) for (g in DSETS) {
  x <- load_cohort(g, sepsis_only = FALSE); e <- x$expr; ph <- x$ph
  if (DRY_RUN) rownames(e) <- sample(rownames(e))
  SUB_SHIFT <- { mn <- min(e, na.rm = TRUE); if (mn < 1) 1 - mn else 0 }
  message(g, ": SUBSPACE shift constant ", round(SUB_SHIFT, 3))
  idx_pt <- which(ph$group == "sepsis" & (g != "GSE185263" | ph$setting == "ICU")); idx_hc <- which(ph$group == "control")
  idx_ed <- if (g == "GSE185263") which(ph$group == "sepsis" & ph$setting == "Emergency Room") else integer(0)
  cts_ok <- !(g %in% CTS_EXCLUDE)
  B <- function(ix, tag, seeds = 1) { r <- run_batch(e[, ix, drop = FALSE], ph$group[ix] == "control", ph$group[ix] == "sepsis", paste(g, tag), cts_ok, seeds); r[, `:=`(cohort = g, scenario = tag)]; r }
  all_res[[paste(g, "C0")]] <- B(idx_pt, "C0", seeds = 10)
  all_res[[paste(g, "C1")]] <- B(c(idx_pt, idx_hc), "C1")
  all_res[[paste(g, "C4")]] <- B(idx_hc, "C4_healthy_alone")
  if (length(idx_ed)) all_res[[paste(g, "C2")]] <- B(c(idx_pt, idx_ed), "C2")
  if (!(g %in% NNLS_INVALID)) {
    cp <- fread(file.path(comp_dir, paste0("composition_", g, ".csv")))[match(ph$sample[idx_pt], sample)]
    ter <- cut(rank(cp$tnlr, ties.method = "first"), 3, labels = c("low", "mid", "high"))
    all_res[[paste(g, "C3hi")]] <- B(idx_pt[ter %in% c("mid", "high")], "C3_with_high")
    all_res[[paste(g, "C3lo")]] <- B(idx_pt[ter %in% c("mid", "low")], "C3_with_low")
    fwrite(data.table(cohort = g, sample = ph$sample[idx_pt], tnlr_tertile = as.character(ter)), file.path(od, paste0("tertiles_", g, ".csv")))
  }
  raw <- subspace_raw(e[, c(idx_pt, idx_hc), drop = FALSE]); is_c <- c(rep(FALSE, length(idx_pt)), rep(TRUE, length(idx_hc)))
  ref_q <- subspace_healthy(raw, is_c)[!is_c]
  ch <- replicate(200, { keep <- c(rep(TRUE, length(idx_pt)), seq_along(idx_hc) %in% sample(seq_along(idx_hc), min(10, length(idx_hc))))
    q <- subspace_healthy(raw[keep, , drop = FALSE], is_c[keep])[seq_along(idx_pt)]; mean(q != ref_q) })
  fwrite(data.table(cohort = g, n_controls = length(idx_hc), draw = seq_along(ch), change = ch), file.path(od, paste0("C5_subspace_controls_", g, ".csv")))
  message(g, " done")
}
if (REUSE) res <- fread(long_file, colClasses = list(character = "value")) else { res <- rbindlist(all_res, fill = TRUE); fwrite(res, long_file) }

kappa_gen <- function(a, b) { lv <- union(a, b); t <- table(factor(a, lv), factor(b, lv)); n <- sum(t); po <- sum(diag(t)) / n; pe <- sum(rowSums(t) * colSums(t)) / n^2; if (pe == 1) NA else (po - pe) / (1 - pe) }
compare <- function(r, g, s1, s2, tl, set = NULL) {
  a <- r[cohort == g & scenario == s1 & tool == tl]; b <- r[cohort == g & scenario == s2 & tool == tl]; tool <- tl
  common <- intersect(a$sample, b$sample); if (!is.null(set)) common <- intersect(common, set); if (!length(common)) return(NULL)
  va <- a[match(common, sample)]$value; vb <- b[match(common, sample)]$value; type <- a$type[1]
  if (type == "label") data.table(cohort = g, contrast = paste(s1, "vs", s2), tool = tool, n = length(common), changed = sum(va != vb), prop_changed = mean(va != vb), kappa = kappa_gen(va, vb))
  else if (type == "continuous") { va <- as.numeric(va); vb <- as.numeric(vb); sdv <- sd(va)
    data.table(cohort = g, contrast = paste(s1, "vs", s2), tool = tool, n = length(common), median_abs_shift_sd = median(abs(vb - va)) / sdv,
               spearman = cor(va, vb, method = "spearman"), median_abs_rank_pct_shift = median(abs(rank(vb) / length(vb) - rank(va) / length(va)))) }
}
tools <- unique(res[type %in% c("label", "continuous") & !grepl("_seed", tool)]$tool)
cmp <- list()
for (g in DSETS) {
  r <- res[cohort == g]; has_c3 <- !(g %in% NNLS_INVALID)
  mid <- if (has_c3) fread(file.path(od, paste0("tertiles_", g, ".csv")))[tnlr_tertile == "mid"]$sample else character(0)
  for (tl in tools) {
    cmp[[length(cmp) + 1]] <- compare(r, g, "C0", "C1", tl)
    cmp[[length(cmp) + 1]] <- compare(r, g, "C0", "C2", tl)
    if (has_c3) cmp[[length(cmp) + 1]] <- compare(r, g, "C3_with_high", "C3_with_low", tl, mid)
    cmp[[length(cmp) + 1]] <- compare(r, g, "C4_healthy_alone", "C1", tl)
  }
  for (i in 2:10) cmp[[length(cmp) + 1]] <- { a <- r[scenario == "C0" & tool == "CTS"]; b <- r[scenario == "C0" & tool == paste0("CTS_seed", i)]
    if (nrow(a) && nrow(b)) data.table(cohort = g, contrast = paste0("seed1 vs seed", i), tool = "CTS", n = nrow(a), changed = sum(a$value != b$value[match(a$sample, b$sample)]), prop_changed = mean(a$value != b$value[match(a$sample, b$sample)])) }
}
cmp <- rbindlist(Filter(Negate(is.null), cmp), fill = TRUE); fwrite(cmp, file.path(od, "batch_comparisons.csv"))

hc_ids <- unique(res[scenario == "C4_healthy_alone", .(cohort, sample)])
hc_dist <- res[scenario %in% c("C4_healthy_alone", "C1") & type == "label"][hc_ids, on = .(cohort, sample), nomatch = 0][, .(n = .N), by = .(cohort, scenario, tool, value)]
fwrite(hc_dist, file.path(od, "label_distribution.csv"))

h4d <- cmp[tool == "SRS_davenport" & contrast == "C0 vs C1"]
H4_THRESHOLD <- 0.05
if (nrow(h4d) >= 2) {
  f <- rma(measure = "PLO", xi = h4d$changed, ni = h4d$n, method = "REML", test = "knha")
  p1 <- pt((f$b[1] - qlogis(H4_THRESHOLD)) / f$se, df = f$k - f$p, lower.tail = FALSE)
  h4 <- data.table(k = nrow(h4d), n = sum(h4d$n), changed = sum(h4d$changed), pooled = plogis(f$b[1]), lo = plogis(f$ci.lb), hi = plogis(f$ci.ub),
                   I2 = f$I2, threshold = H4_THRESHOLD, p_one_sided = p1, supported = plogis(f$ci.lb) > H4_THRESHOLD)
  fwrite(h4, file.path(od, "H4_result.csv")); print(h4)
}
writeLines(capture.output(sessionInfo()), file.path(od, "sessionInfo.txt"))
message("done")
