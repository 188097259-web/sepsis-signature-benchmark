source("scripts/_common.R")
set.seed(SEED)
od <- outdir("composition")
stopifnot(md5_of(LM22_FILE) == LM22_MD5)
lm22 <- as.matrix(fread(LM22_FILE), rownames = 1); rownames(lm22) <- official(rownames(lm22))
lm22 <- lm22[!duplicated(rownames(lm22)), ]
GROUPS <- list(
  neutrophil = "Neutrophils",
  eosinophil = "Eosinophils",
  mono_mac_dc = c("Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2", "Dendritic cells resting", "Dendritic cells activated"),
  t_cell = c("T cells CD8", "T cells CD4 naive", "T cells CD4 memory resting", "T cells CD4 memory activated",
             "T cells follicular helper", "T cells regulatory (Tregs)", "T cells gamma delta"),
  nk = c("NK cells resting", "NK cells activated"),
  b_plasma = c("B cells naive", "B cells memory", "Plasma cells"),
  mast = c("Mast cells resting", "Mast cells activated"))
stopifnot(setequal(unlist(GROUPS), colnames(lm22)))
COVARS <- c("neutrophil", "mono_mac_dc", "t_cell", "nk", "b_plasma")

nnls1 <- function(A, b) {
  f <- function(x) { r <- A %*% x - b; sum(r^2) }; g <- function(x) as.numeric(2 * crossprod(A, A %*% x - b))
  x0 <- rep(1 / ncol(A), ncol(A))
  optim(x0, f, g, method = "L-BFGS-B", lower = 0, control = list(maxit = 1000, factr = 1e5))$par
}
deconv <- function(expr) {
  gg <- intersect(rownames(lm22), rownames(expr)); A <- lm22[gg, ]; Y <- 2^expr[gg, , drop = FALSE]
  sA <- max(A); A <- A / sA; Y <- Y / max(Y)
  P <- t(apply(Y, 2, function(b) { x <- nnls1(A, b); if (sum(x) > 0) x / sum(x) else x }))
  colnames(P) <- colnames(lm22)
  G <- sapply(GROUPS, function(cl) rowSums(P[, cl, drop = FALSE]))
  list(cells = P, groups = G, n_genes = length(gg))
}
marker_sets <- lapply(GROUPS[c("neutrophil", "eosinophil", "mono_mac_dc", "t_cell", "nk", "b_plasma")], function(cl) {
  own <- apply(lm22[, cl, drop = FALSE], 1, max); oth <- apply(lm22[, setdiff(colnames(lm22), cl), drop = FALSE], 1, max)
  names(sort(own / pmax(oth, 1), decreasing = TRUE))[1:10] })
fwrite(data.table(group = rep(names(marker_sets), lengths(marker_sets)), gene = unlist(marker_sets)), file.path(od, "marker_sets_sensitivity.csv"))
marker_scores <- function(expr) sapply(marker_sets, function(gs) { gs <- intersect(gs, rownames(expr)); colMeans(expr[gs, , drop = FALSE]) })
tnlr <- function(G) log((G[, "neutrophil"] + 1e-3) / (G[, "t_cell"] + G[, "nk"] + G[, "b_plasma"] + 1e-3))

COH <- c(MAIN, REPL, SENS, TRAIN)
info <- list()
for (g in COH) {
  x <- load_cohort(g, sepsis_only = FALSE); e <- x$expr
  if (DRY_RUN) rownames(e) <- sample(rownames(e))
  d <- deconv(e); ms <- marker_scores(e)
  out <- data.table(cohort = g, sample = colnames(e), group = x$ph$group, d$groups, tnlr = tnlr(d$groups))
  setnames(ms <- as.data.table(ms), paste0("mk_", colnames(ms))); out <- cbind(out, ms)
  fwrite(out, file.path(od, paste0("composition_", g, ".csv")))
  fwrite(data.table(sample = colnames(e), d$cells), file.path(od, paste0("lm22_cells_", g, ".csv")))
  info[[g]] <- data.table(cohort = g, n_samples = ncol(e), lm22_genes_used = d$n_genes,
                          median_neutrophil = median(d$groups[, "neutrophil"]), median_lymph = median(rowSums(d$groups[, c("t_cell", "nk", "b_plasma")])))
  message(g, ": LM22 genes ", d$n_genes, ", median neutrophil ", round(median(d$groups[, "neutrophil"]), 3))
}
fwrite(rbindlist(info), file.path(od, "deconvolution_summary.csv"))

pd <- fread(file.path(DATA, "cohorts/pdata_GSE54514.csv"))
mcol <- grep("neutrophil", names(pd), value = TRUE, ignore.case = TRUE)[1]
idcol <- intersect(c("geo_accession", "sample", "gsm"), names(pd))[1]
cm <- fread(file.path(od, "composition_GSE54514.csv"))
m <- merge(cm, data.table(sample = pd[[idcol]], measured = suppressWarnings(as.numeric(pd[[mcol]]))), by = "sample")
m <- m[!is.na(measured)]
gate <- data.table(n = nrow(m), rho_nnls = cor(m$neutrophil, m$measured, method = "spearman"),
                   rho_marker = cor(m$mk_neutrophil, m$measured, method = "spearman"))
gate[, pass := rho_nnls >= 0.5]
fwrite(gate, file.path(od, "quality_gate_GSE54514.csv"))
message(sprintf("quality gate: n = %d, Spearman (NNLS) = %.3f, (marker) = %.3f, pass = %s", gate$n, gate$rho_nnls, gate$rho_marker, gate$pass))
writeLines(capture.output(sessionInfo()), file.path(od, "sessionInfo.txt"))
