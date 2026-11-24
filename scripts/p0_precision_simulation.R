suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(limma); library(org.Hs.eg.db); library(metafor) })
set.seed(20261116)
PROJ <- Sys.getenv("SEPSIS_BENCH_ROOT", ".")
DATA <- file.path(PROJ, "data")
od <- file.path(PROJ, "results/power"); dir.create(od, recursive = TRUE, showWarnings = FALSE)
MAIN <- c("GSE185263", "GSE63042", "GSE236713")
REF_PMIDS <- c("35027333","41028542","28864056","36322631","41028543","29449546","26917434","29537985","35847809","33181462")
Z_MODELS <- c("Jiang2022_IRRG9", "Chen2023_IRG7", "Feng2024_ROS37", "Li2026_GMGscore6")
B <- 2000

official <- function(sym) { o <- suppressWarnings(alias2SymbolTable(sym, species = "Hs")); ifelse(is.na(o), sym, o) }
mt <- fread(file.path(DATA, "signatures/model_table.csv"), colClasses = "character")
mt <- mt[!(pmid %in% REF_PMIDS) & tier == "A_coef_complete"]
gn <- fread(file.path(DATA, "signatures/genes.csv"), colClasses = "character")
gn[, coef := suppressWarnings(as.numeric(coef))]; gn[, gene := official(gene_symbol)]
higher <- sapply(mt$model_id, function(id) {
  r <- fromJSON(file.path(DATA, "signatures/extraction_records", paste0(mt[model_id == id]$pmid[1], ".json")), simplifyVector = FALSE)
  m <- Filter(function(z) identical(z$model_id, id), r$models)[[1]]
  h <- if (is.null(m$higher_score_means)) "higher risk" else m$higher_score_means
  if (grepl("lower", h)) -1 else 1
})
scale_mat <- function(m, how) { mu <- rowMeans(m); if (how == "center") m - mu else (m - mu) / pmax(apply(m, 1, sd), 1e-8) }
D <- lapply(setNames(MAIN, MAIN), function(g) {
  x <- readRDS(file.path(DATA, "cohorts", paste0(g, ".rds"))); s <- x$pheno$group == "sepsis"
  e <- x$expr[, s, drop = FALSE]; list(center = scale_mat(e, "center"), z = scale_mat(e, "z"), y = x$pheno$death[s])
})
main_scale <- function(id) if (id %in% Z_MODELS) "z" else "center"
score_of <- function(id, d) {
  g <- gn[model_id == id][!duplicated(gene)]; present <- g$gene %in% rownames(d$center)
  if (mean(present) < 0.8) return(NULL); g <- g[present]
  higher[[id]] * as.numeric(crossprod(g$coef, d[[main_scale(id)]][g$gene, , drop = FALSE]))
}
eligible <- function(id, g) mt[model_id == id][[g]] == "eligible"

delong <- function(s, y) {
  pos <- s[y == 1]; neg <- s[y == 0]; n1 <- length(pos); n0 <- length(neg)
  ra <- rank(c(pos, neg)); v10 <- (ra[seq_len(n1)] - rank(pos)) / n0; v01 <- 1 - (ra[n1 + seq_len(n0)] - rank(neg)) / n1
  c(mean(v10), sqrt(var(v10) / n1 + var(v01) / n0))
}
pool_fixed <- function(a, se) { l <- qlogis(a); w <- 1 / (se / (a * (1 - a)))^2; plogis(sum(w * l) / sum(w)) }

S <- list()
for (id in mt$model_id) for (g in MAIN) if (eligible(id, g)) { s <- score_of(id, D[[g]]); if (!is.null(s)) S[[id]][[g]] <- s }
ids <- names(S); message("models with scores: ", length(ids))

pooled_once <- function(idx) {
  sapply(ids, function(id) { r <- sapply(names(S[[id]]), function(g) delong(S[[id]][[g]][idx[[g]]], D[[g]]$y[idx[[g]]]))
    pool_fixed(r[1, ], r[2, ]) })
}

full <- lapply(setNames(MAIN, MAIN), function(g) seq_along(D[[g]]$y))
p0 <- pooled_once(full)
known <- fread(file.path(DATA, "prior_work/pooled_main.csv"))
cmp <- merge(data.table(model_id = names(p0), auc_here = as.numeric(p0)), known[, .(model_id, auc_known = auc)], by = "model_id")
fwrite(cmp, file.path(od, "check_against_benchmark.csv"))
message(sprintf("check: n=%d, max |diff| = %.2e, median here = %.4f", nrow(cmp), max(abs(cmp$auc_here - cmp$auc_known)), median(p0)))

strat_idx <- function() lapply(setNames(MAIN, MAIN), function(g) { y <- D[[g]]$y
  c(sample(which(y == 1), replace = TRUE), sample(which(y == 0), replace = TRUE)) })
boot_med <- numeric(B); boot_single <- matrix(NA, B, length(ids), dimnames = list(NULL, ids))
for (b in seq_len(B)) { pb <- pooled_once(strat_idx()); boot_med[b] <- median(pb); boot_single[b, ] <- pb }
SE_pub <- sd(boot_med); SE_single <- median(apply(boot_single, 2, sd))
message(sprintf("SE of median pooled AUROC (published) = %.4f; median SE of a single model = %.4f", SE_pub, SE_single))
split_rho <- replicate(200, { h <- sample(ids, length(ids) / 2)
  cor(apply(boot_single[, h, drop = FALSE], 1, median), apply(boot_single[, setdiff(ids, h), drop = FALSE], 1, median)) })
fwrite(data.table(split = seq_along(split_rho), rho = split_rho), file.path(od, "split_half_rho.csv"))
message(sprintf("split-half rho: median %.3f (IQR %.3f-%.3f)", median(split_rho), quantile(split_rho, .25), quantile(split_rho, .75)))

h1 <- CJ(rho = sort(unique(round(c(0, 0.5, 0.8, median(split_rho), quantile(split_rho, 0.25)), 3))), delta_true = c(-0.02, 0, 0.01, 0.02), margin = c(0.03, 0.05))
h1[, se_diff := sqrt(2 * SE_pub^2 * (1 - rho))]
h1[, power := pnorm((margin - delta_true) / se_diff - qnorm(0.975))]
fwrite(h1, file.path(od, "H1_power.csv"))
h2 <- data.table(E1_true = c(0.51, 0.52, 0.53, 0.55)); h2[, power := pnorm((E1_true - 0.5) / SE_pub - qnorm(0.975))]
fwrite(h2, file.path(od, "H2_power.csv"))

hm_se <- function(A, n1, n0) { q1 <- A / (2 - A); q2 <- 2 * A^2 / (1 + A); sqrt((A * (1 - A) + (n1 - 1) * (q1 - A^2) + (n0 - 1) * (q2 - A^2)) / (n1 * n0)) }
rep_se_single <- hm_se(0.57, 60, 101)
rep <- data.table(quantity = c("single model SE (Hanley-McNeil, AUROC 0.57)", "median of models SE (scaled by main-cohort ratio)"),
                  se = c(rep_se_single, rep_se_single * SE_pub / SE_single))
rep[, ci_halfwidth := 1.96 * se]; fwrite(rep, file.path(od, "replication_precision.csv"))

nd <- sapply(MAIN, function(g) c(n = length(D[[g]]$y), d = sum(D[[g]]$y)))
fauc <- function(s, y) { r <- rank(s); n1 <- as.numeric(sum(y)); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * (length(y) - n1)) }
calib <- function(rho, rf, target = 0.568, N = 4e5) {
  C <- rnorm(N); U <- rnorm(N); Sx <- rho * C + sqrt(1 - rho^2) * U; Uu <- runif(N)
  a_adj <- 0.5 + rf * (target - 0.5)
  f <- function(p) { lp <- -1.6 + p[1] * C + p[2] * U; y <- as.integer(Uu < plogis(lp))
    c(fauc(Sx, y) - target, fauc(U, y) - a_adj) }
  o <- optim(c(0.2, 0.1), function(p) sum(f(p)^2) * 1e4, method = "Nelder-Mead", control = list(maxit = 300))
  o$par
}
one_rep <- function(par, rho, nb = 300) {
  dat <- lapply(MAIN, function(g) { n <- nd["n", g]; C <- rnorm(n); U <- rnorm(n); Sx <- rho * C + sqrt(1 - rho^2) * U
    lp <- par[1] * C + par[2] * U; a <- uniroot(function(a) mean(plogis(a + lp)) - nd["d", g] / n, c(-10, 10))$root
    y <- rbinom(n, 1, plogis(a + lp)); res <- resid(lm(Sx ~ C)); list(s = Sx, r = res, y = y) })
  est <- function(ix) { raw <- sapply(seq_along(dat), function(k) delong(dat[[k]]$s[ix[[k]]], dat[[k]]$y[ix[[k]]]))
    adj <- sapply(seq_along(dat), function(k) delong(dat[[k]]$r[ix[[k]]], dat[[k]]$y[ix[[k]]]))
    a0 <- pool_fixed(raw[1, ], raw[2, ]); a1 <- pool_fixed(adj[1, ], adj[2, ]); c(a0, a1, (a1 - 0.5) / (a0 - 0.5), a0 - a1) }
  full_ix <- lapply(dat, function(z) seq_along(z$y)); e0 <- est(full_ix)
  bt <- replicate(nb, est(lapply(dat, function(z) c(sample(which(z$y == 1), replace = TRUE), sample(which(z$y == 0), replace = TRUE)))))
  c(rf = e0[3], rf_hi = quantile(bt[3, ], 0.975, names = FALSE), diff = e0[4], diff_lo = quantile(bt[4, ], 0.025, names = FALSE),
    half_lo = quantile(bt[4, ] - 0.5 * (bt[1, ] - 0.5), 0.025, names = FALSE))
}
h3 <- list()
for (rho in c(0.3, 0.5, 0.7)) for (rf in c(0.2, 0.35, 0.5, 0.7)) {
  par <- calib(rho, rf); r <- t(replicate(300, one_rep(par, rho)))
  h3[[length(h3) + 1]] <- data.table(rho_c = rho, rf_true = rf, bC = par[1], bU = par[2],
    rf_median = median(r[, "rf"]), rf_iqr = IQR(r[, "rf"]),
    power_rf_hi_lt_0.5 = mean(r[, "rf_hi"] < 0.5),
    power_diffform = mean(r[, "half_lo"] > 0))
  message(sprintf("H3 rho=%.1f rf=%.2f done", rho, rf))
}
h3 <- rbindlist(h3); fwrite(h3, file.path(od, "H3_power.csv"))

ns <- c(GSE185263_ICU = 80, GSE236713 = 124, GSE95233 = 51, GSE54514 = 35, GSE65682 = 479)
h4 <- list()
for (p in c(0.12, 0.15, 0.20, 0.30)) for (tau in c(0.3, 0.6)) for (thr in c(0.05, 0.10)) {
  hit <- replicate(1000, { pi <- plogis(qlogis(p) + rnorm(length(ns), 0, tau)); x <- rbinom(length(ns), ns, pi)
    f <- tryCatch(rma(measure = "PLO", xi = x, ni = ns, method = "REML", test = "knha"), error = function(e) NULL)
    if (is.null(f)) NA else plogis(f$ci.lb) > thr })
  h4[[length(h4) + 1]] <- data.table(p_true = p, tau = tau, threshold = thr, power = mean(hit, na.rm = TRUE), failed = mean(is.na(hit)))
}
h4 <- rbindlist(h4); fwrite(h4, file.path(od, "H4_power.csv"))

write_json(list(date = as.character(Sys.Date()), seed = 20261116, B = B, n_models = length(ids),
                median_published = median(p0), SE_median_published = SE_pub, SE_single_model = SE_single,
                split_half_rho_median = median(split_rho), split_half_rho_q25 = unname(quantile(split_rho, 0.25)),
                replication_se = rep$se), file.path(od, "summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
writeLines(capture.output(sessionInfo()), file.path(od, "sessionInfo.txt"))
message("done")
