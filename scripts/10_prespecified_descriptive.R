source("scripts/_common.R")
set.seed(SEED + 10)
IN <- file.path(PROJ, "results")
od1 <- outdir("part1"); od3 <- outdir("batch")
B <- if (DRY_RUN) 20 else 2000
qci <- function(x) unname(quantile(x, c(0.025, 0.975), na.rm = TRUE))

RM <- readRDS(file.path(IN, "part1", "random_models.rds"))
r2 <- data.table(rkey = sapply(RM$R2, `[[`, "key"), matched = sapply(RM$R2, `[[`, "matched")); rm(RM)
p <- readRDS(file.path(IN, "part1", paste0("bootstrap_replication_", REPL, ".rds")))$st0$pooled
stopifnot(all(r2$rkey %in% names(p)))
if (DRY_RUN) p <- setNames(sample(p), names(p))
pub_ids <- sub("^PUB_", "", grep("^PUB_", names(p), value = TRUE))
E3r <- rbindlist(lapply(pub_ids, function(id) { a <- p[[paste0("PUB_", id)]]; r <- p[r2[matched == id]$rkey]
  data.table(model_id = id, pooled_auc = a, random_median = median(r, na.rm = TRUE), percentile = mean(r < a, na.rm = TRUE), n_random = sum(!is.na(r))) }))
fwrite(E3r, file.path(od1, paste0("E3_percentiles_replication_", REPL, ".csv")))
e3m <- fread(file.path(IN, "part1", "E3_percentiles.csv")); e3c <- fread(file.path(IN, "cmaise", "cmaise_E3_percentiles.csv"))
mm <- merge(merge(e3m[, .(model_id, main = percentile)], E3r[, .(model_id, replication = percentile)], by = "model_id"), e3c[, .(model_id, cmaise = percentile)], by = "model_id", all.x = TRUE)
mm[, n_gt95 := (main > 0.95) + (replication > 0.95) + ifelse(is.na(cmaise), 0, cmaise > 0.95)]
sumA <- data.table(item = c("n_scored", "n_gt95", "n_lt05", "n_below_random_median", "median_percentile", "spearman_vs_main", "spearman_vs_cmaise",
                            "n_gt95_main_and_replication", "n_gt95_in_two_or_more_of_three", "n_gt95_in_all_three"),
                   value = c(sum(!is.na(E3r$pooled_auc)), sum(E3r$percentile > 0.95, na.rm = TRUE), sum(E3r$percentile < 0.05, na.rm = TRUE),
                             sum(E3r$pooled_auc < E3r$random_median, na.rm = TRUE), median(E3r$percentile, na.rm = TRUE),
                             cor(mm$main, mm$replication, method = "spearman", use = "complete.obs"), cor(mm$replication, mm$cmaise, method = "spearman", use = "complete.obs"),
                             sum(mm$main > 0.95 & mm$replication > 0.95, na.rm = TRUE), sum(mm$n_gt95 >= 2, na.rm = TRUE), sum(mm$n_gt95 == 3, na.rm = TRUE)))
fwrite(sumA, file.path(od1, paste0("E3_replication_summary_", REPL, ".csv"))); print(sumA)

DSETS <- c("GSE185263", "GSE236713", "GSE95233", "GSE54514", "GSE65682")
res <- fread(file.path(IN, "batch", "batch_outputs_long.csv"), colClasses = list(character = "value"))[type == "label" & !grepl("_seed", tool)]
if (DRY_RUN) res[, value := sample(value), by = .(cohort, scenario, tool)]
death <- rbindlist(lapply(DSETS, function(g) { ph <- load_cohort(g)$ph; data.table(cohort = g, sample = ph$sample, death = as.integer(ph$death)) }))
lab <- function(g, sc, tl, ids) { r <- res[cohort == g & scenario == sc & tool == tl]; r$value[match(ids, r$sample)] }
index_ids <- function(g) unique(res[cohort == g & scenario == "C0"]$sample)
strat <- function(y) c(sample(which(y == 1), replace = TRUE), sample(which(y == 0), replace = TRUE))
TOOLS <- c("SRS_davenport", "SRS_extended", "SRS_davenport_k20", "CTS", "SUBSPACE_quadrant_median", "SUBSPACE_quadrant_healthyref")

ari <- function(a, b) { t <- table(a, b); n <- sum(t); if (n < 2) return(NA_real_)
  sc <- sum(choose(t, 2)); sa <- sum(choose(rowSums(t), 2)); sb <- sum(choose(colSums(t), 2)); e <- sa * sb / choose(n, 2); m <- (sa + sb) / 2
  if (isTRUE(all.equal(m, e))) NA_real_ else (sc - e) / (m - e) }
ariB <- list()
for (g in DSETS) { ids <- index_ids(g); y <- death[cohort == g][match(ids, sample)]$death; stopifnot(!anyNA(y))
  bidx <- replicate(B, strat(y), simplify = FALSE)
  for (sc in c("C0", "C1")) { L <- lapply(setNames(TOOLS, TOOLS), function(tl) lab(g, sc, tl, ids)); L <- L[!sapply(L, function(v) all(is.na(v)))]
    if (length(L) < 2) next
    for (pr in combn(names(L), 2, simplify = FALSE)) { a <- L[[pr[1]]]; b <- L[[pr[2]]]; stopifnot(!anyNA(a), !anyNA(b))
      bt <- sapply(bidx, function(i) ari(a[i], b[i])); ci <- qci(bt)
      ariB[[length(ariB) + 1]] <- data.table(cohort = g, scenario = sc, tool_a = pr[1], tool_b = pr[2], n = length(ids), ari = ari(a, b), lo = ci[1], hi = ci[2]) } } }
ariB <- rbindlist(ariB); fwrite(ariB, file.path(od3, "ari_pairwise.csv"))
ariS <- ariB[, .(k = .N, n_total = sum(n), ari_median = median(ari), ari_min = min(ari), ari_max = max(ari)), by = .(scenario, tool_a, tool_b)]
fwrite(ariS, file.path(od3, "ari_pairwise_summary.csv")); print(ariS, nrows = 100)

cramer_v <- function(l, y) { t <- table(l, y); t <- t[rowSums(t) > 0, colSums(t) > 0, drop = FALSE]; if (nrow(t) < 2 || ncol(t) < 2) return(NA_real_)
  n <- sum(t); ex <- outer(rowSums(t), colSums(t)) / n; sqrt(sum((t - ex)^2 / ex) / (n * (min(dim(t)) - 1))) }
byc <- list(); asc <- list()
for (g in DSETS) {
  mid <- { f <- file.path(IN, "batch", paste0("tertiles_", g, ".csv")); if (file.exists(f)) fread(f)[tnlr_tertile == "mid"]$sample else character(0) }
  CON <- list(list(s1 = "C0", s2 = "C1", ids = index_ids(g)), list(s1 = "C0", s2 = "C2", ids = index_ids(g)), list(s1 = "C3_with_high", s2 = "C3_with_low", ids = mid))
  for (cn in CON) for (tl in TOOLS) { ids <- cn$ids; if (!length(ids)) next
    a <- lab(g, cn$s1, tl, ids); b <- lab(g, cn$s2, tl, ids); if (all(is.na(a)) || all(is.na(b))) next
    stopifnot(!anyNA(a), !anyNA(b)); y <- death[cohort == g][match(ids, sample)]$death; stopifnot(!anyNA(y)); contrast <- paste(cn$s1, "vs", cn$s2)
    for (s in 1:2) { l <- if (s == 1) a else b
      byc[[length(byc) + 1]] <- data.table(cohort = g, contrast = contrast, tool = tl, scenario = c(cn$s1, cn$s2)[s], class = names(table(l)),
                                           n = as.integer(table(l)), deaths = as.integer(tapply(y, l, sum)[names(table(l))])) }
    bt <- replicate(B, { i <- strat(y); c(cramer_v(a[i], y[i]), cramer_v(b[i], y[i])) }); v1 <- cramer_v(a, y); v2 <- cramer_v(b, y)
    c1 <- qci(bt[1, ]); c2 <- qci(bt[2, ]); cd <- qci(bt[2, ] - bt[1, ])
    asc[[length(asc) + 1]] <- data.table(cohort = g, contrast = contrast, tool = tl, n = length(ids), deaths = sum(y), changed = sum(a != b),
                                         v_s1 = v1, v_s1_lo = c1[1], v_s1_hi = c1[2], v_s2 = v2, v_s2_lo = c2[1], v_s2_hi = c2[2],
                                         diff = v2 - v1, diff_lo = cd[1], diff_hi = cd[2]) } }
byc <- rbindlist(byc)[, mortality := deaths / n]; asc <- rbindlist(asc)
fwrite(byc, file.path(od3, "label_death_by_class.csv")); fwrite(asc, file.path(od3, "label_death_association.csv"))
ascS <- asc[, .(k = .N, k_defined = sum(!is.na(diff)), n_total = sum(n), deaths_total = sum(deaths), v_s1_median = median(v_s1, na.rm = TRUE), v_s2_median = median(v_s2, na.rm = TRUE),
                diff_median = median(diff, na.rm = TRUE), diff_min = min(diff, na.rm = TRUE), diff_max = max(diff, na.rm = TRUE),
                n_ci_excludes_0 = sum(diff_lo > 0 | diff_hi < 0, na.rm = TRUE)), by = .(contrast, tool)]
fwrite(ascS, file.path(od3, "label_death_association_summary.csv")); print(asc, nrows = 200); print(ascS, nrows = 100)
writeLines(capture.output(sessionInfo()), file.path(od3, "sessionInfo_10.txt"))
message("done")
