source("scripts/_common.R")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(scales) })
R <- file.path(PROJ, "results"); FIG <- file.path(PROJ, "figures"); TAB <- file.path(PROJ, "tables")
dir.create(FIG, showWarnings = FALSE); dir.create(TAB, showWarnings = FALSE)

C_PUB <- "#eb6834"; C_RND <- "#2a78d6"; C_AQ <- "#1baf7a"
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#8a8985"; GRID <- "#e6e5e1"
theme_pub <- function() theme_classic(base_size = 8, base_family = "Helvetica") +
  theme(axis.line = element_line(linewidth = 0.35, colour = INK2), axis.ticks = element_line(linewidth = 0.35, colour = INK2),
        axis.text = element_text(colour = INK2, size = 8), axis.title = element_text(colour = INK, size = 8),
        plot.title = element_text(colour = INK, face = "bold", size = 8.5, margin = margin(b = 3)),
        plot.tag = element_text(face = "bold", size = 10), legend.title = element_text(size = 8, colour = INK),
        legend.text = element_text(size = 8, colour = INK2), legend.key.size = unit(3.2, "mm"), legend.background = element_blank(),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 8, colour = INK, hjust = 0),
        panel.grid.major.x = element_line(colour = GRID, linewidth = 0.25), plot.margin = margin(4, 6, 4, 4))
save_fig <- function(p, name, h_mm, w_mm = 170, folder = NULL) {
  d <- file.path(FIG, "formal", if (is.null(folder)) name else folder); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  f <- function(ext) file.path(d, paste0(name, ".", ext))
  ggsave(f("pdf"), p, width = w_mm, height = h_mm, units = "mm", device = cairo_pdf)
  ggsave(f("eps"), p, width = w_mm, height = h_mm, units = "mm", device = cairo_ps, fallback_resolution = 600)
  ggsave(f("tiff"), p, width = w_mm, height = h_mm, units = "mm", dpi = 600, device = "tiff", compression = "lzw", type = "cairo", bg = "white")
  ggsave(f("png"), p, width = w_mm, height = h_mm, units = "mm", dpi = 600, bg = "white")
  message("saved ", d)
}
fmt_ci <- function(e, l, h, d = 3) sprintf(paste0("%.", d, "f (%.", d, "f to %.", d, "f)"), e, l, h)

M <- load_models(); pub_ids <- M$mt$model_id
SID <- setNames(sprintf("S%02d", seq_along(pub_ids)), pub_ids)
bm <- readRDS(file.path(R, "part1/bootstrap_main.rds")); pm <- bm$st0$pooled
bc <- readRDS(file.path(R, "cmaise/bootstrap_cmaise.rds")); pc <- bc$st0$pooled
ktype <- function(k) ifelse(grepl("^PUB_", k), "PUB", ifelse(grepl("^R2b_", k), "R2b", ifelse(grepl("^R1_", k), "R1", ifelse(grepl("^R3_", k), "R3",
                     ifelse(grepl("_(lasso|stepwise)$", k), "R2var", "R2")))))
r2_match <- function(k) sub("_[0-9]{3}$", "", sub("^R2_", "", k))
est_m <- fread(file.path(R, "part1/part1_estimates.csv")); est_c <- fread(file.path(R, "cmaise/cmaise_estimates.csv"))
sens_c <- fread(file.path(R, "cmaise/cmaise_sensitivity.csv")); s5 <- fread(file.path(R, "part1/S5_random_effects.csv"))
g_est <- function(d, an, st, sub = "all") { x <- d[analysis == an & subset == sub & stat == st]; c(x$estimate, x$lo, x$hi) }
c_est <- function(st, sub = "all") { x <- est_c[subset == sub & stat == st]; c(x$estimate, x$lo, x$hi) }

dA <- data.table(kk = names(pm), auc = pm)[, type := ktype(kk)][type %in% c("PUB", "R2", "R3")]
r1m <- data.table(kk = names(bm$st0$r1_med), auc = bm$st0$r1_med, type = "R1")
lev <- c("R1", "R3", "R2", "PUB")
lab <- c(PUB = "Published signatures (n = 36)", R2 = "Random-theme pipeline (n = 7,200)",
         R3 = "Whole-gene-set pipeline\n(200 runs, %d distinct signatures)", R1 = "Random genes, no selection\n(matched medians, n = 36)")
n_r3 <- uniqueN(round(dA[type == "R3"]$auc, 10)); lab[["R3"]] <- sprintf(lab[["R3"]], n_r3)
dA <- rbind(dA, r1m); dA[, row := factor(type, lev, lab[lev])]
sumA <- rbindlist(lapply(c("PUB", "R2", "R3", "R1"), function(t) { s <- switch(t, PUB = "pub_median", R2 = "E1", R3 = "R3_median", R1 = "E1_R1")
  v <- g_est(est_m, "main", s); data.table(type = t, est = v[1], lo = v[2], hi = v[3]) }))[, row := factor(type, lev, lab[lev])]
colA <- c(PUB = C_PUB, R2 = C_RND, R3 = C_AQ, R1 = MUTED)
pA <- ggplot(dA, aes(auc, row)) +
  geom_vline(xintercept = 0.5, linetype = "22", linewidth = 0.3, colour = INK2) +
  geom_violin(data = dA[type == "R2"], aes(fill = type), colour = NA, alpha = 0.35, scale = "width", width = 0.8) +
  geom_point(data = unique(dA[type == "R3"][, auc := round(auc, 10)], by = "auc"), aes(colour = type), position = position_jitter(height = 0.13, width = 0, seed = 2), size = 0.9, alpha = 0.8) +
  geom_point(data = dA[type %in% c("PUB", "R1")], aes(colour = type), position = position_jitter(height = 0.13, width = 0, seed = 1), size = 0.9, alpha = 0.8) +
  geom_errorbarh(data = sumA, aes(x = est, xmin = lo, xmax = hi), height = 0, linewidth = 0.6, colour = INK) +
  geom_point(data = sumA, aes(x = est), shape = 21, fill = "white", colour = INK, size = 1.8, stroke = 0.5) +
  geom_text(data = sumA, aes(x = 0.86, label = fmt_ci(est, lo, hi)), hjust = 1, size = 2.85, colour = INK2) +
  scale_fill_manual(values = colA, guide = "none") + scale_colour_manual(values = colA, guide = "none") +
  scale_x_continuous(limits = c(0.38, 0.86), breaks = seq(0.4, 0.8, 0.1)) + scale_y_discrete(limits = unname(lab[lev])) +
  labs(x = "Pooled AUROC for death, three main evaluation cohorts", y = NULL, title = "Discrimination of published and random signatures") + theme_pub()

forest <- rbindlist(list(
  data.table(grp = "Main analysis", lab = "Three main cohorts (36 signatures)", t(g_est(est_m, "main", "Delta"))),
  data.table(grp = "Pre-specified sensitivity", lab = "Trained on GSE65682 only (32)", t(g_est(est_m, "main", "Delta", "trained_GSE65682"))),
  data.table(grp = "Pre-specified sensitivity", lab = "Excluding 9 flagged signatures (27)", t(g_est(est_m, "main", "Delta", "exclude_flagged"))),
  data.table(grp = "Pre-specified sensitivity", lab = "Random weights from LASSO", t(g_est(est_m, "main", "Delta_R2_lasso"))),
  data.table(grp = "Pre-specified sensitivity", lab = "Random weights from stepwise Cox", t(g_est(est_m, "main", "Delta_R2_stepwise"))),
  data.table(grp = "Pre-specified sensitivity", lab = "Random pipeline without DEG step", t(g_est(est_m, "main", "Delta_R2_nodeg"))),
  data.table(grp = "Pre-specified sensitivity", lab = "Adding two sensitivity cohorts", t(g_est(est_m, "S4_main_plus_sensitivity_cohorts", "Delta"))),
  data.table(grp = "Pre-specified sensitivity", lab = "Random-effects pooling (point)", V1 = s5[stat == "Delta_RE"]$estimate, V2 = NA, V3 = NA),
  data.table(grp = "Independent validation", lab = "GSE272769 (USA, ICU)", t(g_est(est_m, "replication_GSE272769", "Delta"))),
  data.table(grp = "Independent validation", lab = "CMAISE (China, 32 centres)", t(c_est("Delta"))),
  data.table(grp = "Independent validation", lab = "CMAISE, within-centre pairs", t(unlist(sens_c[analysis == "centre_stratified" & stat == "Delta", .(estimate, lo, hi)])))
), use.names = FALSE)
setnames(forest, c("grp", "lab", "est", "lo", "hi"))
GRP <- c("Main analysis" = "Main", "Pre-specified sensitivity" = "Sensitivity", "Independent validation" = "Validation")
forest[, y := factor(lab, rev(lab))][, grp := factor(GRP[grp], GRP)]
mtxt <- data.table(grp = factor("Main", GRP), x = 0.054, label = "Margin\n0.05")
pB <- ggplot(forest, aes(est, y)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = MUTED) +
  geom_vline(xintercept = 0.05, linetype = "22", linewidth = 0.4, colour = INK) +
  geom_text(data = mtxt, aes(x = x, y = Inf, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1.15, size = 2.85, colour = INK2, lineheight = 0.9) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.45, colour = INK, na.rm = TRUE) +
  geom_point(aes(shape = grp), size = 1.6, fill = INK, colour = INK) +
  scale_shape_manual(values = c(Main = 18, Sensitivity = 16, Validation = 15), guide = "none") +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(-0.08, 0.10), breaks = seq(-0.08, 0.1, 0.04)) +
  labs(x = "Difference in median AUROC\n(published − random-theme pipeline)", y = NULL, title = "Published minus random signatures") +
  theme_pub() + theme(strip.text.y = element_text(angle = 270, hjust = 0.5, face = "plain", size = 8, colour = INK2))

ladder <- rbindlist(lapply(list(list("Main cohorts\n(3 GEO datasets)", "m", "main"), list("GSE272769\n(USA)", "m", "replication_GSE272769"), list("CMAISE\n(China)", "c", NA)), function(z) {
  rbindlist(lapply(c(R1 = "E1_R1", PUB = "pub_median", R2 = "E1", R3 = "R3_median"), function(s) {
    v <- if (z[[2]] == "m") g_est(est_m, z[[3]], s) else c_est(s); data.table(set = z[[1]], stat = s, est = v[1], lo = v[2], hi = v[3]) }))[, type := c("R1", "PUB", "R2", "R3")] }))
ladder[, set := factor(set, rev(unique(set)))][, type := factor(type, c("R1", "PUB", "R2", "R3"))]
labL <- c(R1 = "Random genes, no selection", PUB = "Published", R2 = "Random-theme pipeline", R3 = "Whole-gene-set pipeline")
pC <- ggplot(ladder, aes(est, set, colour = type)) +
  geom_vline(xintercept = 0.5, linetype = "22", linewidth = 0.3, colour = INK2) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.45, position = position_dodge(width = 0.6)) +
  geom_point(size = 1.6, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c(R1 = MUTED, PUB = C_PUB, R2 = C_RND, R3 = C_AQ), labels = labL, name = NULL) +
  scale_x_continuous(limits = c(0.48, 0.78), breaks = seq(0.5, 0.7, 0.1)) +
  guides(colour = guide_legend(ncol = 1, reverse = TRUE)) +
  labs(x = "Median pooled AUROC (95% CI)", y = NULL, title = "Medians in three data sets") +
  theme_pub() + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
fig2 <- free(pA) / ((pB | pC) + plot_layout(widths = c(1.25, 1))) + plot_layout(heights = c(1, 1.5)) + plot_annotation(tag_levels = "a")
save_fig(fig2, "Fig2_part1_published_vs_random", h_mm = 190, folder = "02_Fig2_part1_published_vs_random")

per_sig <- function(p) { d <- data.table(kk = names(p), auc = p)[ktype(kk) == "R2"][, model_id := r2_match(kk)]
  d[, .(lo = quantile(auc, 0.025, na.rm = TRUE), q1 = quantile(auc, 0.25, na.rm = TRUE), med = median(auc, na.rm = TRUE),
        q3 = quantile(auc, 0.75, na.rm = TRUE), hi = quantile(auc, 0.975, na.rm = TRUE)), by = model_id] }
e3m <- fread(file.path(R, "part1/E3_percentiles.csv")); e3c <- fread(file.path(R, "cmaise/cmaise_E3_percentiles.csv"))
psM <- merge(per_sig(pm), e3m[, .(model_id, pub = pooled_auc, pct = percentile)], by = "model_id")
psC <- merge(per_sig(pc), e3c[, .(model_id, pub = pooled_auc, pct = percentile)], by = "model_id")
ord <- psM[order(pct, pub)]$model_id
cat_plot <- function(d, title, show_y = TRUE) {
  d <- copy(d)[, sid := factor(SID[model_id], SID[ord])][, top := !is.na(pct) & pct > 0.95]
  ggplot(d, aes(y = sid)) +
    geom_vline(xintercept = 0.5, linetype = "22", linewidth = 0.3, colour = INK2) +
    geom_segment(aes(x = lo, xend = hi, yend = sid), linewidth = 0.35, colour = "#9ec5f4") +
    geom_segment(aes(x = q1, xend = q3, yend = sid), linewidth = 1.3, colour = C_RND) +
    geom_point(aes(x = pub, shape = top), colour = C_PUB, fill = C_PUB, size = 1.3, na.rm = TRUE) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    scale_x_continuous(limits = c(0.38, 0.75), breaks = seq(0.4, 0.7, 0.1)) +
    labs(x = "Pooled AUROC", y = if (show_y) "Published signature (ordered by percentile in main cohorts)" else NULL, title = title) +
    theme_pub() + theme(panel.grid.major.y = element_blank(), axis.text.y = if (show_y) element_text(size = 7.5) else element_blank(),
                        axis.ticks.y = if (show_y) element_line() else element_blank())
}
p3a <- cat_plot(psM, "Main evaluation cohorts")
p3b <- cat_plot(psC, "CMAISE (China)", show_y = FALSE)
pp <- merge(psM[, .(model_id, main = pct)], psC[, .(model_id, cm = pct)], by = "model_id")[!is.na(cm)]
rho <- cor(pp$main, pp$cm, method = "spearman")
p3c <- ggplot(pp, aes(main, cm)) +
  annotate("rect", xmin = 0.95, xmax = 1.02, ymin = 0.95, ymax = 1.02, fill = GRID, alpha = 1) +
  geom_hline(yintercept = 0.95, linetype = "22", linewidth = 0.3, colour = INK2) + geom_vline(xintercept = 0.95, linetype = "22", linewidth = 0.3, colour = INK2) +
  geom_point(size = 1.2, colour = INK, alpha = 0.8) +
  geom_text(data = pp[cm > 0.95], aes(label = SID[model_id]), size = 2.85, nudge_x = -0.03, hjust = 1, colour = INK2) +
  geom_text(data = pp[main > 0.95 & cm <= 0.95], aes(label = SID[model_id]), size = 2.85, nudge_y = 0.07, colour = INK2) +

  coord_equal(xlim = c(0, 1.02), ylim = c(0, 1.02)) + scale_x_continuous(labels = percent) + scale_y_continuous(labels = percent) +
  labs(x = "Percentile among matched random\nsignatures, main cohorts", y = "Percentile, CMAISE", title = "Top 5% in both: none", subtitle = sprintf("Spearman ρ = %.2f, n = %d", rho, nrow(pp))) + theme_pub() + theme(plot.subtitle = element_text(size = 8, colour = INK2))
key3 <- ggplot() + theme_void() +
  annotate("segment", x = 0, xend = 0.6, y = 4.3, yend = 4.3, colour = "#9ec5f4", linewidth = 0.35) +
  annotate("segment", x = 0.15, xend = 0.45, y = 4.3, yend = 4.3, colour = C_RND, linewidth = 1.3) +
  annotate("text", x = 0.8, y = 4.3, label = "Matched random signatures:\nIQR (thick), 95% range (thin)", hjust = 0, size = 2.85, colour = INK2, lineheight = 0.9) +
  annotate("point", x = 0.3, y = 3.2, colour = C_PUB, size = 1.6) +
  annotate("text", x = 0.8, y = 3.2, label = "Published signature", hjust = 0, size = 2.85, colour = INK2) +
  annotate("point", x = 0.3, y = 2.3, colour = C_PUB, shape = 17, size = 1.6) +
  annotate("text", x = 0.8, y = 2.3, label = "Published, above the\n95th percentile", hjust = 0, size = 2.85, colour = INK2, lineheight = 0.9) +
  annotate("text", x = 0, y = 1.1, label = "S02 not scored in CMAISE\n(gene coverage 78% < 80%)", hjust = 0, size = 2.85, colour = MUTED, lineheight = 0.9) +
  coord_cartesian(xlim = c(0, 6), ylim = c(0.4, 5), clip = "off")
fig3 <- (p3a | p3b | (p3c / key3 + plot_layout(heights = c(2.2, 1)))) + plot_layout(widths = c(1.1, 0.8, 1.25)) + plot_annotation(tag_levels = list(c("a", "b", "c", "")))
save_fig(fig3, "Fig3_per_signature_percentiles", h_mm = 150, folder = "03_Fig3_per_signature_percentiles")

pd <- fread(file.path(DATA, "cohorts/pdata_GSE54514.csv")); cg <- fread(file.path(R, "composition/composition_GSE54514.csv"))
qg <- merge(cg[, .(sample, group, neutrophil)], data.table(sample = pd$geo_accession, measured = suppressWarnings(as.numeric(pd[["neutrophil proportion"]]))), by = "sample")[!is.na(measured)]
rq <- c(all = cor(qg$neutrophil, qg$measured, method = "spearman"), sepsis = qg[group == "sepsis", cor(neutrophil, measured, method = "spearman")],
        control = qg[group == "control", cor(neutrophil, measured, method = "spearman")])
p4a <- ggplot(qg, aes(measured, neutrophil, colour = group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.3, colour = INK2) +
  geom_point(size = 1.3, alpha = 0.85) +
  scale_colour_manual(values = c(sepsis = C_RND, control = C_AQ), labels = c(sepsis = sprintf("Sepsis (n = %d)", sum(qg$group == "sepsis")), control = sprintf("Healthy (n = %d)", sum(qg$group == "control"))), name = NULL) +
  annotate("text", x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 2.85, colour = INK2,
           label = sprintf("Spearman ρ: all %.2f\nsepsis %.2f, healthy %.2f", rq["all"], rq["sepsis"], rq["control"])) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + scale_x_continuous(labels = percent) + scale_y_continuous(labels = percent) +
  labs(x = "Measured neutrophil fraction (blood count)", y = "Estimated neutrophil fraction (LM22)", title = "Composition estimates vs blood counts (GSE54514)") +
  theme_pub() + theme(legend.position = "inside", legend.position.inside = c(0.22, 0.70))
sv <- fread(file.path(R, "part2/score_vs_composition.csv"))[cohort %in% PART2_MAIN]
NCOH <- c(GSE185263 = 345, GSE63042 = 106); adjr2 <- function(r2, n, p = 5) 1 - (1 - r2) * (n - 1) / (n - p - 1)
sv[, n := NCOH[cohort]][, `:=`(lm22 = adjr2(r2_comp, n), marker = adjr2(r2_marker, n))]
r2s <- sv[, .(lm22 = sum(lm22 * n) / sum(n), marker = sum(marker * n) / sum(n)), by = model_id]
h3 <- fread(file.path(R, "part2/H3_result.csv"))
stopifnot(abs(median(r2s$lm22) - h3$median_adj_r2) < 1e-8)
r2l <- melt(r2s, id.vars = "model_id", variable.name = "method", value.name = "r2")[, method := factor(method, c("marker", "lm22"), c("Marker genes\n(descriptive)", "LM22 + NNLS\n(primary, H3)"))]
h3sum <- data.table(method = factor("LM22 + NNLS\n(primary, H3)", levels(r2l$method)), est = h3$median_adj_r2, lo = h3$lo, hi = h3$hi)
mksum <- r2l[method == levels(method)[1], .(est = median(r2))][, method := factor(levels(r2l$method)[1], levels(r2l$method))]
p4b <- ggplot(r2l, aes(r2, method)) +
  geom_vline(xintercept = 0.5, linetype = "22", linewidth = 0.4, colour = INK) +
  annotate("text", x = 0.51, y = Inf, label = "H3 threshold", hjust = 0, vjust = 1.3, size = 2.85, colour = INK2) +
  geom_point(position = position_jitter(height = 0.15, width = 0, seed = 2), size = 0.9, alpha = 0.75, colour = MUTED) +
  geom_errorbarh(data = h3sum, aes(x = est, xmin = lo, xmax = hi), height = 0, linewidth = 0.6, colour = INK) +
  geom_point(data = rbind(h3sum[, .(method, est)], mksum), aes(x = est), shape = 21, fill = "white", colour = INK, size = 1.8, stroke = 0.5) +
  scale_x_continuous(limits = c(-0.05, 0.85), breaks = seq(0, 0.8, 0.2)) +
  labs(x = "Adjusted R² of signature score on five leukocyte groups\n(weighted over GSE185263 and GSE63042)", y = NULL,
       title = "Variation explained by composition") + theme_pub() + theme(plot.tag.position = c(0.02, 0.985))
h3e <- fread(file.path(R, "part2/H3_estimates.csv")); p2c <- fread(file.path(R, "cmaise/cmaise_part2_descriptive.csv"))
rma <- fread(file.path(R, "part2/random_models_adjusted.csv")); sofa <- fread(file.path(R, "part2/sofa_adjusted_GSE185263.csv"))
tnl <- fread(file.path(R, "part2/tnlr_auroc.csv"))
gv <- function(an, st) h3e[analysis == an & stat == st]$estimate; cv <- function(it) p2c[item == it]$value
db <- data.table(lab = c("Published, main cohorts", "Published, GSE272769", "Published, CMAISE", "Random pipeline, main cohorts", "Random pipeline, CMAISE", "SOFA, GSE185263"),
                 raw = c(gv("main", "median_raw"), gv("replication", "median_raw"), cv("published_median_raw"), median(rma$raw), cv("random_median_raw"), sofa[score == "SOFA"]$auc),
                 adj = c(gv("main", "median_adj"), gv("replication", "median_adj"), cv("published_median_adj"), median(rma$adj), cv("random_median_adj"), sofa[score == "SOFA_adjusted"]$auc))
db[, lab := factor(lab, rev(lab))]
dt <- data.table(lab = c("Transcriptomic NLR alone, main cohorts", "Transcriptomic NLR alone, GSE272769", "Transcriptomic NLR alone, CMAISE"),
                 auc = c(tnl[cohort == "pooled_main"]$auc, tnl[cohort == "GSE272769"]$auc, cv("tnlr_pooled_auc")))
dbl <- rbind(melt(db, id.vars = "lab", variable.name = "kind", value.name = "auc"), dt[, .(lab, kind = "tnlr", auc)])
TOP4 <- c("Published, main cohorts", "Published, GSE272769", "Published, CMAISE", "Random pipeline, main cohorts", "Random pipeline, CMAISE", "SOFA, GSE185263", dt$lab)
LV4 <- rev(TOP4)
dbl[, lab := factor(as.character(lab), LV4)]; db[, lab := factor(as.character(lab), LV4)]
p4c <- ggplot(dbl, aes(auc, lab)) +
  geom_vline(xintercept = 0.5, linetype = "22", linewidth = 0.3, colour = INK2) +
  geom_segment(data = db, aes(x = raw, xend = adj, y = lab, yend = lab), linewidth = 0.6, colour = MUTED) +
  geom_point(aes(shape = kind, colour = kind, fill = kind), size = 1.7) +
  scale_shape_manual(values = c(raw = 21, adj = 21, tnlr = 23), labels = c(raw = "Unadjusted", adj = "Adjusted for composition", tnlr = "Composition only (NLR)"), name = NULL) +
  scale_colour_manual(values = c(raw = INK, adj = C_RND, tnlr = INK2), labels = c(raw = "Unadjusted", adj = "Adjusted for composition", tnlr = "Composition only (NLR)"), name = NULL) +
  scale_fill_manual(values = c(raw = "white", adj = C_RND, tnlr = GRID), labels = c(raw = "Unadjusted", adj = "Adjusted for composition", tnlr = "Composition only (NLR)"), name = NULL) +
  scale_x_continuous(limits = c(0.45, 0.72), breaks = seq(0.45, 0.7, 0.05)) + scale_y_discrete(limits = LV4) +
  labs(x = "Pooled AUROC for death (medians for signatures)", y = NULL, title = "Prognostic information remaining after adjustment") +
  theme_pub() + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
fig4 <- free((p4a | p4b) + plot_layout(widths = c(1, 1.25))) / p4c + plot_layout(heights = c(1.1, 1)) + plot_annotation(tag_levels = "a")
save_fig(fig4, "Fig4_leukocyte_composition", h_mm = 175, folder = "04_Fig4_leukocyte_composition")

bcmp <- fread(file.path(R, "batch/batch_comparisons.csv")); ldist <- fread(file.path(R, "batch/label_distribution.csv")); h4 <- fread(file.path(R, "batch/H4_result.csv"))
TOOLS <- data.table(tool = c("SRS_davenport", "SRS_extended", "SRS_davenport_k20", "CTS", "SUBSPACE_quadrant_healthyref", "SUBSPACE_quadrant_median"),
                    lab = c("SRS, 7 genes", "SRS, 19 genes", "SRS, 7 genes, k = 20", "CTS", "SUBSPACE", "SUBSPACE, median split"),
                    ref = c(rep("External\nreference", 4), "Batch\ncontrols", "Batch\npatients"))
SCEN <- c(`C0 vs C1` = "Add\nhealthy\ncontrols", `C0 vs C2` = "Mix\nICU and\nED", `C3_with_high vs C3_with_low` = "Compo-\nsition\ntitration",
          C5 = "Sub-\nsampled\ncontrols", seed = "Repeat\nrun\n(seed)")
hm <- bcmp[tool %in% TOOLS$tool & !is.na(prop_changed) & !grepl("^seed|C4", contrast), .(cohort, tool, contrast, p = prop_changed)]
seedd <- bcmp[tool == "CTS" & grepl("^seed", contrast), .(p = median(prop_changed)), by = cohort][, `:=`(tool = "CTS", contrast = "seed")]
c5 <- rbindlist(lapply(list.files(file.path(R, "batch"), "^C5_subspace", full.names = TRUE), fread))[, .(p = median(change)), by = cohort][, `:=`(tool = "SUBSPACE_quadrant_healthyref", contrast = "C5")]
hm <- rbind(hm, seedd, c5, use.names = TRUE)
hs <- hm[, .(med = median(p), mn = min(p), mx = max(p), k = .N), by = .(tool, contrast)]
hs <- merge(hs, TOOLS, by = "tool")[, scen := factor(SCEN[contrast], SCEN)]
hs[, row := factor(paste0(lab, "  "), rev(unique(paste0(TOOLS$lab, "  "))))][, ref := factor(ref, unique(TOOLS$ref))]
hs[, txt := ifelse(k > 1, sprintf("%.0f%%\n(%.0f–%.0f)", 100 * med, 100 * mn, 100 * mx), sprintf("%.0f%%", 100 * med))]
p5a <- ggplot(hs, aes(scen, row, fill = med)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = txt, colour = med > 0.3), size = 2.85, lineheight = 0.9) +
  scale_fill_gradient(low = "#cde2fb", high = "#184f95", limits = c(0, 0.6), oob = squish, labels = percent, name = "Patients\nre-labelled") +
  scale_colour_manual(values = c(`FALSE` = INK, `TRUE` = "white"), guide = "none") +
  facet_grid(ref ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = NULL, y = NULL, title = "Patients re-labelled: median across datasets (range)") +
  theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
                      strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "plain", size = 8, colour = INK2),
                      legend.key.height = unit(5, "mm"))
srs1 <- bcmp[tool == "SRS_davenport" & contrast == "C0 vs C1", .(cohort, n, changed)]
CTRY <- c(GSE65682 = "NL", GSE185263 = "CA", GSE95233 = "FR", GSE54514 = "AU", GSE236713 = "UK")
srs1[, c("est", "lo", "hi") := as.list(c(changed / n, binom.test(changed, n)$conf.int)), by = cohort]
srs1[, lab := sprintf("%s (%s), %d/%d", cohort, CTRY[cohort], changed, n)]
srs1 <- rbind(srs1[order(est)], data.table(cohort = "pooled", n = h4$n, changed = h4$changed, est = h4$pooled, lo = h4$lo, hi = h4$hi,
                                            lab = sprintf("Pooled (random effects), %d/%d", h4$changed, h4$n)), fill = TRUE)
srs1[, y := factor(lab, rev(lab))][, pooled := cohort == "pooled"]
p5b <- ggplot(srs1, aes(est, y)) +
  geom_vline(xintercept = 0.05, linetype = "22", linewidth = 0.4, colour = INK) +
  annotate("text", x = 0.055, y = Inf, label = "H4 threshold (5%)", hjust = 0, vjust = 1.2, size = 2.85, colour = INK2) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.45, colour = INK) +
  geom_point(aes(shape = pooled), size = 1.8, colour = INK, fill = INK) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
  scale_x_continuous(labels = percent, limits = c(0, 0.36), breaks = seq(0, 0.3, 0.1)) +
  labs(x = "Patients re-labelled when\nhealthy controls join (%)", y = NULL, title = "H4: SRS, 7 genes") + theme_pub()
hc <- ldist[tool %in% c("CTS", "SRS_davenport") & scenario %in% c("C4_healthy_alone", "C1"), .(n = sum(n)), by = .(tool, scenario, value)]
hc[, value := sub("^SRS", "", value)][, cls := paste0("Type ", value)][, share := n / sum(n), by = .(tool, scenario)]
hc[, bar := factor(paste(ifelse(tool == "CTS", "CTS", "SRS"), ifelse(scenario == "C1", "with patients", "alone"), sep = ", "),
                   rev(c("CTS, alone", "CTS, with patients", "SRS, alone", "SRS, with patients")))]
hcn <- hc[, .(N = sum(n)), by = bar]
p5c <- ggplot(hc, aes(share, bar, fill = cls)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.4, position = position_stack(reverse = TRUE)) +
  geom_text(data = hc[share >= 0.15], aes(label = percent(share, 1)), position = position_stack(vjust = 0.5, reverse = TRUE), size = 2.85, colour = "white") +
  geom_text(data = hcn, aes(x = 1.03, y = bar, label = paste0("n = ", N)), inherit.aes = FALSE, hjust = 0, size = 2.85, colour = INK2) +
  scale_fill_manual(values = c(`Type 1` = C_RND, `Type 2` = C_PUB, `Type 3` = C_AQ), name = NULL) +
  scale_x_continuous(labels = percent, limits = c(0, 1.32), breaks = seq(0, 1, 0.5)) +
  labs(x = "Healthy controls\nassigned to each type (%)", y = NULL, title = "Same healthy controls") +
  theme_pub() + theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
CONT <- data.table(tool = c("Mars1_ratio", "Sweeney2018NC_Duke18_raw", "Sweeney2018NC_Stanford12_raw", "SUBSPACE_myeloid_raw", "SUBSPACE_lymphoid_raw",
                            "SRSq_davenport", "SRSq_extended", "Sweeney2018NC_Duke18_batchz", "Sweeney2018NC_Stanford12_batchz",
                            "Jiang2022_IRRG9_batchz", "Chen2023_IRG7_batchz", "Feng2024_ROS37_batchz", "Li2026_GMGscore6_batchz"),
                   lab = c("Mars1 ratio", "Sweeney Duke, raw", "Sweeney Stanford, raw", "SUBSPACE myeloid, raw", "SUBSPACE lymphoid, raw",
                           "SRSq, 7 genes", "SRSq, 19 genes", "Sweeney Duke, batch z", "Sweeney Stanford, batch z",
                           paste0(SID[c("Jiang2022_IRRG9", "Chen2023_IRG7", "Feng2024_ROS37", "Li2026_GMGscore6")], ", batch z")),
                   grp = c(rep("Single sample", 5), rep("External reference", 2), rep("Batch z score", 6)))
cs <- merge(bcmp[contrast == "C0 vs C1" & tool %in% CONT$tool, .(cohort, tool, shift = median_abs_shift_sd)], CONT, by = "tool")
cs[, lab := factor(lab, rev(CONT$lab))][, grp := factor(grp, unique(CONT$grp))]
p5d <- ggplot(cs, aes(shift, lab)) +
  geom_point(size = 1, colour = MUTED, alpha = 0.9) +
  stat_summary(fun = median, geom = "point", shape = 124, size = 3, colour = INK) +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(limits = c(-0.02, 0.7), breaks = seq(0, 0.7, 0.1)) +
  labs(x = "Median shift of patients' scores (SD) when healthy controls join", y = NULL, title = "Continuous scores: level shift (dots: datasets, bar: median)") +
  theme_pub() + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, face = "plain", size = 8, colour = INK2))
fig5 <- free(p5a) / free((p5b | p5c) + plot_layout(widths = c(1, 1.05))) + plot_layout(heights = c(1.2, 1)) + plot_annotation(tag_levels = "a")
save_fig(fig5, "Fig5_batch_composition", h_mm = 175, folder = "05_Fig5_batch_composition")
save_fig(p5d + labs(title = NULL), "FigS1_continuous_score_shift", h_mm = 105, folder = "S1_FigS1_continuous_score_shift")

hyp <- fread(file.path(R, "hypotheses_final.csv")); Lc <- fread(file.path(R, "cmaise/cmaise_L_result.csv"))
t2 <- rbind(
  hyp[, .(item = hypothesis, criterion = rule, estimate = fmt_ci(estimate, lo, hi), p_one_sided = ifelse(is.na(p_one_sided), "", format.pval(p_one_sided, digits = 2)),
          p_holm = ifelse(is.na(p_holm), "", format.pval(p_holm, digits = 2)), decision = ifelse(hypothesis == "H1", ifelse(ci_rule_met, "Supported", "Not supported"),
          ifelse(supported_after_holm, "Supported", "Not supported")))],
  Lc[item %in% c("L1", "L2"), .(item = paste(item, "(CMAISE, added after main results)"), criterion = rule, estimate = fmt_ci(estimate, lo, hi), p_one_sided = "", p_holm = "",
                                decision = ifelse(met, "Met", "Not met"))])
fwrite(t2, file.path(TAB, "Table2_hypotheses.csv"))
pe <- fread(file.path(PROJ, "docs/pipeline_extraction.csv"))[, .(model_id, theme, theme_size, training_cohort)]
s1 <- data.table(id = SID[pub_ids], model_id = pub_ids, pmid = M$mt$pmid, n_genes = sapply(pub_ids, function(id) nrow(M$gn[model_id == id][!duplicated(gene)])))
s1 <- merge(s1, pe, by = "model_id", all.x = TRUE)
s1 <- merge(s1, psM[, .(model_id, main_auc = round(pub, 3), main_random_median = round(med, 3), main_percentile = round(pct, 3))], by = "model_id", all.x = TRUE)
s1 <- merge(s1, psC[, .(model_id, cmaise_auc = round(pub, 3), cmaise_random_median = round(med, 3), cmaise_percentile = round(pct, 3))], by = "model_id", all.x = TRUE)
s1[, flagged := model_id %in% FLAGGED]; setorder(s1, id)
fwrite(s1, file.path(TAB, "TableS1_signature_key.csv"))

t1 <- data.table(
  dataset = c("GSE65682 (MARS)", "GSE185263", "GSE63042 (CAPSOD)", "GSE236713", "GSE272769", "GSE95233", "GSE54514", "CMAISE"),
  country = c("Netherlands", "Canada (ICU); 4 countries (ED)", "USA", "UK", "USA", "France", "Australia", "China (32 centres)"),
  setting = c("ICU", "ED and ICU", "ED", "ICU", "ICU", "ICU, septic shock", "ICU", "Hospitalised sepsis, 4 severity strata"),
  platform = c("Affymetrix array, whole blood", "RNA-seq, whole blood", "RNA-seq, whole blood", "Agilent array, leukocytes", "Affymetrix array, whole blood",
               "Affymetrix array, whole blood", "Illumina array, whole blood", "RNA-seq, whole blood"),
  outcome = c("Death by day 28", "In-hospital death", "Death by day 28 (sampled by outcome)", "Death (window not reported)", "Death within 30 days",
              "Death by day 28", "Death (window not reported)", "In-hospital death"),
  patients = c(479, 345, 106, 124, 161, 51, 35, 355), deaths = c(114, 52, 28, 30, 60, 17, 9, 47), healthy_controls = c(42, 44, 0, 30, 0, 20, 18, 0),
  part1 = c("Training of random signatures", "Main evaluation", "Main evaluation", "Main evaluation", "Replication", "Sensitivity", "Sensitivity", "Local validation (added after main results)"),
  part2 = c("—", "Main", "Main", "Sensitivity only (gene-centred data)", "Replication", "—", "Accuracy check (measured neutrophils)", "Descriptive"),
  part3 = c("Index patients (not CTS)", "Index patients (ICU, 80); ED patients and controls as partners", "—", "Index patients", "—", "Index patients", "Index patients", "—"))
t1[, deaths := sprintf("%d (%.1f%%)", deaths, 100 * deaths / patients)]
fwrite(t1, file.path(TAB, "Table1_cohorts.csv"))

message("done")
