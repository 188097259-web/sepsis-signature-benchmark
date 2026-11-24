suppressPackageStartupMessages({ library(grid); library(data.table) })
PROJ <- Sys.getenv("SEPSIS_BENCH_ROOT", "."); R <- file.path(PROJ, "results")
OUTD <- file.path(PROJ, "figures/formal/GA_graphical_abstract"); dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
INK <- "#0b0b0b"; INK2 <- "#52514e"; LINE <- "#c9c7c1"; TINT <- "#f4f3ef"; C_PUB <- "#eb6834"; C_RND <- "#2a78d6"; C_AQ <- "#1baf7a"
em <- fread(file.path(R, "part1/part1_estimates.csv")); ec <- fread(file.path(R, "cmaise/cmaise_estimates.csv")); bc <- fread(file.path(R, "batch/batch_comparisons.csv"))
g <- function(an, st) em[analysis == an & subset == "all" & stat == st]$estimate
sets <- data.table(lab = c("3 main GEO cohorts", "USA (replication)", "China (CMAISE)"),
                   pub = c(g("main", "pub_median"), g("replication_GSE272769", "pub_median"), ec[subset == "all" & stat == "pub_median"]$estimate),
                   rnd = c(g("main", "E1"), g("replication_GSE272769", "E1"), ec[subset == "all" & stat == "E1"]$estimate))
med <- function(tl) median(bc[tool == tl & contrast == "C0 vs C1"]$prop_changed, na.rm = TRUE)
bars <- data.table(lab = c("Single-sample scores", "SRS (mNN reference)", "CTS (ComBat reference)"),
                   v = c(max(bc[tool %in% c("Mars1_ratio") & contrast == "C0 vs C1"]$median_abs_shift_sd, 0), med("SRS_davenport"), med("CTS")))
stopifnot(bars$v[1] == 0)
FONT <- "Helvetica"
t <- function(s, x, y, size = 11, face = "plain", col = INK, hj = 0, vj = 0.5) grid.text(s, unit(x, "npc"), unit(y, "npc"), just = c(hj, vj), gp = gpar(fontsize = size, fontface = face, col = col, fontfamily = FONT, lineheight = 1.05))
draw <- function() {
  grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))
  grid.roundrect(unit(0.012, "npc"), unit(0.2, "npc"), unit(0.27, "npc"), unit(0.76, "npc"), just = c("left", "bottom"), r = unit(6, "pt"), gp = gpar(fill = TINT, col = NA))
  t("Do 'X-related gene' sepsis\nprognostic signatures beat\nrandom genes?", 0.028, 0.80, size = 12.5, face = "bold", vj = 1)
  t("36 published signatures vs\n7,200 random signatures built\nwith a common pipeline", 0.028, 0.43, size = 10.5, col = INK2, vj = 0.5)
  x0 <- 0.295; t("Median AUROC for death", x0, 0.91, size = 11, face = "bold")
  sc <- function(v) x0 + 0.155 + (v - 0.5) / 0.2 * 0.13
  for (i in 1:3) { y <- 0.70 - (i - 1) * 0.17
    t(sets$lab[i], x0, y, size = 9.5, col = INK2)
    grid.segments(unit(sc(0.5), "npc"), unit(y, "npc"), unit(sc(0.7), "npc"), unit(y, "npc"), gp = gpar(col = LINE, lwd = 1))
    grid.points(unit(sc(sets$pub[i]), "npc"), unit(y, "npc"), pch = 16, size = unit(7, "pt"), gp = gpar(col = C_PUB))
    grid.points(unit(sc(sets$rnd[i]), "npc"), unit(y, "npc"), pch = 1, size = unit(8, "pt"), gp = gpar(col = C_RND, lwd = 2))
    t(sprintf("%.2f vs %.2f", sets$pub[i], sets$rnd[i]), sc(0.7) + 0.008, y, size = 9.5, col = INK2) }
  t("0.5", sc(0.5), 0.22, size = 8.5, col = INK2, hj = 0.5); t("0.7", sc(0.7), 0.22, size = 8.5, col = INK2, hj = 0.5)
  t("\u25cf published", x0, 0.81, size = 9.5, col = C_PUB); t("\u25cb random", x0 + 0.088, 0.81, size = 9.5, col = C_RND)
  x1 <- 0.70; t("Patients re-labelled when healthy\ncontrols join the batch (median)", x1, 0.88, size = 11, face = "bold", vj = 0.5)
  for (i in 1:3) { y <- 0.66 - (i - 1) * 0.17
    t(bars$lab[i], x1, y + 0.07, size = 9.5, col = INK2)
    grid.rect(unit(x1, "npc"), unit(y - 0.035, "npc"), unit(max(bars$v[i], 0.002) / 0.30 * 0.2, "npc"), unit(0.07, "npc"), just = c("left", "bottom"), gp = gpar(fill = if (i == 1) INK2 else if (i == 2) C_AQ else C_RND, col = NA))
    t(sprintf("%.0f%%", 100 * bars$v[i]), x1 + max(bars$v[i], 0.002) / 0.30 * 0.2 + 0.01, y, size = 10, face = "bold") }
  grid.segments(unit(0.012, "npc"), unit(0.165, "npc"), unit(0.988, "npc"), unit(0.165, "npc"), gp = gpar(col = LINE, lwd = 0.8))
  t("Published sepsis transcriptomic prognostic signatures were no better than random genes run through a common pipeline;\nbatch-dependent endotyping tools re-labelled patients when batch composition changed.",
    0.5, 0.075, size = 9.5, col = INK, hj = 0.5)
}
png(file.path(OUTD, "Graphical_abstract_920x300.png"), width = 920, height = 300, res = 100, type = "cairo"); draw(); invisible(dev.off())
png(file.path(OUTD, "Graphical_abstract_1840x600.png"), width = 1840, height = 600, res = 200, type = "cairo"); draw(); invisible(dev.off())
cairo_pdf(file.path(OUTD, "Graphical_abstract.pdf"), width = 9.2, height = 3, family = FONT); draw(); invisible(dev.off())
fwrite(rbind(sets[, .(item = lab, published = round(pub, 3), random = round(rnd, 3))]), file.path(OUTD, "ga_values_auroc.csv")); fwrite(bars[, .(item = lab, relabelled = round(v, 3))], file.path(OUTD, "ga_values_relabel.csv"))
message("graphical abstract saved")
