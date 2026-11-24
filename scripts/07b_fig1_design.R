suppressPackageStartupMessages({ library(grid); library(data.table) })
PROJ <- Sys.getenv("SEPSIS_BENCH_ROOT", ".")
OUTD <- file.path(PROJ, "figures/formal/01_Fig1_study_design"); dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
W <- 170; H <- 149
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#8a8985"; LINE <- "#c9c7c1"; TINT <- "#f4f3ef"
C_PUB <- "#eb6834"; C_RND <- "#2a78d6"; C_AQ <- "#1baf7a"
RAMP <- c("#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b")
FONT <- "Helvetica"
t1 <- fread(file.path(PROJ, "tables/Table1_cohorts.csv"))
nd <- function(ds) { r <- t1[startsWith(dataset, ds)]; sprintf("%s pts · %s deaths", format(r$patients, big.mark = ","), sub(" .*", "", r$deaths)) }

mm <- function(x) unit(x, "mm"); Y <- function(yt) H - yt
box <- function(x, yt, w, h, fill = "white", col = LINE, lwd = 0.6, r = 1.8, lty = 1)
  grid.roundrect(x = mm(x), y = mm(Y(yt)), width = mm(w), height = mm(h), just = c("left", "top"), r = mm(r),
                 gp = gpar(fill = fill, col = col, lwd = lwd, lty = lty))
txt <- function(s, x, yt, size = 8, face = "plain", col = INK, hjust = 0, vjust = 1, lh = 1.12)
  grid.text(s, x = mm(x), y = mm(Y(yt)), just = c(hjust, vjust), gp = gpar(fontsize = size, fontface = face, col = col, lineheight = lh, fontfamily = FONT))
arr <- function(x0, yt0, x1, yt1, col = INK2, lwd = 0.9)
  grid.segments(mm(x0), mm(Y(yt0)), mm(x1), mm(Y(yt1)), arrow = arrow(length = mm(1.5), type = "closed", angle = 24), gp = gpar(col = col, fill = col, lwd = lwd))
dot <- function(x, yt, r = 0.9, fill = INK2, col = NA, lwd = 0.6) grid.circle(mm(x), mm(Y(yt)), mm(r), gp = gpar(fill = fill, col = col, lwd = lwd))
pill <- function(s, x, yt, fill = INK, col = "white", size = 7.5, w = NULL) {
  w <- if (is.null(w)) nchar(s) * 1.55 + 3 else w
  grid.roundrect(mm(x), mm(Y(yt)), mm(w), mm(4.2), just = c("left", "top"), r = mm(2.1), gp = gpar(fill = fill, col = NA))
  txt(s, x + w / 2, yt + 2.1, size = size, face = "bold", col = col, hjust = 0.5, vjust = 0.5); w
}
sig_card <- function(x, yt, col, w = 6.2, h = 7.2) {
  box(x, yt, w, h, fill = "white", col = col, lwd = 1.1, r = 1)
  for (i in 1:4) grid.segments(mm(x + 1.2), mm(Y(yt + 1.4 + i * 1.15)), mm(x + w - 1.2 - (i %% 2) * 1.2), mm(Y(yt + 1.4 + i * 1.15)), gp = gpar(col = col, lwd = 1.1))
}

draw_fig1 <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  txt("WHAT IS EVALUATED", 2, 2.2, size = 7.5, face = "bold", col = MUTED)
  box(2, 6.5, 52, 35, fill = TINT, col = NA)
  sig_card(5, 10, C_PUB)
  txt("36 published prognostic\nsignatures", 14, 10.1, size = 8, face = "bold")
  txt("2022–2026 · 2–37 genes\n32 trained on GSE65682", 14, 17.9, size = 7.5, col = INK2)
  for (k in 0:2) box(5 + k * 0.9, 27.3 + k * 0.9, 6.2, 7.2, fill = "white", col = INK2, lwd = 0.9, r = 1)
  txt("6 endotyping and\nscoring tools", 14, 27.4, size = 8, face = "bold")
  txt("SRS · CTS · SUBSPACE ·\nMars1 · Sweeney (×2)", 14, 35.0, size = 7.5, col = INK2)

  txt("ADULT SEPSIS COHORTS · BLOOD TRANSCRIPTOMES", 58, 2.2, size = 7.5, face = "bold", col = MUTED)
  chips <- list(
    list(role = "Training", ds = "GSE65682", where = "Netherlands · ICU", r = 1, c = 1),
    list(role = "Main evaluation", ds = "GSE185263", where = "ED/ICU · 4 countries", r = 1, c = 2),
    list(role = "", ds = "GSE63042", where = "USA · ED", r = 1, c = 3),
    list(role = "", ds = "GSE236713", where = "UK · ICU", r = 1, c = 4),
    list(role = "Replication", ds = "GSE272769", where = "USA · ICU", r = 2, c = 1),
    list(role = "Sensitivity", ds = "GSE95233", where = "France · ICU", r = 2, c = 2),
    list(role = "", ds = "GSE54514", where = "Australia · ICU", r = 2, c = 3),
    list(role = "Local validation†", ds = "CMAISE", where = "China · 32 centres", r = 2, c = 4))
  cw <- 26.4; gap <- 1.2; x0 <- 58
  for (ch in chips) {
    x <- x0 + (ch$c - 1) * (cw + gap); yt <- if (ch$r == 1) 10 + 0.5 else 27.5 + 0.5
    if (nzchar(ch$role)) txt(ch$role, x + 0.3, yt - 3.6, size = 7, face = "bold", col = if (ch$ds == "CMAISE") C_RND else INK2)
    local <- ch$ds == "CMAISE"
    box(x, yt, cw, 12.4, fill = if (local) "#eef4fc" else "white", col = if (local) C_RND else LINE, lwd = if (local) 0.9 else 0.6, lty = if (local) 2 else 1)
    txt(ch$ds, x + 1.6, yt + 1.5, size = 8, face = "bold")
    txt(ch$where, x + 1.6, yt + 5.2, size = 7, col = INK2)
    txt(nd(ch$ds), x + 1.6, yt + 8.8, size = 7, col = INK2)
  }

  ty <- 49
  txt("ORDER OF WORK", 2, 44.6, size = 7.5, face = "bold", col = MUTED)
  steps <- c("Protocol, code and\nparameters fixed",
             "7,200 random signatures\nbuilt and archived",
             "Signatures evaluated in\nevaluation cohorts",
             "CMAISE addendum\nfixed, then analysed")
  sx <- c(2, 44.5, 87, 129.5)
  grid.segments(mm(sx[1] + 3), mm(Y(ty + 2.6)), mm(sx[4] + 3), mm(Y(ty + 2.6)), gp = gpar(col = LINE, lwd = 1.4))
  for (i in seq_along(steps)) {
    grid.rect(mm(sx[i] + 5.6), mm(Y(ty - 0.6)), mm(33), mm(6.8), just = c("left", "top"), gp = gpar(fill = "white", col = NA))
    grid.circle(mm(sx[i] + 2.6), mm(Y(ty + 2.6)), mm(2.6), gp = gpar(fill = if (i == 4) C_RND else INK, col = "white", lwd = 1.2))
    txt(as.character(i), sx[i] + 2.6, ty + 2.6, size = 7.5, face = "bold", col = "white", hjust = 0.5, vjust = 0.5)
    txt(steps[i], sx[i] + 6.6, ty - 0.1, size = 7.5, col = INK)
  }

  cy <- 59; ch <- 73; cwid <- 53.3; cx <- c(2, 58.35, 114.7); hh <- 19
  hdr <- list(
    list(p = "PART I", q = "Are published signatures better\nthan random genes put through\na common modelling pipeline?"),
    list(p = "PART II", q = "How much of a signature\nscore is leukocyte\ncomposition?"),
    list(p = "PART III", q = "Does a patient's label\ndepend on who else is\nin the batch?"))
  meth <- list(
    "• 200 random-theme signatures per\n   published signature (7,200)\n• DEG → Cox → LASSO → Cox, with\n   matched theme size and training n\n• pooled AUROC in held-out cohorts;\n   replication: GSE272769, CMAISE",
    "• LM22 deconvolution → 5 groups\n• accuracy vs blood counts (GSE54514)\n• adjusted R² of score on composition\n   (no outcomes used)\n• AUROC after adjustment\n   (descriptive)",
    "• same patients, different partners:\n   + healthy controls, ICU + ED,\n   composition titration\n• tools grouped by reference:\n   sample, batch, external set\n• no outcomes used")
  hyp <- list(
    c("H1 (primary)", "published − random: upper 95% CI < 0.05", "H2", "random: lower 95% CI > 0.5"),
    c("H3", "median adjusted R²:\nlower 95% CI > 0.5", "", ""),
    c("H4", "SRS re-labelling when controls\njoin: lower 95% CI > 5%", "", ""))
  nextq <- c("Where does the prognostic\nsignal come from?", "Is the score also shaped by\nthe other samples in a batch?")
  for (i in 1:3) {
    x <- cx[i]
    box(x, cy, cwid, ch, fill = "white", col = LINE, lwd = 0.7)
    grid.roundrect(mm(x), mm(Y(cy)), mm(cwid), mm(hh), just = c("left", "top"), r = mm(1.8), gp = gpar(fill = TINT, col = NA))
    grid.rect(mm(x), mm(Y(cy + hh - 1.8)), mm(cwid), mm(1.8), just = c("left", "top"), gp = gpar(fill = TINT, col = NA))
    pill(hdr[[i]]$p, x + 2, cy + 2, w = 15)
    txt(hdr[[i]]$q, x + 2, cy + 7.2, size = 8, face = "bold", lh = 1.05)
    iy <- cy + hh + 1.5
    if (i == 1) {
      sig_card(x + 6, iy + 2.2, C_PUB)
      txt("Published", x + 9.1, iy + 11.8, size = 7, col = INK2, hjust = 0.5)
      txt("vs", x + 20.5, iy + 5.8, size = 8, face = "bold", col = INK2, hjust = 0.5)
      for (k in 2:0) sig_card(x + 27 + k * 1.3, iy + 2.2 + k * 1.1, C_RND)
      txt("×200", x + 38.6, iy + 4.4, size = 8, face = "bold", col = C_RND)
      txt("Matched random", x + 32.5, iy + 11.8, size = 7, col = INK2, hjust = 0.5)
    } else if (i == 2) {
      segw <- c(0.46, 0.14, 0.24, 0.08, 0.08) * 30; xs <- x + 3
      for (k in 1:5) { grid.rect(mm(xs), mm(Y(iy + 2.4)), mm(segw[k]), mm(5.5), just = c("left", "top"), gp = gpar(fill = RAMP[k], col = "white", lwd = 0.8)); xs <- xs + segw[k] }
      labs <- c("N", "M", "T", "", ""); xs <- x + 3
      for (k in 1:5) { if (nzchar(labs[k])) txt(labs[k], xs + segw[k] / 2, iy + 5.15, size = 7, face = "bold", col = if (k >= 3) "white" else INK, hjust = 0.5, vjust = 0.5); xs <- xs + segw[k] }
      arr(x + 34.5, iy + 5.15, x + 39.5, iy + 5.15)
      grid.roundrect(mm(x + 40.3), mm(Y(iy + 1.6)), mm(10.5), mm(7), just = c("left", "top"), r = mm(1), gp = gpar(fill = "white", col = C_PUB, lwd = 1.1))
      txt("score", x + 45.55, iy + 5.1, size = 7, face = "bold", col = INK, hjust = 0.5, vjust = 0.5)
      txt("N neutrophils · M monocytes/DC\nT T cells · NK cells · B/plasma cells", x + 3, iy + 9.8, size = 7, col = INK2)
    } else {
      bx <- c(x + 3, x + 28.3); bw <- 22
      for (b in 1:2) {
        box(bx[b], iy + 1, bw, 11.6, fill = "white", col = INK2, lwd = 0.6, r = 1.2)
        txt(if (b == 1) "patients only" else "+ healthy controls", bx[b] + bw / 2, iy + 13.6, size = 7, col = INK2, hjust = 0.5)
        pts <- cbind(c(4, 8, 12, 6, 10, 14), c(4.6, 4, 5.2, 8.6, 9.4, 8.8))
        for (k in 1:6) dot(bx[b] + pts[k, 1], iy + pts[k, 2], r = 0.95, fill = INK2)
        dot(bx[b] + 8, iy + 4, r = 1.75, fill = NA, col = INK, lwd = 1.2)
        if (b == 2) { hc <- cbind(c(17.5, 19.2, 17.5, 19.4), c(3.6, 6.2, 8.8, 10.9)); for (k in 1:4) dot(bx[b] + hc[k, 1], iy + hc[k, 2], r = 0.95, fill = C_AQ) }
      }
      arr(x + 25.3, iy + 6.8, x + 27.9, iy + 6.8)
    }
    txt(meth[[i]], x + 2.4, cy + hh + 18.8, size = 7.5, col = INK, lh = 1.1)
    hb <- cy + ch - 14.2
    box(x + 2, hb, cwid - 4, 12.4, fill = TINT, col = NA, r = 1.2)
    h <- hyp[[i]]
    txt(h[1], x + 3.6, hb + 1.6, size = 7.5, face = "bold")
    txt(h[2], x + 3.6, hb + 4.9, size = 7, col = INK2, lh = 1.05)
    if (nzchar(h[3])) { txt(h[3], x + 3.6, hb + 8.4, size = 7.5, face = "bold"); txt(h[4], x + 8.6, hb + 8.55, size = 7, col = INK2) }
    if (i < 3) {
      gx <- x + cwid + 1.5
      arr(x + cwid + 0.2, cy + 33, x + cwid + 2.9, cy + 33, col = INK, lwd = 1.1)
      txt(nextq[i], gx, cy + ch + 1.6, size = 7, face = "italic", col = INK2, hjust = 0.5)
    }
  }
  fy <- cy + ch + 9.2
  grid.segments(mm(2), mm(Y(fy - 1)), mm(168), mm(Y(fy - 1)), gp = gpar(col = LINE, lwd = 0.6))
  txt("H1 is the primary hypothesis (non-superiority margin 0.05). H2–H4 form one family, Holm-adjusted (one-sided α = 0.025).", 2, fy, size = 7, col = INK2)
  txt("† CMAISE was added after the main results; its criteria (L1, L2) were fixed before any signature–outcome association was computed.", 2, fy + 3.4, size = 7, col = INK2)
}

dev_size <- c(W, H) / 25.4
cairo_pdf(file.path(OUTD, "Fig1_study_design.pdf"), width = dev_size[1], height = dev_size[2], family = FONT); draw_fig1(); invisible(dev.off())
cairo_ps(file.path(OUTD, "Fig1_study_design.eps"), width = dev_size[1], height = dev_size[2], family = FONT, fallback_resolution = 600); draw_fig1(); invisible(dev.off())
png(file.path(OUTD, "Fig1_study_design.png"), width = dev_size[1], height = dev_size[2], units = "in", res = 600, type = "cairo"); draw_fig1(); invisible(dev.off())
tiff(file.path(OUTD, "Fig1_study_design.tiff"), width = dev_size[1], height = dev_size[2], units = "in", res = 600, type = "cairo", compression = "lzw"); draw_fig1(); invisible(dev.off())
message("Fig 1 saved to ", OUTD)
