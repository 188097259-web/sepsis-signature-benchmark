source("scripts/_common.R")
ALPHA <- 0.025
od <- outdir()

b2 <- readRDS(file.path(OUT, "part1", "bootstrap_main.rds"))
e1_bt <- b2$bt["E1", "all", ]; e1_bt <- e1_bt[!is.na(e1_bt)]
h12 <- fread(file.path(OUT, "part1", "H1_H2_result.csv"))
p2 <- (1 + sum(e1_bt <= 0.5)) / (length(e1_bt) + 1)

h3 <- fread(file.path(OUT, "part2", "H3_result.csv"))
h4 <- fread(file.path(OUT, "batch", "H4_result.csv"))

fam <- data.table(hypothesis = c("H2", "H3", "H4"),
                  rule = c("E1 lower 95% CI > 0.50", "median adjusted R2 basic-bootstrap lower 95% CI > 0.5", "SRS label change C0 vs C1 lower 95% CI > 5%"),
                  estimate = c(h12[hypothesis == "H2"]$estimate, h3$median_adj_r2, h4$pooled),
                  lo = c(h12[hypothesis == "H2"]$lo, h3$lo, h4$lo), hi = c(h12[hypothesis == "H2"]$hi, h3$hi, h4$hi),
                  ci_rule_met = c(h12[hypothesis == "H2"]$supported, h3$supported, h4$supported),
                  p_one_sided = c(p2, h3$p_one_sided, h4$p_one_sided))
fam[, p_holm := p.adjust(p_one_sided, "holm")]
fam[, supported_after_holm := p_holm < ALPHA]
fam[, note := fifelse(ci_rule_met & !supported_after_holm, "CI rule met but not after Holm", "")]
h1 <- h12[hypothesis == "H1", .(hypothesis, rule, estimate, lo, hi, ci_rule_met = supported, p_one_sided = NA_real_, p_holm = NA_real_, supported_after_holm = NA, note = "primary; not in Holm family")]
out <- rbind(h1, fam)
fwrite(out, file.path(od, "hypotheses_final.csv")); print(out)
message("done")
