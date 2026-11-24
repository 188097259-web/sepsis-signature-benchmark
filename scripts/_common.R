suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(limma); library(org.Hs.eg.db) })

# Run every script from the repository root, or set SEPSIS_BENCH_ROOT. The LM22 matrix and the CMAISE data are not
# distributed with this repository: point LM22_FILE and CMAISE_DIR to local copies (see README.md).
PROJ <- Sys.getenv("SEPSIS_BENCH_ROOT", ".")
DATA <- file.path(PROJ, "data")
LM22_FILE  <- Sys.getenv("LM22_FILE", file.path(DATA, "external", "LM22.txt"))
CMAISE_DIR <- Sys.getenv("CMAISE_DIR", file.path(DATA, "cmaise"))
LM22_MD5  <- "f315a17039596ef73f83fbb7466677dd"
SEED <- 20261116
MAIN <- c("GSE185263", "GSE63042", "GSE236713")
REPL <- "GSE272769"
SENS <- c("GSE95233", "GSE54514")
TRAIN <- "GSE65682"
NNLS_INVALID <- "GSE236713"
PART2_MAIN <- setdiff(MAIN, NNLS_INVALID)
REF_PMIDS <- c("35027333","41028542","28864056","36322631","41028543","29449546","26917434","29537985","35847809","33181462")
Z_MODELS <- c("Jiang2022_IRRG9", "Chen2023_IRG7", "Feng2024_ROS37", "Li2026_GMGscore6")
FLAGGED <- c("Zhang2024_cuproptosis3","Zheng2022_IRG22","Li2025_immune21","Luo2026_RCDscore6","Luo2025_MRMG4",
             "Zhou2022_7geneLogistic","Xue2026_URS3","Jiang2024_PRG10","Sun2026_NRG5")

DRY_RUN <- identical(Sys.getenv("DRY_RUN"), "1")
OUT <- file.path(PROJ, if (DRY_RUN) "results_dryrun" else "results")
outdir <- function(...) { p <- file.path(OUT, ...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }

official <- function(sym) { o <- suppressWarnings(alias2SymbolTable(sym, species = "Hs")); ifelse(is.na(o), sym, o) }
scale_mat <- function(m, how) { mu <- rowMeans(m); if (how == "center") m - mu else (m - mu) / pmax(apply(m, 1, sd), 1e-8) }

load_cohort <- function(g, sepsis_only = TRUE) {
  if (g == REPL) {
    e <- readRDS(file.path(DATA, "cohorts/GSE272769_expr.rds"))
    st <- fread(file.path(DATA, "cohorts/GSE272769_sample_status.csv"))[dataset == g]
    st <- st[match(colnames(e), sample)]
    ph <- data.table(sample = st$sample, group = st$group, death = as.integer(st$outcome_death), setting = st$setting)
  } else {
    x <- readRDS(file.path(DATA, "cohorts", paste0(g, ".rds"))); e <- x$expr; ph <- as.data.table(x$pheno)
  }
  if (g == TRAIN) ph <- merge_gse65682_time(ph)
  keep <- if (sepsis_only) ph$group == "sepsis" else rep(TRUE, nrow(ph))
  out <- list(expr = e[, keep, drop = FALSE], ph = ph[keep])
  if (DRY_RUN && sepsis_only) out$ph$death <- sample(out$ph$death)
  out
}
merge_gse65682_time <- function(ph) {
  f <- file.path(DATA, "cohorts/GSE65682_series_matrix_header.txt.gz")
  L <- readLines(gzfile(f)); get <- function(key) { l <- grep(paste0("^!Sample_", key), L, value = TRUE); strsplit(l, "\t")[[1]][-1] }
  ids <- gsub('"', "", get("geo_accession"))
  ch <- grep("^!Sample_characteristics_ch1", L, value = TRUE)
  pick <- function(field) { l <- ch[grepl(field, ch)][1]; v <- gsub('"', "", strsplit(l, "\t")[[1]][-1]); sub(".*: ", "", v) }
  tt <- data.table(sample = ids, event28 = suppressWarnings(as.integer(pick("mortality_event_28days"))),
                   time28 = suppressWarnings(as.numeric(pick("time_to_event_28days"))))
  tt <- tt[match(ph$sample, tt$sample)]
  ph[, `:=`(event28 = tt$event28, time28 = tt$time28)]
  ph
}

load_models <- function() {
  mt <- fread(file.path(DATA, "signatures/model_table.csv"), colClasses = "character")
  mt <- mt[!(pmid %in% REF_PMIDS) & tier == "A_coef_complete"]
  gn <- fread(file.path(DATA, "signatures/genes.csv"), colClasses = "character")
  gn[, coef := suppressWarnings(as.numeric(coef))]; gn[, gene := official(gene_symbol)]
  higher <- sapply(mt$model_id, function(id) {
    r <- fromJSON(file.path(DATA, "signatures/extraction_records", paste0(mt[model_id == id]$pmid[1], ".json")), simplifyVector = FALSE)
    m <- Filter(function(z) identical(z$model_id, id), r$models)[[1]]
    h <- if (is.null(m$higher_score_means)) "higher risk" else m$higher_score_means
    if (grepl("lower", h)) -1 else 1 })
  mt[, GSE272769 := "eligible"]
  list(mt = mt, gn = gn, higher = higher)
}
main_scale <- function(id) if (id %in% Z_MODELS) "z" else "center"
eligible <- function(M, id, g) M$mt[model_id == id][[g]] == "eligible"
lin_score <- function(expr, genes, w, how = "center", thr = 0.8) {
  present <- genes %in% rownames(expr); if (mean(present) < thr) return(NULL)
  X <- scale_mat(expr[genes[present], , drop = FALSE], how)
  as.numeric(crossprod(w[present], X))
}
score_published <- function(M, id, expr) {
  g <- M$gn[model_id == id][!duplicated(gene)]
  s <- lin_score(expr, g$gene, g$coef, main_scale(id)); if (is.null(s)) NULL else M$higher[[id]] * s
}

delong <- function(s, y) {
  pos <- s[y == 1]; neg <- s[y == 0]; n1 <- length(pos); n0 <- length(neg)
  ra <- rank(c(pos, neg)); v10 <- (ra[seq_len(n1)] - rank(pos)) / n0; v01 <- 1 - (ra[n1 + seq_len(n0)] - rank(neg)) / n1
  c(auc = mean(v10), se = sqrt(var(v10) / n1 + var(v01) / n0))
}
pool_fixed <- function(a, se) { l <- qlogis(a); w <- 1 / (se / (a * (1 - a)))^2; est <- sum(w * l) / sum(w); sep <- sqrt(1 / sum(w))
  c(auc = plogis(est), lo = plogis(est - 1.96 * sep), hi = plogis(est + 1.96 * sep)) }
strat_idx <- function(Y) lapply(Y, function(y) c(sample(which(y == 1), replace = TRUE), sample(which(y == 0), replace = TRUE)))
full_idx <- function(Y) lapply(Y, seq_along)
md5_of <- function(f) unname(tools::md5sum(f))
