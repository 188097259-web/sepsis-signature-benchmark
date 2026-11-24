source("scripts/_common.R")
suppressPackageStartupMessages({ library(survival); library(glmnet) })
set.seed(SEED)
od <- outdir("part1")
P <- fromJSON(file.path(PROJ, "params/pipeline_params.json"))
N_R2 <- P$n_runs_per_model; N_R1 <- P$n_runs_fixed_size; N_R3 <- P$n_runs_full_universe
if (DRY_RUN) { N_R2 <- 3; N_R1 <- 5; N_R3 <- 2; P$n_runs_nodeg_per_model <- 2 }

tr_all <- load_cohort(TRAIN, sepsis_only = FALSE)
U <- Reduce(intersect, c(list(rownames(tr_all$expr)), lapply(c(MAIN, REPL), function(g) rownames(load_cohort(g)$expr))))
message("gene universe U: ", length(U))
ex <- tr_all$expr[U, ]; ph <- tr_all$ph
sep <- which(ph$group == "sepsis" & !is.na(ph$event28) & !is.na(ph$time28)); ctl <- which(ph$group == "control")
ev <- ph$event28[sep]; tm <- ph$time28[sep]
tm[tm <= 0] <- 0.5
if (DRY_RUN) { gsim <- sample(U, 5); lp <- as.numeric(scale(colSums(ex[gsim, sep]))); h <- 0.01 * exp(0.8 * lp)
  tt0 <- rexp(length(sep), h); ev <- as.integer(tt0 <= 28); tm <- pmin(tt0, 28) }
Xs <- t(ex[, sep]); message("training: ", nrow(Xs), " patients, ", sum(ev), " deaths")

grp <- factor(c(rep("sepsis", length(sep)), rep("control", length(ctl))), levels = c("control", "sepsis"))
fit <- eBayes(lmFit(ex[, c(sep, ctl)], model.matrix(~ grp)))
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
deg <- rownames(tt)[abs(tt$logFC) >= P$deg_abs_log2fc & tt[[if (P$deg_p_type == "adj") "adj.P.Val" else "P.Value"]] < P$deg_p]
message("DEGs: ", length(deg))

unicox_p <- function(X, t, e) apply(X, 2, function(x) summary(coxph(Surv(t, e) ~ x))$coefficients[1, 5])

pipeline_once <- function(S, n_train, seed, settings = P) {
  set.seed(seed)
  theme <- sample(U, min(S, length(U)))
  tr <- if (n_train >= length(ev)) seq_along(ev) else
    c(sample(which(ev == 1), round(n_train * mean(ev))), sample(which(ev == 0), n_train - round(n_train * mean(ev))))
  cand <- if (isTRUE(settings$use_deg)) intersect(theme, deg) else theme
  if (length(cand) < 2) return(list(fail = "deg"))
  pu <- unicox_p(Xs[tr, cand, drop = FALSE], tm[tr], ev[tr]); cand <- cand[pu < settings$unicox_p]
  if (length(cand) < 2) return(list(fail = "unicox"))
  cv <- tryCatch(cv.glmnet(Xs[tr, cand, drop = FALSE], Surv(tm[tr], ev[tr]), family = "cox", nfolds = 10), error = function(e) NULL)
  if (is.null(cv)) return(list(fail = "lasso_error"))
  b <- as.numeric(coef(cv, s = settings$lambda)); sel <- cand[b != 0]
  if (length(sel) < 2) return(list(fail = "lasso_lt2"))
  lasso_coef <- setNames(b[b != 0], sel)
  dfm <- data.frame(t = tm[tr], e = ev[tr], Xs[tr, sel, drop = FALSE], check.names = FALSE); names(dfm)[-(1:2)] <- paste0("g", seq_along(sel))
  f <- tryCatch(coxph(Surv(t, e) ~ ., data = dfm), error = function(e) NULL, warning = function(w) NULL)
  if (is.null(f) || anyNA(coef(f))) return(list(fail = "multicox"))
  fs <- tryCatch(step(f, direction = "both", trace = 0), error = function(e) NULL)
  sg <- if (!is.null(fs) && length(coef(fs)) >= 2) sel[as.integer(sub("^g", "", names(coef(fs))))] else character(0)
  list(genes = sel, coef = unname(coef(f)), lasso_genes = sel, lasso_coef = unname(lasso_coef),
       step_genes = sg, step_coef = if (length(sg)) unname(coef(fs)) else numeric(0), n_cand = length(cand), n_train = length(tr))
}

mp <- as.data.table(P$models)
NCORES <- as.integer(Sys.getenv("NCORES", "10"))
run_model <- function(i, n_target, settings, base, prefix) {
  got <- 0; tries <- 0; fl <- character(0); out <- list()
  while (got < n_target && tries < n_target * P$max_tries_factor) {
    tries <- tries + 1; sd <- base + i * 1e6 + tries
    r <- pipeline_once(mp$theme_size[i], mp$n_train[i], sd, settings)
    if (!is.null(r$fail)) { fl <- c(fl, r$fail); next }
    got <- got + 1; out[[got]] <- c(list(key = sprintf("%s_%s_%03d", prefix, mp$model_id[i], got), matched = mp$model_id[i], seed = sd), r)
  }
  list(models = out, fails = data.table(model_id = mp$model_id[i], theme_size = mp$theme_size[i], n_train = mp$n_train[i], success = got, tries = tries,
       fail_deg = sum(fl == "deg"), fail_unicox = sum(fl == "unicox"), fail_lasso = sum(fl %in% c("lasso_error", "lasso_lt2")), fail_multicox = sum(fl == "multicox")))
}
res2 <- parallel::mclapply(seq_len(nrow(mp)), run_model, n_target = N_R2, settings = P, base = SEED * 10, prefix = "R2", mc.cores = NCORES, mc.preschedule = FALSE)
R2 <- unlist(lapply(res2, `[[`, "models"), recursive = FALSE); fails <- rbindlist(lapply(res2, `[[`, "fails"))
fwrite(fails, file.path(od, "R2_pipeline_failures.csv")); print(fails[, .(model_id, success, tries)])

nodeg <- modifyList(P, list(use_deg = FALSE))
res2b <- parallel::mclapply(seq_len(nrow(mp)), run_model, n_target = P$n_runs_nodeg_per_model, settings = nodeg, base = SEED * 20, prefix = "R2b", mc.cores = NCORES, mc.preschedule = FALSE)
R2b <- unlist(lapply(res2b, `[[`, "models"), recursive = FALSE)
fwrite(rbindlist(lapply(res2b, `[[`, "fails")), file.path(od, "R2b_pipeline_failures.csv"))
seed <- SEED * 30

M <- load_models(); fixed_sizes <- sort(unique(sapply(M$mt$model_id, function(id) nrow(M$gn[model_id == id][!duplicated(gene)]))))
message("R1 sizes: ", paste(fixed_sizes, collapse = ","))
R1 <- list()
for (k in fixed_sizes) for (j in seq_len(N_R1)) {
  seed <- seed + 1; set.seed(seed); gs <- sample(U, k)
  dfm <- data.frame(t = tm, e = ev, Xs[, gs, drop = FALSE], check.names = FALSE); names(dfm)[-(1:2)] <- paste0("g", seq_len(k))
  f <- tryCatch(coxph(Surv(t, e) ~ ., data = dfm), error = function(e) NULL, warning = function(w) NULL)
  if (is.null(f) || anyNA(coef(f))) next
  R1[[length(R1) + 1]] <- list(key = sprintf("R1_k%02d_%04d", k, j), k = k, seed = seed, genes = gs, coef = unname(coef(f)))
}
R3 <- list()
for (j in seq_len(N_R3)) { seed <- seed + 1; r <- pipeline_once(length(U), length(ev), seed)
  if (is.null(r$fail)) R3[[length(R3) + 1]] <- c(list(key = sprintf("R3_%03d", j), seed = seed), r) }

store <- list(U = U, deg = deg, R1 = R1, R2 = R2, R2b = R2b, R3 = R3, params = P, created = as.character(Sys.time()), dry_run = DRY_RUN)
f <- file.path(od, "random_models.rds"); saveRDS(store, f)
writeLines(c(paste("random_models.rds md5:", md5_of(f)), paste("R1:", length(R1), "R2:", length(R2), "R3:", length(R3))), file.path(od, "random_models_checksum.txt"))
writeLines(capture.output(sessionInfo()), file.path(od, "sessionInfo_01.txt"))
message("done: R1 ", length(R1), ", R2 ", length(R2), ", R3 ", length(R3))
