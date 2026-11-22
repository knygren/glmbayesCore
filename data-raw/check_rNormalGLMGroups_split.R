# Check the mode / envelope / simulate split inside rNormalGLMGroups.
#
# rNormalGLMGroups now runs, for each group, the current rNormalGLM step and
# then the three split functions on the same inputs. It stops if the mode or
# any envelope field differs, so a call that returns has passed the exact
# check. It also returns:
#   coefficients_split  the split path's draw per group (not used downstream)
#   timing              seconds in rNormalGLM, and in the split mode,
#                       envelope (setup / EnvelopeEval / finalize) and
#                       simulate phases
#
# This script:
#   1. calls rNormalGLMGroups repeatedly (each call is an exact check);
#   2. compares coefficients and coefficients_split in distribution;
#   3. reports the timing.
# Poisson (47 groups x 4 coefficients) is the main case.
#
# Run after reinstalling glmbayesCore:
#   Rscript data-raw/check_rNormalGLMGroups_split.R            # R = 2000
#   Rscript data-raw/check_rNormalGLMGroups_split.R 500        # fewer reps
#   Rscript data-raw/check_rNormalGLMGroups_split.R 2000 all   # + binomial, Gaussian

library(glmbayesCore)

args <- commandArgs(trailingOnly = TRUE)
r_arg <- suppressWarnings(as.integer(args[1L]))
R_reps <- if (length(args) >= 1L && !is.na(r_arg)) r_arg else 2000L
run_all <- "all" %in% args
alpha_fw <- 0.01

groups_fn <- glmbayesCore:::.rNormalGLMGroups_cpp

# 47 groups, intercept + 3 binary predictors, 3-50 rows per group.
make_case <- function(family, k = 47L, l1 = 4L) {
  n_g <- sample(3:50, k, replace = TRUE)
  grp <- rep(seq_len(k), n_g)
  l2 <- length(grp)
  x <- cbind(1, matrix(rbinom(l2 * (l1 - 1L), 1, 0.4), l2, l1 - 1L))
  colnames(x) <- c("(Intercept)", paste0("x", seq_len(l1 - 1L)))
  b_true <- matrix(rnorm(k * l1, sd = 0.5), k, l1)
  b_true[, 1] <- b_true[, 1] - 1
  eta <- rowSums(x * b_true[grp, , drop = FALSE])
  y <- switch(family$family,
    binomial = rbinom(l2, 1, plogis(eta)),
    poisson  = rpois(l2, exp(eta + 1)),
    gaussian = eta + rnorm(l2, sd = 1.2)
  )
  mu <- matrix(rnorm(l1 * k, sd = 0.3), l1, k)
  mu[1, ] <- mu[1, ] - 1
  P_blocks <- replicate(k, diag(1 / c(0.8, 0.6, 0.4, 0.7)[seq_len(l1)]^2, l1),
                        simplify = FALSE)
  ff <- glmbfamfunc(family)
  list(
    name = paste(family$family, family$link), family = family$family,
    link = family$link, y = y, x = x, offset = rep(0, l2), wt = rep(1, l2),
    dispersion = if (family$family == "gaussian") 1.44 else 1,
    mu = mu, P_blocks = P_blocks, row_blocks = split(seq_len(l2), grp),
    f2 = ff$f2, f3 = ff$f3, k = k, l1 = l1
  )
}

call_case <- function(cs) {
  groups_fn(1L, cs$y, cs$x, cs$offset, cs$wt, cs$dispersion, cs$mu,
            cs$P_blocks, TRUE, cs$row_blocks, cs$f2, cs$f3,
            family = cs$family, link = cs$link)
}

run_case <- function(cs, R) {
  orig <- split <- array(NA_real_, c(R, cs$k, cs$l1))
  tm <- c(rNormalGLM = 0, mode = 0, envelope = 0, env_setup = 0,
          env_eval = 0, env_finalize = 0, simulate = 0)
  for (r in seq_len(R)) {
    out <- call_case(cs)
    orig[r, , ] <- out$coefficients
    split[r, , ] <- out$coefficients_split
    tm <- tm + out$timing
  }
  cat(sprintf("[%s] exact: mode and envelope identical in all %d calls x %d groups\n",
              cs$name, R, cs$k))

  n_tests <- cs$k * cs$l1
  thr <- alpha_fw / n_tests
  p_mean <- p_sd <- p_ks <- matrix(NA_real_, cs$k, cs$l1)
  for (b in seq_len(cs$k)) for (j in seq_len(cs$l1)) {
    a <- orig[, b, j]
    d <- split[, b, j]
    z <- (mean(a) - mean(d)) / sqrt(var(a) / R + var(d) / R)
    p_mean[b, j] <- 2 * pnorm(-abs(z))
    p_sd[b, j] <- stats::var.test(a, d)$p.value
    p_ks[b, j] <- suppressWarnings(stats::ks.test(a, d)$p.value)
  }
  ok <- min(p_mean) > thr && min(p_sd) > thr && min(p_ks) > thr
  cat(sprintf("    draws (R = %d, %d tests each, Bonferroni threshold %.2g): min p mean %.3g, SD %.3g, KS %.3g => %s\n",
              R, n_tests, thr, min(p_mean), min(p_sd), min(p_ks),
              if (ok) "consistent" else "INCONSISTENT"))

  split_tot <- sum(tm[c("mode", "envelope", "simulate")])
  cat(sprintf("    timing per group call: rNormalGLM %.3f ms; split %.3f ms = mode %.0f%%, envelope %.0f%%, simulate %.0f%%\n",
              1000 * tm[["rNormalGLM"]] / (R * cs$k),
              1000 * split_tot / (R * cs$k),
              100 * tm[["mode"]] / split_tot,
              100 * tm[["envelope"]] / split_tot,
              100 * tm[["simulate"]] / split_tot))
  cat(sprintf("    envelope parts: setup %.0f%%, EnvelopeEval %.0f%%, finalize %.0f%% (of envelope time)\n",
              100 * tm[["env_setup"]] / tm[["envelope"]],
              100 * tm[["env_eval"]] / tm[["envelope"]],
              100 * tm[["env_finalize"]] / tm[["envelope"]]))
  ok
}

cases <- list(make_case(poisson()))
if (run_all) {
  cases <- c(cases, list(make_case(binomial(link = "logit")),
                         make_case(gaussian())))
}
reps <- c(R_reps, rep(max(200L, R_reps %/% 4L), length(cases) - 1L))
ok <- mapply(run_case, cases, reps)

cat(sprintf("\nSUMMARY: exact PASS (no stop); distribution %s\n",
            if (all(ok)) "PASS" else "FAIL"))
