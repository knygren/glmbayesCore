#' Summarizing Bayesian mixed model distribution functions
#'
#' These functions are \code{\link{methods}} for class \code{rlmerb} or
#' \code{summary.rlmerb} objects, following the layout of
#' \code{\link[glmbayesCore]{summary.rglmb}} extended to the two-block mixed
#' model: Block~2 (population/level-2) tables are reported per random-effects
#' component, and Block~1 (group-level) effects are summarized in aggregate.
#'
#' @details
#' Like \code{\link[glmbayesCore]{summary.rglmb}}, which fits its own
#' \code{\link[stats]{glm}} maximum-likelihood reference, this method computes
#' its reference comparison internally from the \code{\link{model_setup}}
#' object carried on the fit --- nothing needs to be supplied by the caller.
#'
#' The reference is the \code{lmer} fit built by \code{model_setup()}, whose
#' fixed effects are on a different footing from the Block~2 hyperparameters:
#' \eqn{\gamma} is indexed by (random-effect component, level-2 parameter),
#' whereas \code{lmer} reports one flat coefficient vector. The two are put on
#' a common footing before being tabulated (see
#' \sQuote{Reference correspondence} below), so the \code{lmer} column is the
#' coefficient that estimates the same quantity, not merely one with a similar
#' name.
#'
#' @section Reference correspondence:
#' Writing \eqn{k} for the random-effect component and \eqn{p} for the level-2
#' parameter within that component, \eqn{\gamma_{k,p}} maps to an \code{lmer}
#' fixed effect as follows:
#' \describe{
#'   \item{\eqn{k} and \eqn{p} both the intercept}{the \code{lmer} intercept.}
#'   \item{\eqn{p} the intercept, \eqn{k} a slope}{the main effect \eqn{k}: the
#'     average of the random slope \eqn{k} across groups.}
#'   \item{\eqn{k} the intercept, \eqn{p} a covariate}{the main effect \eqn{p}:
#'     a group-level predictor of the random intercept.}
#'   \item{both \eqn{k} and \eqn{p} slopes}{the cross-level interaction
#'     \eqn{p{:}k}, since \eqn{p} moderates the random slope \eqn{k}. The
#'     interaction column recorded in
#'     \code{model_setup()$popef.moderation} is used when available, and the
#'     \eqn{p{:}k} / \eqn{k{:}p} orderings are tried otherwise.}
#' }
#' Coefficients with no counterpart are reported as \code{NA} rather than
#' matched approximately.
#'
#' @param object An object of class \code{"rlmerb"}. An \code{rglmerb} fit
#'   from \pkg{lmebayesCore} is accepted as well.
#' @param groups Optional character vector of grouping levels for which to
#'   include a per-group Block~1 detail table. When \code{NULL} (default) only
#'   the aggregate \code{groupef_overview} is returned.
#' @param x An object of class \code{"summary.rlmerb"}.
#' @param digits The number of significant digits to use when printing.
#' @param \ldots Additional optional arguments.
#'
#' @return \code{summary.rlmerb} returns an object of class
#'   \code{"summary.rlmerb"}, a list with components:
#'   \item{call, formula}{The originating call and model formula.}
#'   \item{n}{Integer number of stored draws.}
#'   \item{sim_method_used, method, draw_engine}{Character labels for the
#'     engine that produced the draws.}
#'   \item{m_convergence, lambda_star}{Sweeps per stored draw, and the
#'     convergence rate \eqn{\lambda^*}{lambda*} (both \code{NULL} for the
#'     exact iid route).}
#'   \item{n_obs, n_groups, group_name, groupef.names}{Design dimensions.}
#'   \item{varcor}{The reference \code{lmer} variance-correlation object.}
#'   \item{popef}{Named list, one entry per random-effect component, each a
#'     list of matrices \code{coefficients0} (\code{Prior Mean},
#'     \code{Prior.sd}, \code{lmer}, \code{lmer.se}), \code{coefficients}
#'     (\code{Post.Mode}, \code{Post.Mean}, \code{Post.Sd}, \code{MC Error},
#'     \code{Pr(Prior_tail)}) and \code{Percentiles}.}
#'   \item{popef_prior_overview, popef_overview, popef_percentiles_overview}{
#'     The same three tables stacked across components, with rows named
#'     \code{component::parameter}.}
#'   \item{tau2_prior_overview, tau2_overview, tau2_percentiles_overview}{
#'     Block~2 variance components \eqn{\tau^2_k}{tau^2_k}: prior family and
#'     moments against the \code{lmer} estimate, posterior summaries, and
#'     posterior percentiles. A \code{Residual} row for the observation
#'     variance \eqn{\sigma^2}{sigma^2} is appended when it is pooled.}
#'   \item{sigma2_group_overview}{Per-group \eqn{\sigma^2}{sigma^2} table when
#'     the measurement dispersion is per-group, otherwise \code{NULL}.}
#'   \item{groupef_overview}{Distribution of the Block~1 posterior modes
#'     across groups, per component.}
#'   \item{groupef_groups}{Per-group detail for the levels named in
#'     \code{groups}, or \code{NULL}.}
#'   \item{any_non_normal}{\code{TRUE} when any Block~2 component uses a
#'     non-\code{dNormal} prior.}
#'
#'   \code{print.summary.rlmerb} returns its argument invisibly and is called
#'   for the side effect of printing.
#'
#' @seealso \code{\link{rlmerb}}, \code{rglmerb()} in \pkg{lmebayesCore},
#'   \code{\link{model_setup}}, \code{\link{Prior_Setup_GLMM}},
#'   \code{\link[glmbayesCore]{summary.rglmb}}
#' @export
#' @method summary rlmerb
summary.rlmerb <- function(object, groups = NULL, ...) {

  if (!inherits(object, c("rlmerb", "rglmerb"))) {
    stop("'object' must be an rlmerb or rglmerb fit.", call. = FALSE)
  }

  design <- object$design
  prior  <- object$prior
  if (is.null(design) || is.null(prior)) {
    stop(
      "'object' is missing its design/prior components; summary() needs the ",
      "fit as returned by rlmerb()/rglmerb().",
      call. = FALSE
    )
  }

  re_names <- design$groupef.names
  n_draws  <- nrow(object$popef[[re_names[1L]]])
  ci       <- object$convergence_info

  popef_parts <- stats::setNames(
    lapply(re_names, function(k) {
      .lmebayes_popef_component_summary(object, k, n_draws = n_draws)
    }),
    re_names
  )

  res <- list(
    call            = object$call,
    formula         = design$formula,
    n               = n_draws,
    sim_method_used = if (!is.null(ci)) ci$sim_method_used else NULL,
    method          = if (!is.null(ci)) ci$method else NULL,
    draw_engine     = if (!is.null(ci)) ci$draw_engine else NULL,
    m_convergence   = object$m_convergence,
    lambda_star     = if (!is.null(ci)) ci$lambda_star else NULL,
    n_obs           = length(design$y),
    n_groups        = nlevels(design$group),
    group_name      = design$group_name,
    groupef.names   = re_names,
    varcor          = design$varcorr,
    dispersion_mode = prior$dispersion_mode,
    any_non_normal  = isTRUE(object$any_non_normal) ||
      isTRUE(prior$any_non_normal),

    popef                      = popef_parts,
    popef_prior_overview       = .lmebayes_stack_popef(popef_parts,
                                                       "coefficients0"),
    popef_overview             = .lmebayes_stack_popef(popef_parts,
                                                       "coefficients"),
    popef_percentiles_overview = .lmebayes_stack_popef(popef_parts,
                                                       "Percentiles"),

    tau2_prior_overview       = .lmebayes_tau2_prior_overview(object),
    tau2_overview             = .lmebayes_tau2_overview(object, n_draws),
    tau2_percentiles_overview = .lmebayes_tau2_percentiles_overview(object),
    sigma2_group_overview     = .lmebayes_sigma2_group_overview(object),

    groupef_overview = .lmebayes_groupef_overview(object)
  )

  if (!is.null(groups) && length(groups) > 0L) {
    res$groupef_groups <- .lmebayes_groupef_groups_detail(object, groups)
  }

  class(res) <- "summary.rlmerb"
  res
}

#' @rdname summary.rlmerb
#' @export
#' @method print summary.rlmerb
print.summary.rlmerb <- function(x, digits = max(3, getOption("digits") - 3),
                                 ...) {

  cat("Call\n")
  print(x$call)

  cat(sprintf(
    "\nBayesian mixed model fit based on %d draws [%s]\n",
    x$n, .lmebayes_engine_label(x)
  ))
  if (!is.null(x$formula)) {
    cat("Formula:", deparse1(x$formula), "\n")
  }
  cat(sprintf(
    "Number of obs: %d,  groups: %s, %d\n",
    x$n_obs, x$group_name, x$n_groups
  ))

  if (!is.null(x$varcor)) {
    cat("\nRandom effects (lmer reference):\n")
    print(x$varcor, comp = "Std.Dev.", digits = digits)
  }

  if (!is.null(x$tau2_prior_overview)) {
    cat("\n=== Block 2 variance components (tau^2_k) ===\n")
    cat("\nPrior and lmer reference\n\n")
    .lmebayes_print_table(x$tau2_prior_overview, digits = digits)

    ## With every component dNormal and a fixed sigma^2 nothing here is
    ## sampled, so the posterior tables would just repeat the fixed values.
    if (.lmebayes_tau2_is_sampled(x)) {
      cat("\nPosterior summaries\n\n")
      stats::printCoefmat(x$tau2_overview, digits = digits, quote = FALSE)
      cat("\nDistribution percentiles (tau^2)\n\n")
      stats::printCoefmat(x$tau2_percentiles_overview, digits = digits,
                          quote = FALSE)
    } else {
      cat("\n  (All variance components held fixed; nothing sampled.)\n")
    }
  }

  if (!is.null(x$sigma2_group_overview)) {
    cat("\n=== Block 1 measurement dispersion (sigma^2_j) ===\n\n")
    .lmebayes_print_table(x$sigma2_group_overview, digits = digits)
  }

  cat("\n=== Block 2 population effects (level-2 hyperparameters) ===\n")
  cat("\nPrior Estimates with Standard Deviations\n\n")
  stats::printCoefmat(x$popef_prior_overview, digits = digits, quote = FALSE,
                      na.print = "NA")
  cat("\nBayesian Estimates Based on", x$n, "draws\n\n")
  ## Pr(Prior_tail) is a prior/posterior conflict diagnostic, not a
  ## significance test, so it is printed without significance stars.
  stats::printCoefmat(x$popef_overview, digits = digits, quote = FALSE,
                      signif.stars = FALSE)
  cat("\nDistribution Percentiles\n\n")
  stats::printCoefmat(x$popef_percentiles_overview, digits = digits,
                      quote = FALSE)

  cat("\n=== Block 1 group effects ===\n")
  cat("\nDistribution of posterior modes across groups\n\n")
  stats::printCoefmat(x$groupef_overview, digits = digits, quote = FALSE)

  if (!is.null(x$groupef_groups)) {
    cat("\nPer-group detail (requested levels)\n\n")
    print(x$groupef_groups, digits = digits, row.names = FALSE)
  } else {
    cat(
      "\nPer-group effects: inspect fit$groupef.mode or fit$groupef, or call\n",
      "  summary(fit, groups = <level ids>) for selected groups.\n",
      sep = ""
    )
  }
  cat("\n")

  invisible(x)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' @noRd
.lmebayes_engine_label <- function(x) {
  if (identical(x$method, "exact_iid")) {
    return("exact iid draws")
  }
  sweeps <- if (!is.null(x$m_convergence)) {
    sprintf(", %d sweeps/draw", as.integer(x$m_convergence))
  } else {
    ""
  }
  paste0("two-block Gibbs", sweeps)
}

## TRUE when at least one variance component actually varies across draws.
#' @noRd
.lmebayes_tau2_is_sampled <- function(x) {
  if (is.null(x$tau2_overview) || is.null(x$tau2_percentiles_overview)) {
    return(FALSE)
  }
  sds <- x$tau2_overview[, "Post.Sd"]
  any(is.finite(sds) & sds > 0)
}

#' @noRd
.lmebayes_print_table <- function(tab, digits) {
  if (is.null(tab) || nrow(tab) == 0L) {
    return(invisible(tab))
  }
  out <- tab
  num <- vapply(out, is.numeric, logical(1L))
  out[num] <- lapply(out[num], round, digits = digits)
  print(out, right = TRUE)
  invisible(out)
}

#' Map a Block~2 hyperparameter onto the equivalent reference fixed effect
#'
#' Cross-level moderation columns (a level-2 covariate \code{par_name}
#' moderating the random slope \code{re_name}) correspond to the interaction
#' term, not to the level-2 main effect; \code{popef.moderation} records the
#' interaction column name built by \code{\link{model_setup}}.
#'
#' @noRd
.lmebayes_reference_fixef_lookup <- function(fit, re_name, par_name,
                                             moderation = NULL) {

  if (is.null(fit)) {
    return(list(estimate = NA_real_, se = NA_real_))
  }
  fe <- tryCatch(.lmebayes_reference_fixef(fit), error = function(e) NULL)
  if (is.null(fe)) {
    return(list(estimate = NA_real_, se = NA_real_))
  }

  candidates <- if (par_name == "(Intercept)" && re_name == "(Intercept)") {
    "(Intercept)"
  } else if (par_name == "(Intercept)") {
    re_name
  } else if (re_name == "(Intercept)") {
    par_name
  } else {
    from_map <- character(0)
    if (!is.null(moderation) && nrow(moderation) > 0L) {
      hit <- moderation$random_slope == re_name &
        moderation$moderator == par_name
      if (any(hit)) {
        from_map <- unique(moderation$interaction_col[hit])
      }
    }
    c(from_map, paste0(par_name, ":", re_name), paste0(re_name, ":", par_name))
  }

  hit <- candidates[candidates %in% names(fe)]
  if (length(hit) == 0L) {
    return(list(estimate = NA_real_, se = NA_real_))
  }

  nm <- hit[1L]
  se <- tryCatch({
    V <- .lmebayes_reference_vcov(fit)
    if (nm %in% rownames(V)) sqrt(V[nm, nm]) else NA_real_
  }, error = function(e) NA_real_)

  list(estimate = unname(fe[[nm]]), se = se)
}

#' @noRd
.lmebayes_popef_component_summary <- function(object, k, n_draws) {

  design <- object$design
  pl_k   <- object$prior$pop.prior_list[[k]]
  draws  <- object$popef[[k]]
  par    <- colnames(draws)

  prior_mean <- as.numeric(pl_k$mu)
  prior_sd   <- sqrt(diag(as.matrix(pl_k$Sigma)))

  ref <- lapply(par, function(nm) {
    .lmebayes_reference_fixef_lookup(
      design$lmer, k, nm, moderation = design$popef.moderation
    )
  })

  Tab0 <- cbind(
    "Prior Mean" = prior_mean,
    "Prior.sd"   = prior_sd,
    "lmer"       = vapply(ref, `[[`, numeric(1), "estimate"),
    "lmer.se"    = vapply(ref, `[[`, numeric(1), "se")
  )
  rownames(Tab0) <- par

  post_sd <- apply(draws, 2L, stats::sd)
  pval2 <- vapply(seq_along(par), function(j) {
    p1 <- mean(draws[, j] < prior_mean[j])
    min(p1, 1 - p1)
  }, numeric(1))

  TAB <- cbind(
    "Post.Mode"      = as.numeric(object$popef.mode[[k]]),
    "Post.Mean"      = as.numeric(object$popef.means[[k]]),
    "Post.Sd"        = post_sd,
    "MC Error"       = post_sd / sqrt(n_draws),
    "Pr(Prior_tail)" = pval2
  )
  rownames(TAB) <- par

  pct <- t(apply(draws, 2L, stats::quantile,
                 probs = c(0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99)))
  colnames(pct) <- c("1.0%", "2.5%", "5.0%", "Median",
                     "95.0%", "97.5%", "99.0%")
  rownames(pct) <- par

  list(coefficients0 = Tab0, coefficients = TAB, Percentiles = pct)
}

#' @noRd
.lmebayes_stack_popef <- function(parts, what) {
  rows <- lapply(names(parts), function(k) {
    tab <- parts[[k]][[what]]
    if (is.null(tab) || nrow(tab) == 0L) {
      return(NULL)
    }
    rownames(tab) <- paste0(k, "::", rownames(tab))
    tab
  })
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (length(rows) == 0L) {
    return(NULL)
  }
  do.call(rbind, rows)
}

## Prior family and moments for tau^2_k, next to the lmer variance component.
#' @noRd
.lmebayes_tau2_prior_overview <- function(object) {

  ptypes <- object$prior$ptypes
  if (is.null(ptypes)) {
    return(NULL)
  }
  re_names <- object$design$groupef.names
  Psi      <- object$design$Psi

  tab <- do.call(rbind, lapply(re_names, function(k) {
    pf <- object$prior$pfamily_list[[k]]
    pl <- if (!is.null(pf)) pf$prior_list else object$prior$pop.prior_list[[k]]

    if (identical(as.character(ptypes[[k]]), "dIndependent_Normal_Gamma")) {
      shape <- as.numeric(pl$shape[1L])
      rate  <- as.numeric(pl$rate[1L])
      E_tau2 <- if (is.finite(shape) && shape > 1 && is.finite(rate)) {
        rate / (shape - 1)
      } else {
        NA_real_
      }
      d_lo <- suppressWarnings(as.numeric(pl$disp_lower))
      d_hi <- suppressWarnings(as.numeric(pl$disp_upper))
    } else {
      E_tau2 <- as.numeric(object$prior$pop.prior_list[[k]]$dispersion)
      d_lo <- d_hi <- NA_real_
    }

    mer_tau2 <- if (!is.null(Psi) && k %in% names(Psi)) {
      unname(Psi[[k]])
    } else {
      NA_real_
    }

    data.frame(
      Prior           = as.character(ptypes[[k]]),
      `E[tau2]`       = E_tau2,
      `sqrt(E[tau2])` = if (is.finite(E_tau2) && E_tau2 >= 0) {
        sqrt(E_tau2)
      } else {
        NA_real_
      },
      disp_lower      = if (is.finite(d_lo)) d_lo else NA_real_,
      disp_upper      = if (is.finite(d_hi)) d_hi else NA_real_,
      lmer            = mer_tau2,
      `lmer SD`       = if (is.finite(mer_tau2) && mer_tau2 >= 0) {
        sqrt(mer_tau2)
      } else {
        NA_real_
      },
      check.names      = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  rownames(tab) <- re_names

  .lmebayes_append_residual_row(tab, object, kind = "prior")
}

#' @noRd
.lmebayes_tau2_overview <- function(object, n_draws) {

  re_names <- object$design$groupef.names
  td <- object$popef.dispersion
  if (is.null(td)) {
    return(NULL)
  }
  td <- td[, re_names, drop = FALSE]

  post_sd <- apply(td, 2L, stats::sd)
  out <- cbind(
    "Post.Mean" = colMeans(td),
    "Post.Sd"   = post_sd,
    "MC Error"  = post_sd / sqrt(n_draws),
    "Mean SD"   = apply(sqrt(td), 2L, mean)
  )
  rownames(out) <- re_names
  .lmebayes_append_residual_row(out, object, kind = "overview")
}

#' @noRd
.lmebayes_tau2_percentiles_overview <- function(object) {

  re_names <- object$design$groupef.names
  td <- object$popef.dispersion
  if (is.null(td)) {
    return(NULL)
  }
  td <- td[, re_names, drop = FALSE]

  pct <- t(apply(td, 2L, stats::quantile,
                 probs = c(0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99)))
  colnames(pct) <- c("1.0%", "2.5%", "5.0%", "Median",
                     "95.0%", "97.5%", "99.0%")
  rownames(pct) <- re_names
  .lmebayes_append_residual_row(pct, object, kind = "percentiles")
}

## The pooled observation variance sigma^2 shares the tau^2 tables' columns,
## so it is appended as a Residual row. Per-group dispersion has no single
## pooled value and is tabulated by .lmebayes_sigma2_group_overview() instead.
#' @noRd
.lmebayes_append_residual_row <- function(tab, object, kind) {

  mode <- object$prior$dispersion_mode
  if (is.null(tab) || is.null(mode) ||
      !mode %in% c("fixed", "gamma")) {
    return(tab)
  }

  sigma2 <- object$group.dispersion
  draws  <- if (identical(mode, "gamma") && !is.null(sigma2)) {
    as.numeric(sigma2)
  } else {
    NULL
  }
  fixed_val <- as.numeric(object$prior$group.dispersion)[1L]
  mer_sigma2 <- object$design$dispersion

  if (identical(kind, "prior")) {
    df <- data.frame(
      Prior           = if (identical(mode, "gamma")) "dGamma" else "fixed",
      `E[tau2]`       = if (is.null(draws)) fixed_val else mean(draws),
      `sqrt(E[tau2])` = sqrt(if (is.null(draws)) fixed_val else mean(draws)),
      disp_lower      = NA_real_,
      disp_upper      = NA_real_,
      lmer            = mer_sigma2,
      `lmer SD`       = if (is.finite(mer_sigma2) && mer_sigma2 >= 0) {
        sqrt(mer_sigma2)
      } else {
        NA_real_
      },
      check.names      = FALSE,
      stringsAsFactors = FALSE
    )
    rownames(df) <- "Residual"
    return(rbind(tab, df))
  }

  if (identical(kind, "overview")) {
    if (is.null(draws)) {
      out <- cbind("Post.Mean" = fixed_val, "Post.Sd" = 0,
                   "MC Error" = 0, "Mean SD" = sqrt(fixed_val))
    } else {
      sd_d <- stats::sd(draws)
      out <- cbind("Post.Mean" = mean(draws), "Post.Sd" = sd_d,
                   "MC Error" = sd_d / sqrt(length(draws)),
                   "Mean SD" = mean(sqrt(draws)))
    }
    rownames(out) <- "Residual"
    return(rbind(tab, out))
  }

  if (identical(kind, "percentiles")) {
    vals <- if (is.null(draws)) rep(fixed_val, 7L) else {
      unname(stats::quantile(
        draws, probs = c(0.01, 0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
      ))
    }
    out <- matrix(vals, nrow = 1L, dimnames = list("Residual", colnames(tab)))
    return(rbind(tab, out))
  }

  tab
}

## Per-group sigma^2_j, for dispformula models (estimated or known per group).
#' @noRd
.lmebayes_sigma2_group_overview <- function(object) {

  mode <- object$prior$dispersion_mode
  if (is.null(mode) || !mode %in% c("gamma_list", "fixed_vector")) {
    return(NULL)
  }

  mer_sigma2 <- object$design$dispersion
  sigma2 <- object$group.dispersion

  if (identical(mode, "fixed_vector") || !is.matrix(sigma2)) {
    vals <- as.numeric(object$prior$group.dispersion)
    grp  <- names(object$prior$group.dispersion)
    if (is.null(grp)) grp <- levels(object$design$group)
    out <- data.frame(
      `Fixed sigma2` = vals,
      lmer           = mer_sigma2,
      check.names      = FALSE,
      stringsAsFactors = FALSE
    )
    rownames(out) <- grp
    return(out)
  }

  grp <- colnames(sigma2)
  out <- data.frame(
    Post.Mean = colMeans(sigma2),
    Post.Sd   = apply(sigma2, 2L, stats::sd),
    `Mean SD` = colMeans(sqrt(sigma2)),
    `2.5%`    = apply(sigma2, 2L, stats::quantile, probs = 0.025),
    Median    = apply(sigma2, 2L, stats::median),
    `97.5%`   = apply(sigma2, 2L, stats::quantile, probs = 0.975),
    lmer      = mer_sigma2,
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )
  rownames(out) <- grp
  out
}

#' @noRd
.lmebayes_groupef_overview <- function(object) {

  re_names <- object$design$groupef.names
  b_mode   <- object$groupef.mode
  grp_col  <- object$design$group_name

  out <- t(vapply(re_names, function(k) {
    v <- b_mode[, k]
    mcmc_mean <- mean(tapply(object$groupef[[k]],
                             object$groupef[[grp_col]], mean))
    c(Mean = mean(v), SD = stats::sd(v), Min = min(v),
      Q1 = unname(stats::quantile(v, 0.25)),
      Median = stats::median(v),
      Q3 = unname(stats::quantile(v, 0.75)),
      Max = max(v), MCMC.mean = mcmc_mean)
  }, numeric(8)))

  rownames(out) <- re_names
  out
}

#' @noRd
.lmebayes_groupef_groups_detail <- function(object, groups) {

  re_names <- object$design$groupef.names
  grp_col  <- object$design$group_name
  groups   <- as.character(groups)

  known <- rownames(object$groupef.mode)
  bad <- setdiff(groups, known)
  if (length(bad)) {
    stop(
      "Unknown grouping level(s): ", paste(bad, collapse = ", "),
      call. = FALSE
    )
  }

  rows <- lapply(groups, function(lev) {
    out <- data.frame(group = lev, stringsAsFactors = FALSE)
    idx <- object$groupef[[grp_col]] == lev
    for (k in re_names) {
      out[[paste0(k, ".mode")]] <- object$groupef.mode[lev, k]
      out[[paste0(k, ".mean")]] <- mean(object$groupef[idx, k])
      out[[paste0(k, ".sd")]]   <- stats::sd(object$groupef[idx, k])
    }
    out
  })

  do.call(rbind, rows)
}
