## MER prior calibration helpers (from lmebayesCore mixed_rmerb_helpers.R; originals retained there).

#' Expand a scalar-or-length-\eqn{N} argument into a named length-\eqn{N} vector
#'
#' Shared "scalar or one value per named unit" resolver used for arguments
#' that can be supplied either as a single number (recycled to every unit in
#' \code{names_ref}) or as a length-\eqn{N} numeric vector (named, matching
#' \code{names_ref} in any order, or positional). Generalizes the \code{expand()}
#' closure in \code{.lmebayes_resolve_disp_prior()} and the \code{check_w()} /
#' recycling logic in \code{.lmebayes_resolve_measurement_disp_prior_group()}
#' so both Block~1 (per-group) and Block~2 (per-RE-component) "scalar or
#' vector" arguments -- and \code{dGamma_list()}'s override -- share one
#' validation path.
#' @noRd
.lmebayes_expand_scalar_or_vector <- function(x, names_ref, what,
                                               range = c(0.5, 1)) {
  n <- length(names_ref)
  check_range <- function(v) {
    if (!is.numeric(v) || anyNA(v) || any(v <= range[1]) || any(v >= range[2])) {
      stop(
        sprintf(
          "'%s' must be numeric with all values in (%s, %s).",
          what, range[1], range[2]
        ),
        call. = FALSE
      )
    }
  }

  if (length(x) == 1L) {
    check_range(x)
    v <- rep(as.numeric(x), n)
  } else if (length(x) == n) {
    check_range(x)
    v <- as.numeric(x)
    nms <- names(x)
    if (!is.null(nms) && any(nzchar(nms))) {
      if (!setequal(nms, names_ref)) {
        stop(
          sprintf(
            "Names of '%s' must match: %s.",
            what, paste(names_ref, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      names(v) <- nms
      v <- v[names_ref]
    }
  } else {
    stop(
      sprintf("'%s' must have length 1 or %d.", what, n),
      call. = FALSE
    )
  }

  stats::setNames(v, names_ref)
}

#' Resolve Block~1 \eqn{\sigma^2} prior weight into observation-scale \code{n_prior}.
#'
#' \eqn{n_{\mathrm{prior}} = w/(1-w)\times n} with \eqn{w =} \code{pwt_measurement}.
#' Independent of Block~2 fixef \code{pwt} and Block~2 \eqn{\tau^2}
#' \code{pwt_dispersion}.
#' @noRd
.lmebayes_resolve_measurement_disp_prior <- function(
    pwt_measurement,
    n_prior_measurement,
    n_obs
) {
  if (!is.null(pwt_measurement) && !is.null(n_prior_measurement)) {
    stop(
      "Supply at most one of 'pwt_measurement' and 'n_prior_measurement'.",
      call. = FALSE
    )
  }
  if (!is.numeric(n_obs) || length(n_obs) != 1L || !is.finite(n_obs) ||
      n_obs <= 0) {
    stop("'n_obs' must be a positive finite scalar.", call. = FALSE)
  }

  if (!is.null(pwt_measurement)) {
    if (!is.numeric(pwt_measurement) || length(pwt_measurement) != 1L) {
      stop(
        "'pwt_measurement' for the pooled Block~1 path must be a scalar; ",
        "supply a length-J vector for per-group calibration via ",
        "dGamma_list() only.",
        call. = FALSE
      )
    }
    if (is.na(pwt_measurement) || pwt_measurement <= 0 || pwt_measurement >= 1) {
      stop(
        "'pwt_measurement' must be a scalar in (0, 1).",
        call. = FALSE
      )
    }
    w <- as.numeric(pwt_measurement)
    n_prior <- w / (1 - w) * n_obs
    src <- "user-supplied (group.dispersion.pwt)"
  } else if (!is.null(n_prior_measurement)) {
    if (!is.numeric(n_prior_measurement) || length(n_prior_measurement) != 1L ||
        is.na(n_prior_measurement) || n_prior_measurement <= 0 ||
        !is.finite(n_prior_measurement)) {
      stop(
        "'n_prior_measurement' must be a positive finite scalar.",
        call. = FALSE
      )
    }
    n_prior <- as.numeric(n_prior_measurement)
    w <- n_prior / (n_prior + n_obs)
    src <- "user-supplied (group.dispersion.nprior)"
  } else {
    w <- 0.01
    n_prior <- w / (1 - w) * n_obs
    src <- "default (group.dispersion.pwt = 0.01)"
  }

  if (n_prior > n_obs) {
    stop(
      "Measurement dispersion prior requires n_prior <= n (equivalently ",
      "pwt_measurement <= 0.5); got n_prior = ", signif(n_prior, 4),
      ", n = ", n_obs, ".",
      call. = FALSE
    )
  }

  list(
    pwt_measurement     = w,
    n_prior_measurement = n_prior,
    source              = src
  )
}

#' Resolve per-group Block~1 \eqn{\sigma^2} prior weights into \code{n_prior_j}
#'
#' \eqn{n_{\mathrm{prior},j} = w_j/(1-w_j)\times n_j} for each group level.
#' When \code{pwt_measurement} is a scalar, the same weight applies to every
#' group.  When it is a length-\eqn{J} vector, names must match
#' \code{group_levels} if supplied.
#' @noRd
.lmebayes_resolve_measurement_disp_prior_group <- function(
    pwt_measurement,
    n_prior_measurement,
    n_j,
    group_levels
) {
  if (!is.null(pwt_measurement) && !is.null(n_prior_measurement)) {
    stop(
      "Supply at most one of 'pwt_measurement' and 'n_prior_measurement'.",
      call. = FALSE
    )
  }

  J <- length(group_levels)
  n_j <- as.integer(n_j)
  if (length(n_j) != J || anyNA(n_j) || any(n_j <= 0L)) {
    stop(
      "'n_j' must be a positive integer vector of length J (one per group level).",
      call. = FALSE
    )
  }
  names(n_j) <- group_levels

  check_w <- function(v, what) {
    if (!is.numeric(v) || anyNA(v) || any(v <= 0) || any(v >= 1)) {
      stop(sprintf("%s must be numeric with all values in (0, 1).", what),
           call. = FALSE)
    }
  }

  if (!is.null(pwt_measurement)) {
    if (length(pwt_measurement) == 1L) {
      check_w(pwt_measurement, "'pwt_measurement'")
      w <- rep(as.numeric(pwt_measurement), J)
    } else {
      if (length(pwt_measurement) != J) {
        stop(
          sprintf(
            "'pwt_measurement' vector must have length J = %d (number of group levels).",
            J
          ),
          call. = FALSE
        )
      }
      check_w(pwt_measurement, "'pwt_measurement'")
      w <- as.numeric(pwt_measurement)
      nms <- names(pwt_measurement)
      if (!is.null(nms) && any(nzchar(nms))) {
        if (!setequal(nms, group_levels)) {
          stop(
            "Names of 'pwt_measurement' must match group levels: ",
            paste(group_levels, collapse = ", "),
            call. = FALSE
          )
        }
        names(w) <- nms
        w <- w[group_levels]
      } else {
        names(w) <- group_levels
      }
    }
    n_prior <- w / (1 - w) * n_j
    src <- if (length(pwt_measurement) == 1L) {
      "user-supplied scalar (group.dispersion.pwt)"
    } else {
      "user-supplied vector (group.dispersion.pwt)"
    }
  } else {
    w <- rep(0.01, J)
    n_prior <- w / (1 - w) * n_j
    src <- if (!is.null(n_prior_measurement)) {
      "default per group (group.dispersion.pwt = 0.01; scalar group.dispersion.nprior applies to pooled path only)"
    } else {
      "default (group.dispersion.pwt = 0.01 per group)"
    }
  }

  names(w) <- group_levels
  names(n_prior) <- group_levels

  if (any(n_prior > n_j)) {
    bad <- names(n_prior)[n_prior > n_j]
    stop(
      "Per-group measurement dispersion prior requires n_prior_j <= n_j for every group; ",
      "failed for: ", paste(bad, collapse = ", "),
      call. = FALSE
    )
  }

  list(
    pwt_measurement     = w,
    n_prior_measurement = n_prior,
    source              = src
  )
}

#' Within-group Block~1 formula from random-coefficient names
#'
#' Per-group \eqn{\sigma^2} calibration (\code{\link[glmbayesCore]{Prior_Setup}} parity) fits
#' only predictors that enter the within-group likelihood---the population-mean
#' structure aligned with \code{design$groupef.names}.  Level-2 hyper covariates
#' and cross-level moderation terms in the full mixed-model formula are excluded.
#' @noRd
.lmebayes_block_formula_from_re <- function(formula, groupef.names) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (length(groupef.names) < 1L || anyNA(groupef.names)) {
    stop(
      "'groupef.names' must be a non-empty character vector.",
      call. = FALSE
    )
  }

  resp <- all.vars(formula)[1L]
  slope_terms <- setdiff(groupef.names, "(Intercept)")
  rhs <- if (length(slope_terms) == 0L) {
    "1"
  } else {
    paste(c("1", slope_terms), collapse = " + ")
  }

  stats::as.formula(paste(resp, "~", rhs))
}

#' Prior mean vector for block-formula Gaussian calibration (Prior_Setup parity)
#'
#' Matches \code{\link[glmbayesCore]{Prior_Setup}} defaults on a group subset:
#' intercept from an intercept-only \code{lm()} when
#' \code{intercept_source = "null_model"}, slopes zero when
#' \code{effects_source = "null_effects"}.
#' @noRd
.lmebayes_block_formula_prior_mu <- function(
    block_formula,
    dat_j,
    intercept_source = c("null_model", "full_model"),
    effects_source = c("null_effects", "full_model")
) {
  intercept_source <- match.arg(intercept_source)
  effects_source   <- match.arg(effects_source)

  X         <- stats::model.matrix(block_formula, data = dat_j)
  var_names <- colnames(X)
  mu        <- rep(0, length(var_names))
  names(mu) <- var_names

  if ("(Intercept)" %in% var_names) {
    if (intercept_source == "null_model") {
      resp <- all.vars(block_formula)[1L]
      null_fit <- stats::lm(
        stats::as.formula(paste(resp, "~ 1")),
        data = dat_j
      )
      mu["(Intercept)"] <- unname(stats::coef(null_fit)["(Intercept)"])
    } else {
      full_fit <- stats::lm(block_formula, data = dat_j)
      mu["(Intercept)"] <- unname(stats::coef(full_fit)["(Intercept)"])
    }
  }

  if (effects_source == "full_model") {
    full_fit <- stats::lm(block_formula, data = dat_j)
    for (nm in setdiff(var_names, "(Intercept)")) {
      mu[nm] <- unname(stats::coef(full_fit)[nm])
    }
  }

  matrix(mu, ncol = 1L, dimnames = list(var_names, "mu"))
}

#' Per-group Gaussian measurement-dispersion calibration (Block~1 dGamma density)
#'
#' Within-group glm inputs for per-group measurement-dispersion ING calibration.
#' @noRd
.lmebayes_ing_prior_measurement_group_glm_inputs <- function(
    lev,
    dat_j,
    block_formula,
    sd_tau,
    family = gaussian(),
    intercept_source = c("null_model", "full_model"),
    effects_source = c("null_effects", "full_model")
) {
  intercept_source <- match.arg(intercept_source)
  effects_source   <- match.arg(effects_source)

  mf <- stats::model.frame(block_formula, data = dat_j)
  X  <- stats::model.matrix(block_formula, data = dat_j)
  Y  <- stats::model.response(mf)
  var_names <- colnames(X)
  nvar <- ncol(X)
  n_j  <- nrow(X)
  weights <- rep(1, n_j)
  offset  <- rep(0, n_j)

  glm_full <- stats::glm.fit(
    x = X,
    y = Y,
    weights = weights,
    family = family
  )
  glm_full$weights <- weights
  class(glm_full) <- c("glm", "lm")

  V0 <- stats::vcov(glm_full)
  if (anyNA(V0)) {
    XtW <- sweep(X, 1, weights, `*`)
    Gm  <- crossprod(XtW, X)
    Ginv <- tryCatch(
      solve(Gm),
      error = function(e) {
        stop(
          "Group '", lev, "': vcov(glm) is NA and (X'WX) is singular.",
          call. = FALSE
        )
      }
    )
    res <- Y - X %*% coef(glm_full)
    rss <- sum(weights * res^2)
    if (n_j <= nvar || !is.finite(rss) || rss <= 0) {
      stop(
        "Group '", lev, "': cannot recover vcov for rank-deficient glm fit.",
        call. = FALSE
      )
    }
    d_v0 <- rss / (n_j - nvar)
    V0 <- d_v0 * Ginv
    dimnames(V0) <- list(var_names, var_names)
  }

  V0_diag <- diag(V0)
  if (any(V0_diag <= 0)) {
    stop(
      "Group '", lev, "': diagonal entries of V0 must be positive.",
      call. = FALSE
    )
  }

  sd_vec <- sd_tau[var_names]
  if (anyNA(sd_vec)) {
    stop(
      "Group '", lev, "': block_formula coefficients must align with sd_tau names.",
      call. = FALSE
    )
  }

  bhat <- coef(glm_full)
  res  <- residuals(glm_full, type = "response")
  rss  <- sum(weights * res^2)
  if (n_j <= nvar || !is.finite(rss) || rss <= 0) {
    stop(
      "Group '", lev, "': Gaussian dispersion requires n_j > p.",
      call. = FALSE
    )
  }
  dispersion_classical <- rss / (n_j - nvar)
  mu <- .lmebayes_block_formula_prior_mu(
    block_formula    = block_formula,
    dat_j            = dat_j,
    intercept_source = intercept_source,
    effects_source   = effects_source
  )

  list(
    X                    = X,
    Y                    = Y,
    weights              = weights,
    offset               = offset,
    V0                   = V0,
    bhat                 = bhat,
    dispersion_classical = dispersion_classical,
    mu_vec               = as.numeric(mu),
    var_names            = var_names,
    nvar                 = nvar,
    n_j                  = n_j,
    sd_vec               = sd_vec
  )
}

#' Pack \code{compute_gaussian_prior()} output for measurement-dispersion lists.
#' @noRd
.lmebayes_ing_prior_list_from_cal <- function(
    cal,
    n_prior_j,
    n_j,
    p_re,
    pwt_record,
    pwt_group_j
) {
  sh <- cal$shape_ING
  rt <- cal$rate_gamma
  list(
    sigma2_hat  = cal$dispersion,
    shape       = cal$shape,
    shape_ING   = sh,
    rate        = cal$rate,
    rate_gamma  = rt,
    E_sigma2    = if (is.finite(sh) && sh > 1 && is.finite(rt) && rt > 0) {
      rt / (sh - 1)
    } else {
      NA_real_
    },
    inv_E       = if (is.finite(sh) && sh > 0 && is.finite(rt) && rt > 0) {
      rt / sh
    } else {
      NA_real_
    },
    n_prior     = n_prior_j,
    n_j         = n_j,
    n_combined  = n_prior_j + n_j,
    p_re        = p_re,
    pwt         = pwt_record,
    pwt_group   = pwt_group_j
  )
}

#' \code{compute_gaussian_prior()} on within-group glm inputs and \code{Sigma}.
#' @noRd
.lmebayes_compute_ing_prior_cal_from_sigma <- function(inp, Sigma, n_prior_j) {
  Sigma_0 <- Sigma / inp$dispersion_classical
  glmbayesCore::compute_gaussian_prior(
    X           = inp$X,
    Y           = inp$Y,
    weights     = inp$weights,
    offset      = inp$offset,
    dispersion  = NULL,
    n_effective = inp$n_j,
    bhat        = inp$bhat,
    mu          = inp$mu_vec,
    Sigma_0     = Sigma_0,
    Sigma       = Sigma,
    n_prior     = n_prior_j,
    k           = 1
  )
}

#' Per-group Part-0 \code{Sigma_j} + Part VI \code{Omega_j}, then
#' \code{compute_gaussian_prior()} (shared by pooled aggregator and group path).
#' @noRd
.lmebayes_measurement_group_smarg_cal <- function(
    lev,
    dat_j,
    design,
    block_formula,
    sd_tau,
    prior_list,
    n_prior_j,
    family = gaussian(),
    intercept_source = "null_model",
    effects_source = "null_effects"
) {
  re_names <- design$groupef.names
  inp <- .lmebayes_ing_prior_measurement_group_glm_inputs(
    lev              = lev,
    dat_j            = dat_j,
    block_formula    = block_formula,
    sd_tau           = sd_tau,
    family           = family,
    intercept_source = intercept_source,
    effects_source   = effects_source
  )

  pwt_j <- diag(inp$V0)
  pwt_j <- pwt_j / (pwt_j + inp$sd_vec^2)
  names(pwt_j) <- inp$var_names

  if (length(pwt_j) == 1L) {
    Sigma <- ((1 - pwt_j) / pwt_j) * inp$V0
  } else {
    scale_vec <- sqrt((1 - pwt_j) / pwt_j)
    Sigma <- inp$V0 * outer(scale_vec, scale_vec)
  }

  Omega_j <- matrix(
    0, nrow = length(inp$var_names), ncol = length(inp$var_names),
    dimnames = list(inp$var_names, inp$var_names)
  )
  for (k in re_names) {
    Wk_row <- design$W[[k]][lev, , drop = FALSE]
    Sigma_k <- prior_list[[k]]$Sigma
    Omega_j[k, k] <- as.numeric(Wk_row %*% Sigma_k %*% t(Wk_row))
  }

  cal <- .lmebayes_compute_ing_prior_cal_from_sigma(
    inp, Sigma + Omega_j, n_prior_j
  )
  list(inp = inp, cal = cal, Sigma = Sigma, Omega_j = Omega_j, pwt_j = pwt_j)
}

#' Aggregate per-group \eqn{S_{\mathrm{marg},j}} into a pooled \eqn{\hat\sigma^2}.
#'
#' \eqn{\hat\sigma^2_{\mathrm{pool}} =
#' (\sum_j S_{\mathrm{marg},j}) / (n - J p_{\mathrm{re}})} under the
#' block-diagonal joint (independent group priors). Uses a dummy
#' \code{n_prior,j} only to call \code{compute_gaussian_prior()}; the
#' returned center does not depend on that weight.
#' @noRd
.lmebayes_pooled_measurement_smarg_from_groups <- function(
    design,
    data,
    block_formula,
    sd_tau,
    prior_list,
    group_levels,
    family = gaussian(),
    intercept_source = "null_model",
    effects_source = "null_effects"
) {
  p_re <- length(design$groupef.names)
  if (p_re < 1L) {
    stop(
      "Pooled measurement S_marg aggregation requires at least one random coefficient.",
      call. = FALSE
    )
  }
  J <- length(group_levels)
  n <- length(design$y)

  S_marg <- stats::setNames(numeric(J), group_levels)
  sigma2_j <- stats::setNames(numeric(J), group_levels)
  n_j_vec <- stats::setNames(integer(J), group_levels)
  df_j <- stats::setNames(numeric(J), group_levels)

  for (lev in group_levels) {
    idx <- design$group == lev
    dat_j <- data[idx, , drop = FALSE]
    n_j <- sum(idx)
    n_j_vec[[lev]] <- as.integer(n_j)
    ## Dummy n_prior: S_marg / sigma2_hat do not depend on it.
    n_prior_j <- max(1e-8, 0.01 / 0.99 * n_j)
    piece <- tryCatch(
      .lmebayes_measurement_group_smarg_cal(
        lev              = lev,
        dat_j            = dat_j,
        design           = design,
        block_formula    = block_formula,
        sd_tau           = sd_tau,
        prior_list       = prior_list,
        n_prior_j        = n_prior_j,
        family           = family,
        intercept_source = intercept_source,
        effects_source   = effects_source
      ),
      error = function(e) NULL
    )
    if (is.null(piece)) {
      ## Within-group RE design singular (e.g. slope constant in group):
      ## fall back to intercept-only residual SS for that group's contribution.
      Y <- as.numeric(design$y[idx])
      if (n_j < 2L) {
        stop(
          "Pooled measurement S_marg aggregation: group '", lev,
          "' has n_j < 2 and a singular within-group design.",
          call. = FALSE
        )
      }
      S_marg[[lev]] <- sum((Y - mean(Y))^2)
      df_j[[lev]] <- n_j - 1
      sigma2_j[[lev]] <- S_marg[[lev]] / df_j[[lev]]
    } else {
      df_j[[lev]] <- piece$inp$n_j - p_re
      sigma2_j[[lev]] <- piece$cal$dispersion
      S_marg[[lev]] <- piece$cal$dispersion * df_j[[lev]]
    }
  }

  den <- sum(df_j)
  if (!is.finite(den) || den <= 0) {
    stop(
      "Pooled measurement S_marg aggregation produced a non-positive residual df.",
      call. = FALSE
    )
  }

  sigma2_pool <- sum(S_marg) / den
  if (!is.finite(sigma2_pool) || sigma2_pool <= 0) {
    stop(
      "Pooled measurement S_marg aggregation produced a non-positive center.",
      call. = FALSE
    )
  }

  list(
    sigma2_pool = as.numeric(sigma2_pool),
    S_marg      = S_marg,
    sigma2_j    = sigma2_j,
    n_j         = n_j_vec,
    df_j        = df_j,
    n           = n,
    J           = J,
    p_re        = p_re,
    den         = den
  )
}

#' Block~2 hyper-regression marginal \eqn{\hat\tau^2_k} via
#' \code{compute_gaussian_prior()} (A12 Sec. 3.3.4).
#'
#' Treats group-level coefficients \eqn{b_{\cdot k}} as the response for
#' \eqn{b_{\cdot k}\sim N(W_k\gamma_k,\tau^2_k I_J)} with prior
#' \eqn{\gamma_k\sim N(\mu_k,\Sigma_k)}.
#' @noRd
.lmebayes_compute_ing_prior_cal_tau2_hyper <- function(
    Y,
    X,
    mu,
    Sigma,
    n_prior,
    k_arg = 1
) {
  Y <- as.numeric(Y)
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  J <- length(Y)
  p_k <- ncol(X)
  if (nrow(X) != J) {
    stop(
      "Block~2 tau2 marginal: nrow(X) must equal length(Y) (= J).",
      call. = FALSE
    )
  }
  if (J <= p_k) {
    stop(
      "Block~2 tau2 marginal requires J > p_k; got J = ", J,
      ", p_k = ", p_k, ".",
      call. = FALSE
    )
  }
  if (!is.numeric(mu) || length(mu) != p_k || any(!is.finite(mu))) {
    stop("Block~2 tau2 marginal: 'mu' must be a finite length-p_k vector.",
         call. = FALSE)
  }
  if (!is.matrix(Sigma) || nrow(Sigma) != p_k || ncol(Sigma) != p_k ||
      anyNA(Sigma)) {
    stop("Block~2 tau2 marginal: 'Sigma' must be a finite [p_k x p_k] matrix.",
         call. = FALSE)
  }
  if (!is.numeric(n_prior) || length(n_prior) != 1L || !is.finite(n_prior) ||
      n_prior <= 0) {
    stop("Block~2 tau2 marginal: 'n_prior' must be a positive finite scalar.",
         call. = FALSE)
  }

  weights <- rep(1, J)
  offset <- rep(0, J)
  fit <- stats::lm.fit(x = X, y = Y)
  bhat <- as.numeric(fit$coefficients)
  if (any(!is.finite(bhat))) {
    stop("Block~2 tau2 marginal: OLS fit for b.k on W_k failed.",
         call. = FALSE)
  }
  names(bhat) <- colnames(X)
  rss <- sum((Y - as.numeric(X %*% bhat))^2)
  if (!is.finite(rss) || rss <= 0) {
    stop(
      "Block~2 tau2 marginal: classical RSS for b.k on W_k must be positive.",
      call. = FALSE
    )
  }
  dispersion_classical <- rss / (J - p_k)
  Sigma_0 <- Sigma / dispersion_classical

  glmbayesCore::compute_gaussian_prior(
    X           = X,
    Y           = Y,
    weights     = weights,
    offset      = offset,
    dispersion  = NULL,
    n_effective = J,
    bhat        = bhat,
    mu          = as.numeric(mu),
    Sigma_0     = Sigma_0,
    Sigma       = Sigma,
    n_prior     = n_prior,
    k           = k_arg
  )
}

#' Extract group-level coefficient vector \eqn{b_{\cdot k}} from a reference fit.
#' @noRd
.lmebayes_reference_b_component <- function(
    fit,
    group_name,
    component,
    group_levels
) {
  co <- .lmebayes_reference_coef(fit)[[group_name]]
  if (is.null(co)) {
    stop(
      "Could not find grouping factor '", group_name, "' in coef(fit).",
      call. = FALSE
    )
  }
  if (!component %in% colnames(co)) {
    stop(
      "coef(fit)[[\"", group_name, "\"]] is missing column '", component, "'.",
      call. = FALSE
    )
  }
  rn <- rownames(co)
  if (is.null(rn)) {
    stop(
      "coef(fit)[[\"", group_name, "\"]] must have row names (group levels).",
      call. = FALSE
    )
  }
  if (!setequal(rn, group_levels)) {
    stop(
      "coef(fit) group levels do not match design group levels.",
      call. = FALSE
    )
  }
  stats::setNames(as.numeric(co[group_levels, component]), group_levels)
}

#' Per-group Block~1 measurement-dispersion calibration for \code{dGamma_list()}.
#'
#' \code{sigma2_hat}, \code{shape_ING}, and \code{rate_gamma} from shared
#' population \code{sd_tau} coefficient shrinkage (\eqn{V_0} scaled by
#' per-coefficient \code{pwt_j}). Also stores \code{rate} (A12 3.3.4
#' \eqn{S_{\mathrm{marg}}}) for dev comparison only.
#'
#' Also folds in the Part VI extension of
#' \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md} -- "also integrating out
#' the prior mean \code{mu_j}" -- with a model-derived \code{Omega_j}: each
#' RE component's \code{mu_j[k]} stands in for the model's own conditional
#' prior mean \code{W_j[[k]] \%*\% gamma_k} (the Block~2 hyper-regression),
#' and \code{gamma_k}'s own calibrated uncertainty (\code{prior_list[[k]]$Sigma},
#' already computed above in \code{Prior_Setup_GLMM()}) propagates
#' through group \code{j}'s own hyper-design row,
#' \code{Omega_j[k, k] = W_j[[k]] \%*\% Sigma_k \%*\% t(W_j[[k]])}
#' (diagonal across RE components -- each \code{gamma_k} is calibrated
#' independently). \code{Sigma_j' = Sigma_j + Omega_j} is then used in place
#' of \code{Sigma_j} for \code{compute_gaussian_prior()}, so the resulting
#' \code{rate}/\code{sigma2_hat} integrate out both \code{b_j} (random
#' effects, as before) and \code{gamma} (fixed effects, via \code{Omega_j}).
#' This is now the permanent default (not opt-in) -- see
#' \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md} Part VI.
#'
#' \code{disp_lower}/\code{disp_upper} are computed here too, as literal
#' quantiles of this same \code{Gamma(shape_ING, rate)} marginal (via
#' \code{\link{.lmebayes_ing_prior_quantile_window}}) -- the same construction
#' already used for the pooled \code{ing_prior_measurement} case -- rather
#' than \code{dGamma_list()}'s former, decoupled \code{n_combined}-based
#' window.
#' @noRd
.lmebayes_calibrate_ing_prior_measurement_group <- function(
    design,
    data,
    block_formula,
    sd_tau,
    pwt_group,
    n_prior_group,
    group_levels,
    prior_list,
    max_disp_perc_group,
    family = gaussian(),
    intercept_source = c("null_model", "full_model"),
    effects_source = c("null_effects", "full_model")
) {
  intercept_source <- match.arg(intercept_source)
  effects_source   <- match.arg(effects_source)
  p_re <- length(design$groupef.names)
  if (p_re < 1L) {
    stop(
      "Per-group measurement dispersion calibration requires at least one random coefficient.",
      call. = FALSE
    )
  }
  if (length(sd_tau) != p_re || anyNA(sd_tau) || any(sd_tau <= 0)) {
    stop(
      "'sd_tau' must be a named numeric vector of positive RE standard deviations.",
      call. = FALSE
    )
  }
  stats::setNames(
    lapply(group_levels, function(lev) {
      idx   <- design$group == lev
      dat_j <- data[idx, , drop = FALSE]
      n_prior_j <- unname(n_prior_group[[lev]])

      piece <- .lmebayes_measurement_group_smarg_cal(
        lev              = lev,
        dat_j            = dat_j,
        design           = design,
        block_formula    = block_formula,
        sd_tau           = sd_tau,
        prior_list       = prior_list,
        n_prior_j        = n_prior_j,
        family           = family,
        intercept_source = intercept_source,
        effects_source   = effects_source
      )
      inp <- piece$inp
      cal <- piece$cal

      .ing_stop_if_prior_exceeds_data(
        shape       = cal$shape_ING,
        p           = inp$nvar,
        n_w         = inp$n_j,
        detail      = paste0("group '", lev, "' has n_j = ", inp$n_j),
        limit_label = "n_j",
        prefix      = "Per-group measurement dispersion: "
      )

      ## The window must bound the POSTERIOR spread the sampler's own
      ## envelope machinery actually draws sigma2_j from (EnvelopeDispersionBuild.cpp:
      ## shape2 = Shape + n_w/2), not the prior alone -- shape_ING,j (cal$shape_ING)
      ## is the PRIOR shape (n_prior,j-only) fed to the sampler as-is (unchanged
      ## below); only the window's own shape/rate are inflated by n_j/2, mean-matched
      ## at the same sigma2_hat,j. See inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md Part II
      ## (algebraically shape_post,j == shape_w,j = (n_combined,j+1)/2 + p_re/2).
      mdp_j <- unname(max_disp_perc_group[[lev]])
      shape_post_j <- cal$shape_ING + inp$n_j / 2
      rate_post_j  <- cal$dispersion * (shape_post_j - 1)
      win <- .lmebayes_ing_prior_quantile_window(
        shape_post_j, rate_post_j, mdp_j
      )

      out <- .lmebayes_ing_prior_list_from_cal(
        cal         = cal,
        n_prior_j   = n_prior_j,
        n_j         = inp$n_j,
        p_re        = p_re,
        pwt_record  = piece$pwt_j,
        pwt_group_j = unname(pwt_group[[lev]])
      )
      out$disp_lower    <- win$disp_lower
      out$disp_upper    <- win$disp_upper
      out$max_disp_perc <- mdp_j
      out$omega_j       <- piece$Omega_j
      out
    }),
    group_levels
  )
}

#' Print \code{rate_gamma} (A12 3.3.5, downstream) vs \code{rate} (A12 3.3.4
#' \eqn{S_{\mathrm{marg}}}) from the same per-group calibration.
#' @noRd
.lmebayes_print_ing_prior_measurement_group_compare <- function(
    existing,
    digits = 4
) {
  if (is.null(existing)) {
    return(invisible(NULL))
  }
  grp <- names(existing)
  if (is.null(grp) || length(grp) < 1L) {
    return(invisible(NULL))
  }

  pct_diff <- function(new, old) {
    if (!is.finite(old) || old == 0) {
      return(NA_real_)
    }
    100 * (new - old) / old
  }

  inv_E_rate <- function(sh, rt) {
    if (is.finite(sh) && sh > 0 && is.finite(rt) && rt > 0) {
      rt / sh
    } else {
      NA_real_
    }
  }

  tab <- do.call(rbind, lapply(grp, function(g) {
    ex <- existing[[g]]
    rt334 <- ex$rate
    data.frame(
      group          = g,
      n_j            = ex$n_j,
      n_prior        = ex$n_prior,
      shape_ING      = ex$shape_ING,
      rate_gamma     = ex$rate_gamma,
      rate           = rt334,
      pct_rate       = pct_diff(rt334, ex$rate_gamma),
      inv_E_gamma    = ex$inv_E,
      inv_E_rate     = inv_E_rate(ex$shape_ING, rt334),
      pct_inv_E      = pct_diff(
        inv_E_rate(ex$shape_ING, rt334),
        ex$inv_E
      ),
      stringsAsFactors = FALSE
    )
  }))
  rownames(tab) <- NULL

  cat(
    "\n--- Per-group Block~1 gamma: rate_gamma (A12 3.3.5) vs rate (A12 3.3.4 S_marg) ---\n",
    "  dGamma_list downstream uses rate (S_marg); rate_gamma retained for comparison.\n",
    "  sigma2_hat and truncation bounds unchanged.\n\n",
    sep = ""
  )
  num_cols <- vapply(tab, is.numeric, logical(1))
  if (any(num_cols)) {
    tab[num_cols] <- lapply(tab[num_cols], round, digits = digits)
  }
  print(tab, row.names = FALSE)
  invisible(tab)
}

#' Central 98% prior-mass \eqn{\sigma^2}/\eqn{\tau^2} window from calibrated precision prior
#'
#' Precision \eqn{1/\sigma^2 \sim \mathrm{Gamma}(\code{shape}, \code{rate})};
#' bounds are 0.01/0.99 quantiles inverted to the variance scale.
#' @noRd
.lmebayes_ing_prior_quantile_window <- function(shape, rate, max_disp_perc = 0.99) {
  if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
    stop(
      "ING prior quantile window requires positive finite shape and rate.",
      call. = FALSE
    )
  }
  list(
    disp_lower = 1 / stats::qgamma(max_disp_perc,       shape = shape, rate = rate),
    disp_upper = 1 / stats::qgamma(1 - max_disp_perc,   shape = shape, rate = rate)
  )
}

#' Prospective \code{dGamma()} measurement \eqn{\sigma^2} calibration from setup
#'
#' Mean-matched inverse-Gamma hyperparameters for Block~1 ING (same algebra as
#' \code{ing_prior} for \eqn{\tau^2_k}, with \eqn{\hat\sigma^2} =
#' \code{group.dispersion}, \eqn{p = p_{\mathrm{re}}}, and
#' \eqn{n_{\mathrm{prior}} = \mathrm{pwt\_measurement}/(1-\mathrm{pwt\_measurement})\times n} on the total
#' observation count).  Truncation bounds are the central 98% prior-mass
#' interval from the same \code{shape}/\code{rate}.
#' @noRd
.lmebayes_calibrate_ing_prior_measurement <- function(
    design,
    group.dispersion,
    n_prior,
    max_disp_perc = 0.99
) {
  p_re <- length(design$groupef.names)
  n    <- length(design$y)
  if (p_re < 1L) {
    stop(
      "Measurement dispersion calibration requires at least one random coefficient.",
      call. = FALSE
    )
  }

  if (!is.numeric(n_prior) || length(n_prior) != 1L || !is.finite(n_prior) ||
      n_prior <= 0) {
    stop(
      "'n_prior' must be a positive finite scalar for measurement dispersion calibration.",
      call. = FALSE
    )
  }
  if (n_prior > n) {
    stop(
      "Measurement dispersion prior requires n_prior <= n; got n_prior = ",
      signif(n_prior, 4), ", n = ", n, ".",
      call. = FALSE
    )
  }

  shape <- (n_prior + 1) / 2 + p_re / 2
  rate  <- as.numeric(group.dispersion) * (n_prior + p_re - 1) / 2
  if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
    stop(
      "Measurement dispersion ING calibration produced non-positive shape/rate.",
      call. = FALSE
    )
  }

  win <- .lmebayes_ing_prior_quantile_window(shape, rate, max_disp_perc)

  list(
    sigma2_hat    = as.numeric(group.dispersion),
    shape         = shape,
    rate          = rate,
    disp_lower    = win$disp_lower,
    disp_upper    = win$disp_upper,
    max_disp_perc = max_disp_perc,
    n_prior       = n_prior,
    n_effective = n,
    p_re        = p_re
  )
}

#' Limiting-posterior \eqn{\sigma^2}/\eqn{\tau^2} truncation window (lmebayes default)
#'
#' Central 98% mass of \code{Gamma((J+1)/2, d_hat*(J-1)/2)} inverted to the
#' variance scale; see \code{inst/ING_TRUNCATION_WINDOW.md}.
#' @noRd
.lmebayes_ing_limiting_posterior_window <- function(d_hat, J, max_disp_perc = 0.99) {
  if (!is.numeric(d_hat) || length(d_hat) != 1L || !is.finite(d_hat) ||
      d_hat <= 0) {
    stop(
      "'d_hat' must be a positive finite scalar (classical variance plug-in).",
      call. = FALSE
    )
  }
  J <- as.integer(J[1L])
  if (!is.finite(J) || J < 1L) {
    stop("'J' must be a positive integer (number of groups).", call. = FALSE)
  }
  a_inf <- (J + 1) / 2
  b_inf <- as.numeric(d_hat) * (J - 1) / 2
  if (b_inf <= 0) {
    stop(
      "Limiting-posterior ING window requires J >= 2 (got J = ", J, ").",
      call. = FALSE
    )
  }
  list(
    disp_lower = 1 / stats::qgamma(max_disp_perc,     shape = a_inf, rate = b_inf),
    disp_upper = 1 / stats::qgamma(1 - max_disp_perc, shape = a_inf, rate = b_inf)
  )
}
