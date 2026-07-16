## Internal helpers for matrix-level mixed-model samplers (rlmerb / rglmerb).
## lmebayes calls these via glmbayesCore::: from lmerb() / glmerb() only.

#' @noRd
.lmebayes_resolve_dispersion_ranef <- function(
    dispersion_ranef,
    family,
    design = NULL,
    fn_name = "lmerb"
) {
  has_dispersion <- family$family %in%
    c("gaussian", "Gamma", "quasipoisson", "quasibinomial")

  if (!has_dispersion) {
    if (!is.null(dispersion_ranef)) {
      stop(
        "'dispersion_ranef' must be NULL for family = ", family$family,
        "() (no observation-level dispersion).",
        call. = FALSE
      )
    }
    return(list(
      mode                   = "none",
      dispersion_fix         = NULL,
      dispersion_prior_list  = NULL,
      dispersion_pfamily     = NULL
    ))
  }

  if (is.list(dispersion_ranef) && !inherits(dispersion_ranef, "pfamily")) {
    return(.lmebayes_resolve_dispersion_ranef_group_list(
      dispersion_ranef = dispersion_ranef,
      design           = design,
      fn_name          = fn_name
    ))
  }

  if (inherits(dispersion_ranef, "pfamily")) {
    if (!identical(dispersion_ranef$pfamily, "dGamma")) {
      stop(
        fn_name, "(): 'dispersion_ranef' pfamily must be dGamma(); got ",
        dispersion_ranef$pfamily, ". RE priors belong in 'pfamily_list'.",
        call. = FALSE
      )
    }
    pl <- dispersion_ranef$prior_list
    if (!isTRUE(pl$Inv_Dispersion)) {
      stop(
        fn_name, "(): dGamma() observation-dispersion prior requires ",
        "Inv_Dispersion = TRUE.",
        call. = FALSE
      )
    }
    shape <- as.numeric(pl$shape[1L])
    rate  <- as.numeric(pl$rate[1L])
    if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
      stop(
        fn_name, "(): dGamma() dispersion_ranef prior_list requires positive ",
        "'shape' and 'rate'.",
        call. = FALSE
      )
    }
    return(list(
      mode                  = "gamma",
      dispersion_fix        = shape / rate,
      dispersion_prior_list = pl,
      dispersion_pfamily    = dispersion_ranef
    ))
  }

  if (is.null(dispersion_ranef) || !is.numeric(dispersion_ranef) ||
      length(dispersion_ranef) != 1L || !is.finite(dispersion_ranef) ||
      dispersion_ranef <= 0) {
    stop(
      "'dispersion_ranef' must be a single positive number, a dGamma() ",
      "pfamily, or a named list of dGamma() pfamilies (one per group) for ",
      "family = ", family$family, "().",
      call. = FALSE
    )
  }
  list(
    mode                  = "fixed",
    dispersion_fix        = as.numeric(dispersion_ranef),
    dispersion_prior_list = NULL,
    dispersion_pfamily    = NULL
  )
}

#' Resolve a per-group list of \code{dGamma()} pfamilies for \code{dispersion_ranef}
#'
#' Third \code{dispersion_ranef} option alongside a fixed scalar and a single
#' (pooled) \code{dGamma()}: a named list with one \code{dGamma()} pfamily per
#' group level. Each entry keeps its own \code{shape}/\code{rate}/\code{disp_lower}/
#' \code{disp_upper} -- there is no requirement that groups share the same
#' \code{shape}/\code{rate} (\code{Prior_Setup_lmebayes()} may choose to build
#' them from a shared hyperprior, but the engine itself is agnostic).
#' Requires every group to be full column rank (rank-deficient groups are not
#' yet supported for this option).
#' @noRd
.lmebayes_resolve_dispersion_ranef_group_list <- function(
    dispersion_ranef,
    design,
    fn_name = "lmerb"
) {
  if (is.null(design) || is.null(design$groups)) {
    stop(
      fn_name, "(): a list of dGamma() priors for 'dispersion_ranef' requires ",
      "'design' with grouping information.",
      call. = FALSE
    )
  }
  group_levels <- levels(design$groups)
  J <- length(group_levels)

  if (length(dispersion_ranef) != J) {
    stop(
      fn_name, "(): 'dispersion_ranef' is a list of length ",
      length(dispersion_ranef), " but there are ", J, " group level(s) (",
      paste(group_levels, collapse = ", "), "). Supply exactly one dGamma() ",
      "pfamily per group.",
      call. = FALSE
    )
  }
  nms <- names(dispersion_ranef)
  if (is.null(nms) || any(!nzchar(nms)) || !setequal(nms, group_levels)) {
    stop(
      fn_name, "(): names(dispersion_ranef) must match the group levels (",
      paste(group_levels, collapse = ", "), ") exactly.",
      call. = FALSE
    )
  }
  window_diagnostics <- attr(dispersion_ranef, "window_diagnostics")
  dispersion_ranef <- dispersion_ranef[group_levels]

  if (is.null(design$re_rank) || !all(design$re_rank[group_levels])) {
    deficient <- if (!is.null(design$re_rank)) {
      group_levels[!design$re_rank[group_levels]]
    } else {
      group_levels
    }
    stop(
      fn_name, "(): a list of per-group dGamma() dispersion priors currently ",
      "requires every group to be full column rank; rank-deficient group(s): ",
      paste(deficient, collapse = ", "), ". Use a single dGamma() or a fixed ",
      "scalar 'dispersion_ranef' for models with rank-deficient groups.",
      call. = FALSE
    )
  }

  shape_group      <- stats::setNames(numeric(J), group_levels)
  rate_group       <- stats::setNames(numeric(J), group_levels)
  disp_lower_group <- stats::setNames(numeric(J), group_levels)
  disp_upper_group <- stats::setNames(numeric(J), group_levels)

  for (lev in group_levels) {
    pf <- dispersion_ranef[[lev]]
    if (!inherits(pf, "pfamily") || !identical(pf$pfamily, "dGamma")) {
      stop(
        fn_name, "(): dispersion_ranef[[\"", lev, "\"]] must be a dGamma() ",
        "pfamily.",
        call. = FALSE
      )
    }
    pl <- pf$prior_list
    if (!isTRUE(pl$Inv_Dispersion)) {
      stop(
        fn_name, "(): dispersion_ranef[[\"", lev, "\"]] dGamma() prior ",
        "requires Inv_Dispersion = TRUE.",
        call. = FALSE
      )
    }
    shape <- as.numeric(pl$shape[1L])
    rate  <- as.numeric(pl$rate[1L])
    if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
      stop(
        fn_name, "(): dispersion_ranef[[\"", lev, "\"]] must have positive, ",
        "finite 'shape' and 'rate'.",
        call. = FALSE
      )
    }
    lo <- pl$disp_lower
    hi <- pl$disp_upper
    if (is.null(lo) || is.null(hi) ||
        !is.numeric(lo) || !is.numeric(hi) ||
        length(lo) != 1L || length(hi) != 1L ||
        !is.finite(lo) || !is.finite(hi) ||
        lo <= 0 || hi <= lo) {
      stop(
        fn_name, "(): dispersion_ranef[[\"", lev, "\"]] must supply finite ",
        "'disp_lower' and 'disp_upper' with 0 < disp_lower < disp_upper -- a ",
        "list of dGamma() priors requires explicit per-group truncation bounds.",
        call. = FALSE
      )
    }
    shape_group[[lev]]      <- shape
    rate_group[[lev]]       <- rate
    disp_lower_group[[lev]] <- as.numeric(lo)
    disp_upper_group[[lev]] <- as.numeric(hi)
  }

  list(
    mode                  = "gamma_list",
    dispersion_fix        = mean(shape_group / rate_group),
    dispersion_prior_list = list(
      shape_group      = shape_group,
      rate_group       = rate_group,
      disp_lower_group = disp_lower_group,
      disp_upper_group = disp_upper_group
    ),
    dispersion_pfamily    = dispersion_ranef,
    window_diagnostics    = window_diagnostics
  )
}

#' @noRd
.lmebayes_validate_dispersion_ranef <- function(
    dispersion_ranef,
    family,
    fn_name = "lmerb"
) {
  resolved <- .lmebayes_resolve_dispersion_ranef(
    dispersion_ranef = dispersion_ranef,
    family           = family,
    design           = NULL,
    fn_name          = fn_name
  )
  resolved$dispersion_fix
}

#' @noRd
.lmebayes_priors_from_pfamily_list <- function(pfamily_list,
                                               dispersion_ranef,
                                               design,
                                               family,
                                               fn_name = "lmerb") {

  re_names <- design$re_coef_names
  p_re     <- length(re_names)

  ## --- dispersion_ranef (Block 1 measurement dispersion) -------------------
  disp_res <- .lmebayes_resolve_dispersion_ranef(
    dispersion_ranef = dispersion_ranef,
    family           = family,
    design           = design,
    fn_name          = fn_name
  )
  dispersion_ranef <- disp_res$dispersion_fix

  ## --- pfamily_list ---------------------------------------------------------
  if (!is.list(pfamily_list) || length(pfamily_list) != p_re) {
    stop(
      "'pfamily_list' must be a list with one pfamily per random-effect ",
      "component (", p_re, " expected: ", paste(re_names, collapse = ", "),
      "). Build it with pfamily_list(Prior_Setup_lmebayes(...)).",
      call. = FALSE
    )
  }
  if (is.null(names(pfamily_list)) || !setequal(names(pfamily_list), re_names)) {
    stop(
      "Names of 'pfamily_list' must match the random-effect coefficient ",
      "names: ", paste(re_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  pfamily_list <- pfamily_list[re_names]

  prior_list <- stats::setNames(vector("list", p_re), re_names)
  tau2   <- stats::setNames(numeric(p_re), re_names)
  ptypes <- stats::setNames(character(p_re), re_names)

  for (k in re_names) {
    pf <- pfamily_list[[k]]
    if (!inherits(pf, "pfamily")) {
      stop("pfamily_list[[\"", k, "\"]] must be a pfamily object.",
           call. = FALSE)
    }
    if (!pf$pfamily %in% c("dNormal", "dIndependent_Normal_Gamma")) {
      stop(
        fn_name, "() supports only dNormal and dIndependent_Normal_Gamma ",
        "pfamilies in 'pfamily_list'; component \"", k, "\" is ",
        pf$pfamily, ".",
        call. = FALSE
      )
    }
    ptypes[[k]] <- pf$pfamily

    par_names <- colnames(design$X_hyper[[k]])
    q_k <- length(par_names)

    mu_k <- as.numeric(pf$prior_list$mu)
    if (length(mu_k) != q_k) {
      stop(
        "pfamily_list[[\"", k, "\"]]$prior_list$mu has length ",
        length(mu_k), " but the hyper design has ", q_k, " column(s): ",
        paste(par_names, collapse = ", "), ".",
        call. = FALSE
      )
    }
    mu_nms <- rownames(pf$prior_list$mu)
    if (!is.null(mu_nms) && all(nzchar(mu_nms))) {
      if (!setequal(mu_nms, par_names)) {
        stop(
          "Parameter names of pfamily_list[[\"", k, "\"]] (",
          paste(mu_nms, collapse = ", "), ") do not match the hyper design ",
          "columns (", paste(par_names, collapse = ", "), ").",
          call. = FALSE
        )
      }
      ord <- match(par_names, mu_nms)
      mu_k <- mu_k[ord]
      Sigma_k <- as.matrix(pf$prior_list$Sigma)[ord, ord, drop = FALSE]
    } else {
      Sigma_k <- as.matrix(pf$prior_list$Sigma)
    }
    names(mu_k) <- par_names
    dimnames(Sigma_k) <- list(par_names, par_names)

    ## Keep the pfamily object itself aligned with the hyper-design column
    ## order: it is passed straight to the v2 sampler as the Block 2 source
    ## of truth, so its mu/Sigma must match x_hyper[[k]].
    pfamily_list[[k]]$prior_list$mu <-
      matrix(mu_k, ncol = 1L, dimnames = list(par_names, NULL))
    pfamily_list[[k]]$prior_list$Sigma <- Sigma_k

    if (identical(pf$pfamily, "dNormal")) {
      d_k <- pf$prior_list$dispersion
      if (isTRUE(pf$prior_list$ddef)) {
        warning(
          fn_name, ": pfamily_list[[\"", k, "\"]] uses the default ",
          "dispersion = 1 (none was supplied to dNormal()); the Block 1 ",
          "random-effect variance tau^2 for \"", k, "\" is therefore 1.",
          call. = FALSE
        )
      }
    } else {
      ## ING: disp_lower/disp_upper fix the truncation window and lambda*
      ## calibration only; ICM plug-in tau^2 comes from the pfamily spec.
      d_k <- pf$prior_list$disp_lower
      if (is.null(d_k) || !is.numeric(d_k) || length(d_k) != 1L ||
          !is.finite(d_k) || d_k <= 0) {
        stop(
          fn_name, "(): pfamily_list[[\"", k, "\"]] is ",
          "dIndependent_Normal_Gamma and must supply a positive scalar ",
          "'disp_lower' (lower dispersion truncation) for lambda* calibration.",
          call. = FALSE
        )
      }
      u_k <- pf$prior_list$disp_upper
      if (is.null(u_k) || !is.numeric(u_k) || length(u_k) != 1L ||
          !is.finite(u_k) || u_k <= as.numeric(d_k)) {
        stop(
          fn_name, "(): pfamily_list[[\"", k, "\"]] is ",
          "dIndependent_Normal_Gamma and must supply a finite scalar ",
          "'disp_upper' > 'disp_lower' (upper dispersion truncation), so ",
          "the tau^2 truncation window is fixed across Gibbs sweeps. ",
          "pfamily_list(Prior_Setup_lmebayes(...)) sets both bounds to the ",
          "0.01/0.99 prior dispersion quantiles by default.",
          call. = FALSE
        )
      }
    }

    tau2_k <- .two_block_tau2_plug_in_from_pfamily(pf)
    tau2[[k]] <- tau2_k
    prior_list[[k]] <- list(
      mu_fixef         = mu_k,
      Sigma_fixef      = Sigma_k,
      dispersion_fixef = tau2_k
    )
  }

  Sigma_ranef <- diag(unname(tau2), nrow = p_re, ncol = p_re)
  dimnames(Sigma_ranef) <- list(re_names, re_names)

  list(
    pfamily_list          = pfamily_list,
    dispersion_ranef      = dispersion_ranef,
    dispersion_mode       = disp_res$mode,
    dispersion_pfamily    = disp_res$dispersion_pfamily,
    dispersion_prior_list = disp_res$dispersion_prior_list,
    window_diagnostics    = disp_res$window_diagnostics,
    Sigma_ranef           = Sigma_ranef,
    prior_list            = prior_list,
    ptypes         = ptypes,
    any_non_normal = any(ptypes != "dNormal")
  )
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
    src <- "user-supplied (pwt_measurement)"
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
    src <- "user-supplied (n_prior_measurement)"
  } else {
    w <- 0.01
    n_prior <- w / (1 - w) * n_obs
    src <- "default (pwt_measurement = 0.01)"
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
        w <- w[group_levels]
      } else {
        names(w) <- group_levels
      }
    }
    n_prior <- w / (1 - w) * n_j
    src <- if (length(pwt_measurement) == 1L) {
      "user-supplied scalar (pwt_measurement)"
    } else {
      "user-supplied vector (pwt_measurement)"
    }
  } else {
    w <- rep(0.01, J)
    n_prior <- w / (1 - w) * n_j
    src <- if (!is.null(n_prior_measurement)) {
      "default per group (pwt_measurement = 0.01; scalar n_prior_measurement applies to pooled path only)"
    } else {
      "default (pwt_measurement = 0.01 per group)"
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
#' Per-group \eqn{\sigma^2} calibration (\code{\link{Prior_Setup}} parity) fits
#' only predictors that enter the within-group likelihood---the population-mean
#' structure aligned with \code{design$re_coef_names}.  Level-2 hyper covariates
#' and cross-level moderation terms in the full mixed-model formula are excluded.
#' @noRd
.lmebayes_block_formula_from_re <- function(formula, re_coef_names) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (length(re_coef_names) < 1L || anyNA(re_coef_names)) {
    stop(
      "'re_coef_names' must be a non-empty character vector.",
      call. = FALSE
    )
  }

  resp <- all.vars(formula)[1L]
  slope_terms <- setdiff(re_coef_names, "(Intercept)")
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
  compute_gaussian_prior(
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

#' \code{compute_gaussian_prior()} on within-group glm inputs, widened by a
#' fixed \code{Omega} to also integrate out uncertainty about the prior mean
#' (\code{mu_j -> mu_bar}); see
#' \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md} Part VI. \code{Omega} is
#' fixed and absolute (same coefficient-scale units as \code{Sigma}), so it is
#' divided by the same group-specific \code{dispersion_classical} that
#' \code{Sigma} already is -- no rescaling by the (unknown) working
#' dispersion is involved. Dev-only: not called from any exported path.
#' @noRd
.lmebayes_compute_full_marginal_ing_prior_cal_from_sigma <- function(
    inp, Sigma, n_prior_j, mu_bar, Omega
) {
  Sigma_0  <- Sigma / inp$dispersion_classical
  Omega_0  <- Omega / inp$dispersion_classical
  compute_gaussian_prior(
    X           = inp$X,
    Y           = inp$Y,
    weights     = inp$weights,
    offset      = inp$offset,
    dispersion  = NULL,
    n_effective = inp$n_j,
    bhat        = inp$bhat,
    mu          = mu_bar,
    Sigma_0     = Sigma_0 + Omega_0,
    Sigma       = Sigma + Omega,
    n_prior     = n_prior_j,
    k           = 1
  )
}

#' Per-group Block~1 measurement-dispersion calibration for \code{dGamma_list()}.
#'
#' \code{sigma2_hat}, \code{shape_ING}, and \code{rate_gamma} from shared
#' population \code{sd_tau} coefficient shrinkage (\eqn{V_0} scaled by
#' per-coefficient \code{pwt_j}). Also stores \code{rate} (A12 3.3.4
#' \eqn{S_{\mathrm{marg}}}) for dev comparison only.
#' @noRd
.lmebayes_calibrate_ing_prior_measurement_group <- function(
    design,
    data,
    block_formula,
    sd_tau,
    pwt_group,
    n_prior_group,
    group_levels,
    family = gaussian(),
    intercept_source = c("null_model", "full_model"),
    effects_source = c("null_effects", "full_model")
) {
  intercept_source <- match.arg(intercept_source)
  effects_source   <- match.arg(effects_source)
  p_re <- length(design$re_coef_names)
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
      idx   <- design$groups == lev
      dat_j <- data[idx, , drop = FALSE]
      n_prior_j <- unname(n_prior_group[[lev]])

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

      cal <- .lmebayes_compute_ing_prior_cal_from_sigma(
        inp, Sigma, n_prior_j
      )

      .ing_stop_if_prior_exceeds_data(
        shape       = cal$shape_ING,
        p           = inp$nvar,
        n_w         = inp$n_j,
        detail      = paste0("group '", lev, "' has n_j = ", inp$n_j),
        limit_label = "n_j",
        prefix      = "Per-group measurement dispersion: "
      )

      .lmebayes_ing_prior_list_from_cal(
        cal         = cal,
        n_prior_j   = n_prior_j,
        n_j         = inp$n_j,
        p_re        = p_re,
        pwt_record  = pwt_j,
        pwt_group_j = unname(pwt_group[[lev]])
      )
    }),
    group_levels
  )
}

#' Part VI ("full marginal") per-group measurement-dispersion calibration:
#' same as \code{.lmebayes_calibrate_ing_prior_measurement_group()} but
#' additionally integrates out \code{fixef} against its own \code{dNormal}
#' prior (Block~2), via a fixed \code{Omega} widening \code{Sigma_0j} and
#' \code{mu_j -> mu_bar}. \code{mu_bar}/\code{Omega} are read directly off
#' \code{prior_list} (\code{Prior_Setup_lmebayes()}'s own Block~2 fixef prior,
#' one \code{mu_fixef}/\code{Sigma_fixef} pair per RE component \code{k}) --
#' population-level and shared across all groups, not re-derived from pooled
#' data. Assembled block-diagonally onto \code{design$re_coef_names}; requires
#' each RE component's own \code{X_hyper[[k]]} to be intercept-only (the usual
#' random-intercept/random-slope case with no Block~2-level covariates).
#' \code{shape_ING} is unchanged (depends only on \code{n_prior_j}, \code{p_re}).
#' Called from \code{Prior_Setup_lmebayes()}, which stores its result as the
#' production \code{ing_prior_measurement_group} consumed by
#' \code{dGamma_list()} (falling back to
#' \code{.lmebayes_calibrate_ing_prior_measurement_group()}'s Part~I result on
#' failure); see \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md} Part VI.
#' @noRd
.lmebayes_calibrate_ing_prior_measurement_group_full_marginal <- function(
    design,
    data,
    block_formula,
    sd_tau,
    pwt_group,
    n_prior_group,
    group_levels,
    prior_list,
    family = gaussian(),
    intercept_source = c("null_model", "full_model"),
    effects_source = c("null_effects", "full_model")
) {
  intercept_source <- match.arg(intercept_source)
  effects_source   <- match.arg(effects_source)
  re_names <- design$re_coef_names
  p_re <- length(re_names)
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
  if (!identical(names(prior_list), re_names) && !setequal(names(prior_list), re_names)) {
    stop(
      "'prior_list' must have one element per RE component in ",
      "design$re_coef_names.",
      call. = FALSE
    )
  }

  ## Anchor: mu_bar / Omega are fixef's own dNormal prior (Block~2),
  ## population-level and shared across every group -- not re-derived from
  ## pooled data. Requires each RE component's own Block~2 hyper-design to be
  ## intercept-only (X_hyper[[k]] has exactly one column); there is no
  ## established mapping from block_formula's per-group coefficient k onto a
  ## multi-column Sigma_fixef[[k]] otherwise.
  mu_bar <- stats::setNames(numeric(p_re), re_names)
  Omega  <- matrix(0, p_re, p_re, dimnames = list(re_names, re_names))
  for (k in re_names) {
    Xh_k <- design$X_hyper[[k]]
    if (is.null(Xh_k) || ncol(Xh_k) != 1L) {
      stop(
        "Part VI full-marginal calibration requires an intercept-only ",
        "Block~2 hyper-design for RE component '", k, "' (no ",
        "Block~2-level covariates); ncol(X_hyper[['", k, "']]) != 1.",
        call. = FALSE
      )
    }
    pk <- prior_list[[k]]
    mu_bar[[k]] <- unname(pk$mu_fixef[1L])
    Omega[k, k] <- unname(pk$Sigma_fixef[1L, 1L])
  }

  stats::setNames(
    lapply(group_levels, function(lev) {
      idx   <- design$groups == lev
      dat_j <- data[idx, , drop = FALSE]
      n_prior_j <- unname(n_prior_group[[lev]])

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

      cal <- .lmebayes_compute_full_marginal_ing_prior_cal_from_sigma(
        inp, Sigma, n_prior_j, mu_bar, Omega
      )

      .ing_stop_if_prior_exceeds_data(
        shape       = cal$shape_ING,
        p           = inp$nvar,
        n_w         = inp$n_j,
        detail      = paste0("group '", lev, "' has n_j = ", inp$n_j),
        limit_label = "n_j",
        prefix      = "Per-group measurement dispersion (full marginal): "
      )

      .lmebayes_ing_prior_list_from_cal(
        cal         = cal,
        n_prior_j   = n_prior_j,
        n_j         = inp$n_j,
        p_re        = p_re,
        pwt_record  = pwt_j,
        pwt_group_j = unname(pwt_group[[lev]])
      )
    }),
    group_levels
  )
}

#' Print \code{rate_gamma} (A12 3.3.5) vs \code{rate} (A12 3.3.4
#' \eqn{S_{\mathrm{marg}}}) from the Part~I ("classical") per-group
#' calibration -- diagnostic only; \code{dGamma_list()} normally uses the
#' Part~VI "full marginal" calibration instead (see
#' \code{.lmebayes_print_ing_prior_measurement_group_full_marginal_compare()}).
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
    "\n--- Part I 'classical' calibration: rate_gamma (A12 3.3.5) vs\n",
    "    rate (A12 3.3.4 S_marg) ---\n",
    "  Both retained for comparison only; dGamma_list() normally uses the\n",
    "  Part VI 'full marginal' rate/sigma2_hat below instead (falls back to\n",
    "  this Part I 'rate' on failure). sigma2_hat/bounds shown here are\n",
    "  Part I's own (pre-Part VI) values.\n\n",
    sep = ""
  )
  num_cols <- vapply(tab, is.numeric, logical(1))
  if (any(num_cols)) {
    tab[num_cols] <- lapply(tab[num_cols], round, digits = digits)
  }
  print(tab, row.names = FALSE)
  invisible(tab)
}

#' Print the Part VI ("full marginal") \code{sigma2_hat}/\code{rate} vs the
#' Part~I (\code{"classical"}, A12 3.3.4) per-group calibration, plus the
#' implied vector of Part VI \code{sigma2_hat} values -- the dispersion point
#' estimates \code{dGamma_list()}'s window is now centered on (Part VI is the
#' production calibration as of the swap described in
#' \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md} Part VI). Diagnostic
#' print only; does not itself affect \code{Prior_Setup_lmebayes()}'s
#' returned object (the swap happens separately, in
#' \code{Prior_Setup_lmebayes()}'s own body).
#' @noRd
.lmebayes_print_ing_prior_measurement_group_full_marginal_compare <- function(
    existing,
    full_marginal,
    digits = 4
) {
  if (is.null(existing) || is.null(full_marginal)) {
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

  tab <- do.call(rbind, lapply(grp, function(g) {
    ex <- existing[[g]]
    fm <- full_marginal[[g]]
    data.frame(
      group             = g,
      n_j               = ex$n_j,
      n_prior           = ex$n_prior,
      sigma2_hat        = ex$sigma2_hat,
      sigma2_hat_full   = fm$sigma2_hat,
      pct_sigma2_hat    = pct_diff(fm$sigma2_hat, ex$sigma2_hat),
      rate              = ex$rate,
      rate_full         = fm$rate,
      pct_rate          = pct_diff(fm$rate, ex$rate),
      stringsAsFactors  = FALSE
    )
  }))
  rownames(tab) <- NULL

  cat(
    "\n--- Per-group Block~1 gamma: Part I 'classical' rate/sigma2_hat\n",
    "    (A12 3.3.4) vs Part VI 'full marginal' (mu_j -> mu_bar, integrates\n",
    "    out fixef via Block~2's own dNormal prior) ---\n",
    "  dGamma_list() now uses Part VI (falls back to Part I on failure); see\n",
    "  inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md Part VI.\n\n",
    sep = ""
  )
  num_cols <- vapply(tab, is.numeric, logical(1))
  if (any(num_cols)) {
    tab[num_cols] <- lapply(tab[num_cols], round, digits = digits)
  }
  print(tab, row.names = FALSE)

  implied_disp <- stats::setNames(
    vapply(grp, function(g) full_marginal[[g]]$sigma2_hat, numeric(1)),
    grp
  )
  cat("\nImplied Part VI dispersion values (sigma2_hat_full, per group):\n")
  print(round(implied_disp, digits))

  invisible(list(compare = tab, sigma2_hat_full = implied_disp))
}

#' BLUP/OLS residual RSS inflation ratio per group (for dGamma upper bounds)
#' @noRd
.lmebayes_group_blup_rss_inflation <- function(
    data,
    block_formula,
    fit_ref,
    groups,
    group_levels,
    group_name
) {
  beta_blup <- .lmebayes_reference_coef(fit_ref)[[group_name]]
  if (is.null(beta_blup)) {
    stop(
      "Reference fit has no random-effect levels for group '", group_name, "'.",
      call. = FALSE
    )
  }

  infl <- vapply(group_levels, function(lev) {
    idx     <- groups == lev
    dat_lev <- data[idx, , drop = FALSE]
    Xg      <- stats::model.matrix(block_formula, data = dat_lev)
    yg      <- stats::model.response(stats::model.frame(block_formula, data = dat_lev))
    beta_lev <- as.numeric(beta_blup[lev, colnames(Xg), drop = FALSE])
    resid_g  <- yg - Xg %*% beta_lev
    RSS_blup <- sum(resid_g^2)
    RSS_ols  <- sum(stats::residuals(stats::lm(block_formula, data = dat_lev))^2)
    if (!is.finite(RSS_blup) || !is.finite(RSS_ols) || RSS_ols <= 0) {
      stop(
        "Group '", lev, "': non-finite or non-positive OLS RSS for BLUP inflation.",
        call. = FALSE
      )
    }
    if (RSS_blup < RSS_ols - sqrt(.Machine$double.eps) * max(RSS_ols, 1)) {
      stop(
        "Group '", lev, "': BLUP RSS (", RSS_blup, ") < OLS RSS (", RSS_ols,
        "); inflation ratio must be >= 1.",
        call. = FALSE
      )
    }
    RSS_blup / RSS_ols
  }, numeric(1))

  stats::setNames(infl, group_levels)
}

#' Per-group envelope-centering dispersion estimate (for \code{dGamma_list()}
#' \code{disp_center = "dispersion2"})
#'
#' Pure-R replica of the \code{EnvelopeCentering()} trace-correction fixed
#' point (\code{src/EnvelopeCentering.cpp}), run \strong{without} a dispersion
#' prior contribution so the group's own \code{n_j} observations are used
#' exactly once (see "double-counting pitfall",
#' \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md} Part III). At a working
#' dispersion \code{dispersion2}, the posterior mean of \code{b_j} under
#' \code{b_j ~ N(mu_j, Sigma_ranef)} and the Gaussian likelihood is the
#' ridge/GLS estimator \code{b2 = (X'X/dispersion2 + P)^{-1} (X'Y/dispersion2
#' + P mu_j)}; the expected RSS under its posterior
#' (\code{RSS_precomputed = ||Y - X b2||^2 + tr(X'X Cov(b2))}) updates
#' \code{dispersion2 <- RSS_precomputed / (n_j - p)} to a 10-iteration (by
#' default) fixed point.
#' @noRd
.lmebayes_group_dispersion2_envelope_centering <- function(
    data,
    block_formula,
    Sigma_ranef,
    groups,
    group_levels,
    intercept_source = c("null_model", "full_model"),
    effects_source   = c("null_effects", "full_model"),
    n_rss_iter = 10L
) {
  intercept_source <- match.arg(intercept_source)
  effects_source   <- match.arg(effects_source)

  P  <- solve(Sigma_ranef)
  RA <- chol(P)

  out <- vapply(group_levels, function(lev) {
    idx   <- groups == lev
    dat_j <- data[idx, , drop = FALSE]
    X <- stats::model.matrix(block_formula, data = dat_j)
    Y <- stats::model.response(stats::model.frame(block_formula, data = dat_j))
    n_j <- nrow(X)
    p   <- ncol(X)

    if (n_j <= p) {
      stop(
        "Group '", lev, "': disp_center = \"dispersion2\" requires n_j > p ",
        "(n_j = ", n_j, ", p = ", p, ").",
        call. = FALSE
      )
    }

    mu_j <- .lmebayes_block_formula_prior_mu(
      block_formula    = block_formula,
      dat_j            = dat_j,
      intercept_source = intercept_source,
      effects_source   = effects_source
    )
    mu_vec <- as.numeric(mu_j)

    fit0 <- stats::lm.fit(X, Y)
    rss0 <- sum(fit0$residuals^2)
    dispersion2 <- rss0 / (n_j - p)

    z_bot <- RA %*% mu_vec
    XtX   <- t(X) %*% X

    for (iter in seq_len(n_rss_iter)) {
      s <- 1 / sqrt(dispersion2)
      W <- rbind(s * X, RA)
      z <- c(s * Y, z_bot)
      Sigma_post  <- solve(t(W) %*% W)
      b2          <- Sigma_post %*% (t(W) %*% z)
      r           <- Y - X %*% b2
      rss_at_mean <- sum(r^2)
      trace_term  <- sum(diag(XtX %*% Sigma_post))
      RSS_precomputed <- rss_at_mean + trace_term
      dispersion2 <- RSS_precomputed / (n_j - p)
    }

    dispersion2
  }, numeric(1))

  stats::setNames(out, group_levels)
}

#' \eqn{\sigma^2} CDF under precision \eqn{1/\sigma^2 \sim \mathrm{Gamma}(\code{shape}, \code{rate})}
#' @noRd
.lmebayes_dgamma_sigma2_pct_under_rate <- function(sigma2, shape, rate) {
  if (!is.finite(sigma2) || sigma2 <= 0 ||
      !is.finite(shape) || shape <= 0 ||
      !is.finite(rate) || rate <= 0) {
    return(NA_real_)
  }
  1 - stats::pgamma(1 / sigma2, shape = shape, rate = rate)
}

#' Cross-percentiles for asymmetric per-group \code{dGamma()} truncation windows
#'
#' \code{disp_lower} is the nominal lower quantile under the OLS-matched rate;
#' \code{disp_upper} under the BLUP-scaled upper rate. Returns relative-tail
#' ratios \code{R_lo = lo_pct_BLUP/lo_pct_OLS} and
#' \code{R_hi = hi_pct_OLS/hi_pct_BLUP}.
#' @noRd
.lmebayes_dgamma_window_cross_percentiles <- function(
    shape,
    rate_w,
    rate_u,
    max_disp_perc = 0.99,
    blup_infl = NA_real_,
    sigma2_hat = NA_real_
) {
  win <- .lmebayes_ing_prior_quantile_window_asymmetric(
    shape         = shape,
    rate_lower    = rate_w,
    rate_upper    = rate_u,
    max_disp_perc = max_disp_perc
  )
  lo <- win$disp_lower
  hi <- win$disp_upper
  lo_pct_ols  <- 100 * (1 - max_disp_perc)
  hi_pct_blup <- 100 * max_disp_perc
  lo_pct_blup <- 100 * .lmebayes_dgamma_sigma2_pct_under_rate(lo, shape, rate_u)
  hi_pct_ols  <- 100 * .lmebayes_dgamma_sigma2_pct_under_rate(hi, shape, rate_w)
  R_lo <- lo_pct_blup / lo_pct_ols
  R_hi <- hi_pct_ols / hi_pct_blup
  list(
    disp_lower  = lo,
    disp_upper  = hi,
    lo_pct_OLS  = lo_pct_ols,
    lo_pct_BLUP = lo_pct_blup,
    hi_pct_BLUP = hi_pct_blup,
    hi_pct_OLS  = hi_pct_ols,
    R_lo        = R_lo,
    R_hi        = R_hi,
    blup_infl   = blup_infl,
    sigma2_hat  = sigma2_hat
  )
}

#' Flag asymmetric per-group \code{dGamma()} truncation windows
#' @noRd
.lmebayes_dgamma_window_asymmetric_flag <- function(
    R_lo,
    R_hi,
    asymmetric_R_lo = 0.25,
    asymmetric_R_hi = 4
) {
  (is.finite(R_lo) && R_lo < asymmetric_R_lo) ||
    (is.finite(R_hi) && R_hi > asymmetric_R_hi)
}

#' Warn when per-group \code{dGamma()} window cross-percentiles indicate asymmetry
#' @noRd
.lmebayes_warn_dgamma_window_asymmetry <- function(
    diag,
    asymmetric_R_lo = 0.25,
    asymmetric_R_hi = 4,
    print_table = FALSE
) {
  if (is.null(diag) || !is.data.frame(diag) || nrow(diag) < 1L) {
    return(invisible(diag))
  }
  flagged <- diag[
    !is.na(diag$asymmetric_window) & diag$asymmetric_window,
    ,
    drop = FALSE
  ]
  if (!nrow(flagged)) {
    return(invisible(diag))
  }
  warning(
    sprintf(
      paste0(
        "%d group(s) with asymmetric dGamma truncation window ",
        "(R_lo = lo_pct_BLUP/lo_pct_OLS < %g or ",
        "R_hi = hi_pct_OLS/hi_pct_BLUP > %g): %s."
      ),
      nrow(flagged),
      asymmetric_R_lo,
      asymmetric_R_hi,
      paste(flagged$group, collapse = ", ")
    ),
    call. = FALSE
  )
  if (isTRUE(print_table)) {
    show_cols <- c(
      "group", "blup_infl", "R_lo", "R_hi", "lo_pct_BLUP", "hi_pct_OLS"
    )
    show_cols <- intersect(show_cols, names(flagged))
    out <- flagged[, show_cols, drop = FALSE]
    num_cols <- vapply(out, is.numeric, logical(1L))
    out[num_cols] <- lapply(out[num_cols], function(x) round(x, 3))
    cat("\n--- dGamma window asymmetry (flagged groups) ---\n")
    print(out, row.names = FALSE)
  }
  invisible(diag)
}

#' Approximate-posterior truncation window for per-group Block~1 \eqn{\sigma^2}
#' @noRd
.lmebayes_ing_prior_quantile_window_asymmetric <- function(
    shape,
    rate_lower,
    rate_upper,
    max_disp_perc = 0.99
) {
  win_lo <- .lmebayes_ing_prior_quantile_window(shape, rate_lower, max_disp_perc)
  list(
    disp_lower = win_lo$disp_lower,
    disp_upper = 1 / stats::qgamma(
      1 - max_disp_perc,
      shape = shape,
      rate  = rate_upper
    )
  )
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
#' \code{dispersion_ranef}, \eqn{p = p_{\mathrm{re}}}, and
#' \eqn{n_{\mathrm{prior}} = \mathrm{pwt\_measurement}/(1-\mathrm{pwt\_measurement})\times n} on the total
#' observation count).  Truncation bounds are the central 98% prior-mass
#' interval from the same \code{shape}/\code{rate}.
#' @noRd
.lmebayes_calibrate_ing_prior_measurement <- function(
    design,
    dispersion_ranef,
    n_prior,
    max_disp_perc = 0.99
) {
  p_re <- length(design$re_coef_names)
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
  rate  <- dispersion_ranef * (n_prior + p_re - 1) / 2
  if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
    stop(
      "Measurement dispersion ING calibration produced non-positive shape/rate.",
      call. = FALSE
    )
  }

  win <- .lmebayes_ing_prior_quantile_window(shape, rate, max_disp_perc)

  list(
    sigma2_hat    = dispersion_ranef,
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
#' variance scale; see \code{inst/ING_TRUNCATION_WINDOW.md} in \pkg{lmebayes}.
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

#' Build shared ING Block~1 measurement \code{prior_list} for lmebayes glue
#' @noRd
.lmebayes_ing_measurement_prior_list <- function(prior, disp_info, design) {
  re_names <- design$re_coef_names
  p_re     <- length(re_names)
  pl       <- disp_info$dispersion_prior_list
  if (is.null(pl$shape) || is.null(pl$rate)) {
    stop(
      "dGamma() dispersion_ranef prior_list must contain 'shape' and 'rate'.",
      call. = FALSE
    )
  }
  mu <- matrix(
    0,
    nrow = p_re,
    ncol = 1L,
    dimnames = list(re_names, NULL)
  )
  Sigma <- as.matrix(prior$Sigma_ranef)
  if (nrow(Sigma) != p_re || ncol(Sigma) != p_re) {
    stop(
      "prior$Sigma_ranef must be ", p_re, " x ", p_re, ".",
      call. = FALSE
    )
  }
  out <- list(
    mu            = mu,
    Sigma         = Sigma,
    shape         = pl$shape,
    rate          = pl$rate,
    max_disp_perc = if (!is.null(pl$max_disp_perc)) pl$max_disp_perc else 0.99
  )
  if (!is.null(pl$disp_lower)) out$disp_lower <- pl$disp_lower
  if (!is.null(pl$disp_upper)) out$disp_upper <- pl$disp_upper
  if (is.null(out$disp_lower) || is.null(out$disp_upper)) {
    win <- .lmebayes_ing_prior_quantile_window(
      shape = as.numeric(pl$shape[1L]),
      rate  = as.numeric(pl$rate[1L])
    )
    if (is.null(out$disp_lower)) out$disp_lower <- win$disp_lower
    if (is.null(out$disp_upper)) out$disp_upper <- win$disp_upper
  }
  if (out$disp_upper <= out$disp_lower) {
    stop(
      "dGamma() measurement prior: implied disp_upper must exceed disp_lower.",
      call. = FALSE
    )
  }
  out
}

#' Build per-group ING Block~1 measurement \code{prior_list} for lmebayes glue
#'
#' Third \code{dispersion_ranef} option: a named list of \code{dGamma()}
#' pfamilies, one per group. Unlike \code{.lmebayes_ing_measurement_prior_list()}
#' (single pooled \code{shape}/\code{rate}), each group keeps its own
#' \code{shape}/\code{rate}/\code{disp_lower}/\code{disp_upper}.
#' @noRd
.lmebayes_ing_measurement_prior_list_group <- function(prior, disp_info, design) {
  re_names <- design$re_coef_names
  p_re     <- length(re_names)
  pl       <- disp_info$dispersion_prior_list
  req <- c("shape_group", "rate_group", "disp_lower_group", "disp_upper_group")
  miss <- req[!req %in% names(pl)]
  if (length(miss)) {
    stop(
      "Internal error: per-group dGamma() dispersion_ranef prior_list must ",
      "contain ", paste(req, collapse = ", "), ".",
      call. = FALSE
    )
  }
  mu <- matrix(
    0,
    nrow = p_re,
    ncol = 1L,
    dimnames = list(re_names, NULL)
  )
  Sigma <- as.matrix(prior$Sigma_ranef)
  if (nrow(Sigma) != p_re || ncol(Sigma) != p_re) {
    stop(
      "prior$Sigma_ranef must be ", p_re, " x ", p_re, ".",
      call. = FALSE
    )
  }
  list(
    mu               = mu,
    Sigma            = Sigma,
    shape_group      = pl$shape_group,
    rate_group       = pl$rate_group,
    disp_lower_group = pl$disp_lower_group,
    disp_upper_group = pl$disp_upper_group,
    max_disp_perc    = if (!is.null(pl$max_disp_perc)) pl$max_disp_perc else 0.99
  )
}

#' @noRd
.lmebayes_add_fixef_summaries <- function(x) {
  if (!is.null(x$fixef)) {
    x$fixef.means <- lapply(x$fixef, colMeans)
  }
  if (!is.null(x$fixef.dispersion)) {
    x$fixef.dispersion.mean <- colMeans(x$fixef.dispersion)
  }
  if (!is.null(x$fixef.iters) && !is.null(x$m_convergence)) {
    x$fixef.iters.mean <- colMeans(x$fixef.iters) / x$m_convergence
  }
  if (!is.null(x$ranef.iters) && !is.null(x$m_convergence)) {
    x$ranef.iters.mean <- mean(x$ranef.iters) / x$m_convergence
  }
  if (!is.null(x$sigma2.iters) && !is.null(x$m_convergence)) {
    x$sigma2.iters.mean <- colMeans(x$sigma2.iters) / x$m_convergence
  }
  x
}

#' @noRd
.lmebayes_block2_icm_labels <- function(prior, family = gaussian()) {
  any_ing <- isTRUE(prior$any_non_normal)
  is_gauss <- is.null(family) || identical(family$family, "gaussian")
  ref_label <- "prior mean"
  if (any_ing) {
    icm_label   <- "gamma @ lmer tau2"
    icm_verbose <- "Block 2 start at lmer tau^2 plug-in"
    conv_label  <- "Plug-in fixed point"
  } else if (is_gauss) {
    icm_label   <- "ICM mean"
    icm_verbose <- "ICM posterior mean"
    conv_label  <- "ICM"
  } else {
    icm_label   <- "ICM mode"
    icm_verbose <- "ICM posterior mode"
    conv_label  <- "ICM"
  }
  list(
    ref_label   = ref_label,
    icm_label   = icm_label,
    icm_verbose = icm_verbose,
    conv_label  = conv_label
  )
}

#' Per-group full column-rank flag for Block~1 \code{Z_j} (same rule as
#' \code{model_setup()$re_rank}).
#' @noRd
.lmebayes_re_rank_from_Z <- function(Z, groups, group_levels = NULL) {
  Z <- as.matrix(Z)
  g_chr <- as.character(groups)
  levs <- if (is.null(group_levels)) {
    unique(g_chr)
  } else {
    as.character(group_levels)
  }
  p_re <- ncol(Z)
  stats::setNames(
    vapply(
      levs,
      function(lev) {
        rows <- which(g_chr == lev)
        Z_j  <- Z[rows, , drop = FALSE]
        nrow(Z_j) >= p_re &&
          Matrix::rankMatrix(Z_j, method = "qr")[1L] == p_re
      },
      logical(1L)
    ),
    levs
  )
}

#' @noRd
.lmebayes_print_icm_fixef_table <- function(
    prior_list,
    re_names,
    fixef_icm,
    icm_info,
    ref_label,
    icm_label,
    conv_label = "ICM",
    header,
    verbose
) {
  if (!isTRUE(verbose) || is.null(fixef_icm)) {
    return(invisible(NULL))
  }
  fixef_ref <- lapply(prior_list, `[[`, "mu_fixef")
  names(fixef_ref) <- re_names
  hdr <- sprintf("  %-18s  %-30s  %14s  %18s",
                 "RE component", "parameter", ref_label, icm_label)
  sep <- paste0("  ", strrep("-", nchar(hdr) - 2L))
  cat(header, "\n")
  cat(hdr, "\n")
  cat(sep, "\n")
  for (k in re_names) {
    nms_k  <- names(fixef_ref[[k]])
    ref_v  <- fixef_ref[[k]]
    icm_v  <- fixef_icm[[k]]
    for (nm in nms_k) {
      cat(sprintf("  %-18s  %-30s  %14.4f  %18.4f\n",
                  k, nm, ref_v[[nm]], icm_v[[nm]]))
    }
  }
  if (!is.null(icm_info)) {
    cat(sprintf("  (%s converged: %s, %d iter, delta = %.2e)\n\n",
                conv_label,
                icm_info$converged, icm_info$iterations, icm_info$delta))
  } else {
    cat("\n")
  }
  invisible(NULL)
}

#' @noRd
.lmebayes_print_ranef_mode_reference <- function(
    ranef_mode,
    re_names,
    group_levels,
    verbose
) {
  invisible(NULL)
}

#' @noRd
.lmebayes_print_fixef_init <- function(
    fixef_init,
    re_names,
    verbose,
    header = "--- main-stage fixef.init (pilot colMeans) ---"
) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }
  cat(header, "\n")
  for (k in re_names) {
    for (nm in names(fixef_init[[k]])) {
      cat(sprintf("  %-18s  %-30s  %12.4f\n",
                  k, nm, fixef_init[[k]][[nm]]))
    }
  }
  cat("\n")
  invisible(NULL)
}
