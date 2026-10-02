test_that("normalize_group covers factor, integer, and list partitions", {
  l2 <- 30L
  school <- rep(1:3, each = 10L)
  blk_factor <- normalize_group(factor(school), l2)
  blk_int    <- normalize_group(as.integer(school), l2)
  blk_list   <- normalize_group(split(seq_len(l2), school), l2)
  expect_identical(blk_factor$k, blk_int$k)
  expect_identical(blk_factor$l2_blocks, blk_int$l2_blocks)
  expect_equal(blk_factor$k, 3L)
  expect_equal(sum(blk_list$l2_blocks), l2)
})

test_that("rNormal_reg_group returns k x l1 coefficients and finite coef.mode", {
  set.seed(42)
  n_groups <- 3L
  n_per    <- 10L
  school   <- rep(seq_len(n_groups), each = n_per)
  x        <- cbind(1, rnorm(n_groups * n_per))
  l1       <- ncol(x)
  sigma2   <- 1.5
  y        <- rnorm(nrow(x), sd = sqrt(sigma2))
  prior_list <- list(
    mu         = rep(0, l1),
    Sigma      = diag(100, l1),
    dispersion = sigma2,
    ddef       = FALSE
  )

  out <- rNormal_reg_group(
    n = 1L, y = y, x = x, group = school, prior_list = prior_list
  )
  expect_s3_class(out, "rNormal_reg_group")
  expect_equal(dim(out$coefficients), c(n_groups, l1))
  expect_true(all(is.finite(out$coef.mode)))
  expect_equal(out$group_info$k, n_groups)
})

test_that("pfamily_list.Prior_SetupGroup builds dNormal pfamilies per block", {
  data("iris", package = "datasets")
  ps <- Prior_SetupGroup(
    Sepal.Length ~ Sepal.Width,
    group = "Species",
    data = iris,
    family = gaussian()
  )
  pf <- pfamily_list(ps)
  expect_s3_class(pf, "pfamily_list")
  expect_equal(length(pf), length(ps))
  expect_identical(names(pf), names(ps))
  expect_true(all(vapply(pf, function(p) inherits(p, "pfamily"), logical(1L))))
})

test_that("pfamily_list.Prior_SetupGroup matches lmebayesCore Ex example (3 covariates)", {
  data("iris", package = "datasets")
  ps_block <- Prior_SetupGroup(
    Sepal.Length ~ Sepal.Width + Petal.Length,
    group = "Species",
    data = iris,
    family = gaussian()
  )
  pf <- suppressMessages(pfamily_list(ps_block))
  expect_equal(names(pf), names(ps_block))
  pf_ng <- suppressMessages(pfamily_list(ps_block, ptypes = "dNormal_Gamma"))
  expect_equal(attr(pf_ng, "ptypes"), stats::setNames(rep("dNormal_Gamma", 3L), names(ps_block)))
})
