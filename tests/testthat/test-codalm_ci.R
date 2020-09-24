test_that("bootstrap CI works with sequential evaluation", {
    require(gtools)
    require(future)
    set.seed(123)
    x <- rdirichlet(100, rep(1, 2))
    y <- rdirichlet(100, rep(1, 3))
    ci <- codalm_ci(y, x, nboot = 100)
    B_est <- codalm(y, x)
    B_ci_L <- ci$ci_L
    expect_true(is.matrix(B_ci_L))
    expect_true(mean(B_ci_L >= 0) == 1)
    expect_true(all.equal(dim(B_ci_L), c(2,3)))
    expect_true(mean(B_ci_L <= B_est) == 1)

    B_ci_U <- ci$ci_U
    expect_true(is.matrix(B_ci_U))
    expect_true(mean(B_ci_U >= 0) == 1)
    expect_true(all.equal(dim(B_ci_U), c(2,3)))
    expect_true(mean(B_ci_U >= B_est) == 1)
    oplan <- plan("list")
    expect_true(identical(plan("list"), oplan))
})
test_that("bootstrap CI works with multisession evaluation", {
    require(gtools)
    require(future)
    set.seed(123)
    x <- rdirichlet(100, rep(1, 2))
    y <- rdirichlet(100, rep(1, 3))
    ci <- codalm_ci(y, x, nboot = 100, parallel = TRUE, strategy = 'multisession',
                    ncpus = 2)
    B_est <- codalm(y, x)
    B_ci_L <- ci$ci_L
    expect_true(is.matrix(B_ci_L))
    expect_true(mean(B_ci_L >= 0) == 1)
    expect_true(all.equal(dim(B_ci_L), c(2,3)))
    expect_true(mean(B_ci_L <= B_est) == 1)

    B_ci_U <- ci$ci_U
    expect_true(is.matrix(B_ci_U))
    expect_true(mean(B_ci_U >= 0) == 1)
    expect_true(all.equal(dim(B_ci_U), c(2,3)))
    expect_true(mean(B_ci_U >= B_est) == 1)
    oplan <- plan("list")
    expect_true(identical(plan("list"), oplan))
})
