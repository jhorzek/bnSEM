test_that("MIMIC works", {
  library(mxsem)
  library(lavaan)
  library(bnSEM)
  set.seed(73)
  sim_model <- "
  eta =~ 1*y1 + 1*y2 + 1*y3 + 1*y4
  eta ~ .3*x1 + .5*x2 + .2*x3 + .1*x4 + .1*x5
  x1 ~~ -.1*x2 + .1*x3 + .2*x4 + .1*x5
  x2 ~~ .1*x3 + -.1*x4 + .1*x5
  x3 ~~ -.2*x4 + .1*x5
  x4 ~~ .1*x5
  "

  suppressWarnings(data <- lavaan::simulateData(sim_model))

  model <- "
  eta =~ y1 + y2 + y3 + y4
  eta ~  x1 + x2 + x3 + x4 + x5
  "
  fit_mx <- mxsem(model,
                  data = data) |>
    mxTryHard()

  bn <- bnSEM::bnSEM(mx_model = fit_mx)

  # simulate data from the network and refit SEM to check if the estimates align:
  sim <- bnlearn::rbn(x = bn$bayes_net, n = 100000)

  fit_sim <- mxsem(model,
                 data = sim[,fit_mx$manifestVars]) |>
    mxTryHard()
  testthat::expect_true(all(abs(coef(fit_sim) -
                                  coef(fit_mx)) < .1))
})
