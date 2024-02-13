test_that("MIMIC works", {
  library(mxsem)
  library(lavaan)
  library(banSEM)
  model <- "
  eta =~ y1 + y2 + y3 + y4
  eta ~  x1 + x2 + x3 + x4 + x5
  "

  suppressWarnings(data <- lavaan::simulateData(model))

  fit_mx <- mxsem(model,
                  data = data) |>
    mxTryHard()

  bn <- banSEM::banSEM(mx_model = fit_mx)

  # simulate data from the network and refit SEM to check if the estimates align:
  sim <- bnlearn::rbn(x = bn$bayes_net, n = 100000)

  fit_sim <- mxsem(model,
                 data = sim[,fit_mx$manifestVars]) |>
    mxTryHard()
  testthat::expect_true(all(abs(coef(fit_sim) -
                                  coef(fit_mx)) < .1))
})
