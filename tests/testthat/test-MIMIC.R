test_that("MIMIC works", {
  library(lavaan)
  library(banSEM)
  model <- "
  eta =~ y1 + y2 + y3 + y4
  eta ~  x1 + x2 + x3 + x4 + x5
  "

  suppressWarnings(data <- lavaan::simulateData(model))

  fit_lavaan <- cfa(model,
                    data = data,
                    missing = "ml",
                    fixed.x = FALSE,
                    meanstructure = TRUE)

  bn <- banSEM::banSEM(lavaan_model = fit_lavaan)

  # simulate data from the network and refit SEM to check if the estimates align:
  sim <- bnlearn::rbn(x = bn$bayes_net, n = 100000)

  fit_sim <- sem(model,
                 data = sim[,fit_lavaan@Data@ov.names[[1]]],
                 missing = "ml",
                 fixed.x = FALSE,
                 meanstructure = TRUE)
  testthat::expect_true(all(abs(coef(fit_sim) -
                                  coef(fit_lavaan)) < .1))
})
