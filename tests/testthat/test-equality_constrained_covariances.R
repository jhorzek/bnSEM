test_that("Equality Constrained Covariances", {
  library(mxsem)
  library(bnSEM)
  set.seed(45098)
  model <- '
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ f*y5
    y2 ~~ f*y4 + f*y6
    y3 ~~ f*y7
    y4 ~~ f*y8
    y6 ~~ f*y8
'

  mx_model <- mxsem(model,
                    data = OpenMx::Bollen) |>
    mxTryHard()

  network <- bnSEM::bnSEM(mx_model = mx_model)

  # simulate data from the network and refit SEM to check if the estimates align:
  sim <- bnlearn::rbn(x = network$bayes_net, n = 10000)

  fit_sim <- mxsem(model,
                   data = sim[,mx_model$manifestVars]) |>
    mxTryHard()

  testthat::expect_true(all(abs(coef(fit_sim) -
                                  coef(mx_model)) / abs(coef(mx_model)) < .1))

})
