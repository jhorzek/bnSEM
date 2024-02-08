#' banSEM
#'
#' banSEM (Bayesian network SEM) translates lavaan models to Bayesian networks
#' fitted with bnlearn. The resutling network can be used to investigate conditional
#' distributions of the SEM.
#'
#' @param lavaan_model fitted lavaan model
#' @param phantom_free what to free in the phantom variables. Currently only supports "variance"
#' @returns list with
#' \itemize{
#'  \item{bayes_net: }{A fitted Bayesian network of class bn.fit}
#'  \item{dag: }{The underlying directed acyclical graph}
#'  \item{internal: }{Internal elements}
#' }
#' @import lavaan
#' @import bnlearn
#' @export
#' @examples
#' library(lavaan)
#' library(banSEM)
#' model <- '
#'   # latent variable definitions
#'      ind60 =~ x1 + x2 + x3
#'      dem60 =~ y1 + a*y2 + b*y3 + c*y4
#'      dem65 =~ y5 + a*y6 + b*y7 + c*y8
#'
#'   # regressions
#'     dem60 ~ ind60
#'     dem65 ~ ind60 + dem60
#'
#'   # residual correlations
#'     y1 ~~ y5
#'     y2 ~~ y4 + y6
#'     y3 ~~ y7
#'     y4 ~~ y8
#'     y6 ~~ y8
#' '
#'
#' lavaan_model <- sem(model,
#'                     data = PoliticalDemocracy,
#'                     meanstructure = TRUE)
#'
#' network <- banSEM::banSEM(lavaan_model = lavaan_model)
#'
#' # plot network
#' plot(network$dag)
#'
#' # Check conditional distribution
#' # probability that dem65 in (1,2) given dem60 > 1:
#' bnlearn::cpquery(fitted = network$bayes_net,
#'                  event = (dem65 > 1 & dem65 < 2),
#'                  evidence = (dem60 > 1))
#'
#' # Get distribution under this assumption:
#' dist <- bnlearn::cpdist(fitted = network$bayes_net,
#'                        node = "dem65",
#'                        evidence = (dem60 > 1))
#' hist(dist$dem65)
#'
#' # simulate data from the network and refit SEM to check if the estimates align:
#' sim <- bnlearn::rbn(x = network$bayes_net, n = 10000)
#'
#' fit_sim <- sem(model,
#'                data = sim[,lavaan_model@Data@ov.names[[1]]],
#'                meanstructure = TRUE)
#' round(abs(coef(fit_sim) -
#'             coef(lavaan_model)) / abs(coef(lavaan_model)), 3)
banSEM <- function(lavaan_model,
                   phantom_free = "variance"){

  ##### Setup lavaan model & parameters ####

  check_lavaan_model(lavaan_model = lavaan_model,
                     phantom_free = phantom_free)

  # Extract parameter table
  parameter_table <- lavaan_model@ParTable |>
    as.data.frame() |>
    add_labels()

  # We have to replace all covariances with directed effects of phantom variables
  if(any((parameter_table$op == "~~") &
         (parameter_table$lhs != parameter_table$rhs))){

    lavaan_model_int <- cov_to_phantom(parameter_table,
                                       lavaan_model,
                                       phantom_free)

  }else{

    lavaan_model_int <- lavaan_model

  }

  # extract new parameter table
  parameter_table <- lavaan_model_int@ParTable |>
    as.data.frame()

  # remove equality constraints and simplify
  parameter_table <- parameter_table[!parameter_table$op %in% c("<", ">", "=="),]
  parameter_table <- parameter_table[,c("lhs", "op", "rhs", "label", "est")]

  ##### Translate to Bayesian network ####
  parameter_table_bn <- pt_to_pt_bn(parameter_table = parameter_table)
  bn_model <- create_bn_model(parameter_table_bn = parameter_table_bn)

  # create dag
  dag <- bnlearn::model2network(string = bn_model)
  # set up parameters
  bn_fit = bnlearn::custom.fit(dag,
                               dist = parameter_table_bn)

  return(list(bayes_net = bn_fit,
              dag = dag,
              internal = list(parameter_table_bn = parameter_table_bn,
                              internal_model = lavaan_model_int)))
}
