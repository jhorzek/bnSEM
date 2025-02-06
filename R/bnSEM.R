#' bnSEM
#'
#' bnSEM (Bayesian network SEM) translates OpenMx models to Bayesian networks
#' fitted with bnlearn. The resulting network can be used to investigate conditional
#' distributions of the SEM.
#'
#' @param mx_model fitted OpenMx model of type MxRAMModel
#' @param phantom_free what to free in the phantom variables. Currently only supports "variance"
#' @param optimize should the substitute model in case of covariances be optimized? If
#' the model is not optimized, the Bayesian Network may deviate from the SEM.
#' @returns list with
#'  \item{bayes_net: }{A fitted Bayesian network of class bn.fit}
#'  \item{dag: }{The underlying directed acyclical graph}
#'  \item{internal: }{Internal elements}
#' @import mxsem
#' @import OpenMx
#' @import bnlearn
#' @export
#' @examples
#' library(mxsem)
#' library(bnSEM)
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
#' mx_model <- mxsem(model,
#'                   data = OpenMx::Bollen) |>
#'   mxTryHard()
#'
#' network <- bnSEM::bnSEM(mx_model = mx_model)
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
#'                         node = "dem65",
#'                         evidence = (dem60 > 1))
#' hist(dist$dem65)
#'
#' # simulate data from the network and refit SEM to check if the estimates align:
#' sim <- bnlearn::rbn(x = network$bayes_net, n = 10000)
#'
#' fit_sim <- mxsem(model,
#'                  data = sim[,mx_model$manifestVars]) |>
#'   mxTryHard()
#' round(abs(coef(fit_sim) -
#'             coef(mx_model)) / abs(coef(mx_model)), 3)
bnSEM <- function(mx_model,
                  phantom_free = "variance",
                  optimize = TRUE){

  ##### Setup model & parameters ####

  check_mx_model(mx_model = mx_model,
                 phantom_free = phantom_free,
                 optimize = optimize)

  # Extract parameter table
  parameter_table <- mx_model |>
    OpenMx::omxLocateParameters()

  # We have to replace all covariances with directed effects of phantom variables
  if(any((parameter_table$matrix == "S") &
         (parameter_table$row != parameter_table$col))){

    mx_model_int <- cov_to_phantom(parameter_table,
                                   mx_model,
                                   phantom_free,
                                   optimize)

  }else{

    mx_model_int <- mx_model

  }

  # extract new parameter table
  parameter_table <- mx_model_int |>
    OpenMx::omxLocateParameters()

  # simplify
  parameter_table <- parameter_table[,c("label", "matrix", "row", "col", "value")]

  ##### Translate to Bayesian network ####
  parameter_table_bn <- pt_to_pt_bn(parameter_table = parameter_table,
                                    mx_model_int = mx_model_int)

  bn_model <- create_bn_model(parameter_table_bn = parameter_table_bn)

  # create dag
  dag <- bnlearn::model2network(string = bn_model)
  # set up parameters
  bn_fit = bnlearn::custom.fit(dag,
                               dist = parameter_table_bn)

  # also return an extended data set that could be used with bnlearn
  data_set <- mx_model@data$observed[,mx_model@manifestVars] |>
    as.data.frame()

  # add variables
  for(v in names(dag$nodes)){
    if(!v %in% colnames(data_set))
      data_set[[v]] <- NA_real_
  }

  return(list(bayes_net = bn_fit,
              dag = dag,
              internal = list(parameter_table_bn = parameter_table_bn,
                              internal_model = mx_model_int,
                              extended_data = data_set)))
}
