#' bnSEM
#'
#' bnSEM (Bayesian network SEM) translates OpenMx models to Bayesian networks
#' fitted with bnlearn. The resulting network can be used to investigate conditional
#' distributions of the SEM.
#'
#' In Bayesian networks, all dependencies between variables must be expressed with
#' directed "effects". There is no covariance in the same sense as in OpenMx.
#' Instead, we must replace covariances with effects of a latent phantom variable.
#' For instance x1 <-> x2 is replaced with ph -> y1; ph -> y2; ph <-> ph. The
#' resulting model is identical to the initial OpenMx model in terms of fit, but the residual
#' variance estimates will change.
#'
#' bnSEM currently supports three different approaches:
#'
#' First, when using phantom_type = "refit-free-loadings", bnSEM will add phantom
#' variables only for covariances and re-estimate the SEM. For "refit-free-loadings",
#' the loadings of the phantom on the observed variables is (1) constraint to
#' equality if the covariance is positive and (2) constraint to have equal values,
#' but opposite signs for covariances that are negative. The variance is constraint
#' to one. This is based on a discussion with Prof. Marcel Paulssen for a bifactor
#' modeling approach.
#'
#' Second, when using phantom_type = "refit-free-single-loading", bnSEM will also add phantom
#' variables only for covariances and re-estimate the SEM. However, in contrast
#' to "refit-free-loadings" only one of the loadings is estimated, while the
#' other one is constraint. This results in exactly the same fit as the first approach
#' and is based on a suggestion by Dr. Christian Gische.
#'
#' The advantage of the refitting approaches is that each manifest and latent
#' variable of the original model will still have a residual
#' variance in the new model and the Bayesian Network. This is necessary for cpdist
#' to work as expected. The main disadvantages are that (1) refitting may fail and (2)
#' constraints on the variance and covariance parameters cannot be accounted for.
#'
#' Third, when using phantom_type = "cholesky", bnSEM will replace the full residual
#' (co-)variance matrix with a Cholesky decomposition (S = DD^t). The latent and
#' manifest variables (v) of the original model are now given by v = Du, where u
#' is a vector of the same size as v with standard-normally distributed items.
#' The main advantage of the Cholesky decomposition approach is that it does not
#' require refitting the model and therefore also does not result in non-convergence.
#' The main disadvantage is that each variable in the original model is now determined
#' fully by the new variables in u. That is, none of the original variables has any
#' residual variances in the new model. This can result in unexpected behavior
#' when using cpdist on the model, where distributions for deterministic variables
#' do not work. This approach has been suggested by Prof. Manuel C. Voelkle.
#'
#' @param mx_model fitted OpenMx model of type MxRAMModel
#' @param phantom_type type of phantom variable approach to use, There are two approaches:
#' "refit-free-loadings" and "cholesky". See details.
#' @param optimize should the substitute model in case of phantom_type = "refit-free-loadings" be optimized?
#' @returns list with
#' \itemize{
#'  \item bayes_net: A fitted Bayesian network of class bn.fit
#'  \item dag: The underlying directed acyclical graph
#'  \item internal: Internal elements
#' }
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
                  phantom_type = "refit-free-loadings",
                  optimize = TRUE){

  ##### Setup model & parameters ####
  check_mx_model(mx_model = mx_model,
                 optimize = optimize)

  # Extract parameter table
  parameter_table <- mx_model |>
    OpenMx::omxLocateParameters()

  # We have to replace all covariances with directed effects of phantom variables
  if(any((parameter_table$matrix == "S") &
         (parameter_table$row != parameter_table$col))){

    mx_model_int <- cov_to_phantom(parameter_table,
                                   mx_model,
                                   phantom_type = phantom_type,
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
