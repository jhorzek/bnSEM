#' cov_to_phantom
#'
#' In Bayesian networks, all dependencies between variables must be expressed with
#' directed "effects". There is no covariance in the same sense as in lavaan.
#' Instead, we must replace covariances with effects of a latent phantom variable.
#' For instance x1 ~~ x2 is replaced with x1 =~ 1*ph; x2 =~ 1*ph; ph ~~ ph. The
#' resulting model is identical to the initial lavaan model in terms of fit, but the residual
#' variance estimates will change.
#' @param parameter_table parameter table of a lavaan model
#' @param lavaan_model fitted lavaan model
#' @param phantom_free what to free in the phantom variables. Currently only supports "variance"
#' @returns fitted lavaan_model with phantom variables
#' @importFrom methods is
#' @keywords internal
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
#' lavaan_model <- sem(model,
#'                     data = PoliticalDemocracy,
#'                     meanstructure = TRUE)
#' lavaan_model@ParTable |>
#'   as.data.frame() |>
#'   banSEM:::add_labels() |>
#'   banSEM:::cov_to_phantom(lavaan_model)
cov_to_phantom <- function(parameter_table,
                           lavaan_model,
                           phantom_free = "variance"){

  if(!is(parameter_table, "data.frame"))
    stop("parameter_table must be a data.frame. Use as.data.frame(parameter_table).")

  message("Found covariances in your model. The model will be translated to an",
          " equivalent model with phantom variables. The model with phantom ",
          "variables will be returned in the internal list as 'internal_model'")

  # We constrain the estimates of all parameters we already know and only estimate
  # the variances of variables that also covary
  parameter_table$free[] <- 0

  for(i in which((parameter_table$op == "~~") & (parameter_table$lhs != parameter_table$rhs))){

    cov_at <- parameter_table[i,]
    # we need to add a latent phantom variable
    current_value <- parameter_table$est[i]
    # Add intercept of 0
    parameter_table <- rbind(parameter_table,
                             add_parameter(id = max(parameter_table$id) + 1,
                                           lhs = paste0("ph_", parameter_table$label[i]),
                                           op = "~1",
                                           rhs = "",
                                           free = 0,
                                           label = paste0(cov_at$label, "~1"),
                                           plabel = paste0(".p",nrow(parameter_table) + 1,"."),
                                           start = 0,
                                           est = 0,
                             ))

    # Add a freely estimated variance
    parameter_table <- rbind(parameter_table,
                             add_parameter(id = max(parameter_table$id) + 1,
                                           lhs = paste0("ph_", parameter_table$label[i]),
                                           op = "~~",
                                           rhs = paste0("ph_", parameter_table$label[i]),
                                           free = ifelse(phantom_free == "variance", max(parameter_table$free) + 1, 0),
                                           label = parameter_table$label[i],
                                           plabel = paste0(".p",nrow(parameter_table) + 1,"."),
                                           start = ifelse(phantom_free == "variance", abs(current_value), 1),
                                           est = ifelse(phantom_free == "variance", abs(current_value), 1),
                             ))

    # Additionally, we need to specify loadings of 1 / -1 on the covarying items
    # Note: For positive covariances the loadings are all 1. For negative covariances,
    # the loadings on one item are 1 and on the other -1.
    parameter_table <- rbind(parameter_table,
                             add_parameter(id = max(parameter_table$id) + 1,
                                           lhs = paste0("ph_", parameter_table$label[i]),
                                           op = "=~",
                                           rhs = parameter_table$lhs[i],
                                           free = 0,
                                           label = paste0(parameter_table$lhs[i], "=~",
                                                          paste0("ph_", parameter_table$label[i])),
                                           plabel = paste0(".p",nrow(parameter_table) + 1,"."),
                                           start = sign(current_value)*1,
                                           est = sign(current_value)*1,
                             ))


    parameter_table <- rbind(parameter_table,
                             add_parameter(id = max(parameter_table$id) + 1,
                                           lhs = paste0("ph_", parameter_table$label[i]),
                                           op = "=~",
                                           rhs = parameter_table$rhs[i],
                                           free = 0,
                                           label = paste0(parameter_table$rhs[i], "=~",
                                                          paste0("ph_", parameter_table$label[i])),
                                           plabel = paste0(".p",nrow(parameter_table) + 1,"."),
                                           start = 1,
                                           est = 1,
                             ))

    # set variances free
    parameter_table$free[(parameter_table$rhs == parameter_table$rhs[i]) &
                           (parameter_table$lhs == parameter_table$rhs[i]) &
                           (parameter_table$op == "~~")] <- max(parameter_table$free) + 1

    parameter_table$free[(parameter_table$rhs == parameter_table$lhs[i]) &
                           (parameter_table$lhs == parameter_table$lhs[i]) &
                           (parameter_table$op == "~~")] <- max(parameter_table$free) + 1

  }

  # Remove all covariances. These are now captured by the phantom variables.

  parameter_table <- parameter_table[!((parameter_table$op == "~~") &
                                         (parameter_table$lhs != parameter_table$rhs)),]


  # in case of phantom variables, we now have to refit the model
  message("Refitting the model with phantom variables. The fit will ",
          "be the same, but the variances of the residuals will change.")
  call_lavaan <- list(model = parameter_table,
                      data = lavaan_model@Data)
  call_lavaan <- c(call_lavaan,
                   lavaan_model@Options[names(lavOptions())])

  lavaan_model_int <- do.call(what = lavaan::lavaan,
                              args = call_lavaan)

  # check fit
  if(abs(logLik(lavaan_model_int) - logLik(lavaan_model)) / logLik(lavaan_model)  > .01)
    stop("Refactoring the model failed")


  return(lavaan_model_int)

}
