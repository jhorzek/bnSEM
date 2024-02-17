#' cov_to_phantom
#'
#' In Bayesian networks, all dependencies between variables must be expressed with
#' directed "effects". There is no covariance in the same sense as in OpenMx.
#' Instead, we must replace covariances with effects of a latent phantom variable.
#' For instance x1 <-> x2 is replaced with ph -> y1; ph -> y2; ph <-> ph. The
#' resulting model is identical to the initial OpenMx model in terms of fit, but the residual
#' variance estimates will change.
#' @param parameter_table parameter table of an OpenMx model (see ?OpenMx::omxLocateParameters)
#' @param mx_model fitted OpenMx model
#' @param phantom_free what to free in the phantom variables. Currently only supports "variance"
#' @param optimize should the substitute model in case of covariances be optimized?
#' @returns fitted mx_model with phantom variables
#' @importFrom methods is
#' @importFrom stats logLik
#' @keywords internal
cov_to_phantom <- function(parameter_table,
                           mx_model,
                           phantom_free = "variance",
                           optimize){

  mx_model_int <- mx_model

  if(!is(parameter_table, "data.frame"))
    stop("parameter_table must be a data.frame. Use as.data.frame(parameter_table).")

  message("Found covariances in your model. The model will be translated to an",
          " equivalent model with phantom variables. The model with phantom ",
          "variables will be returned in the internal list as 'internal_model'")

  for(i in which((parameter_table$matrix == "S") & (parameter_table$row != parameter_table$col))){

    # covariance matrix is symmetric -> skip if element on opposite side was already
    # visited
    if(which((parameter_table$matrix == "S") &
             (parameter_table$row == parameter_table$col[i]) &
             (parameter_table$col == parameter_table$row[i])) < i)
      next

    cov_at <- parameter_table[i,]
    # we need to add a latent phantom variable
    current_value <- parameter_table$value[i]
    # Add intercept of 0
    new_latent <- ifelse(paste0("ph_", parameter_table$label[i]) %in% mx_model_int$latentVars,
                         paste0("ph_", paste0(sample(x = LETTERS, size = 3), collapse = ""),
                                "_", parameter_table$label[i]),
                         paste0("ph_", parameter_table$label[i]))
    mx_model_int <- mxModel(mx_model_int,
                            latentVars = new_latent,
                            # Add intercept of 0
                            mxPath(from = "one",
                                   to = new_latent,
                                   values = 0,
                                   free = FALSE,
                                   labels = paste0("one", mxsem::unicode_directed(),
                                                   new_latent)),
                            # Add a freely estimated variance
                            mxPath(from = new_latent,
                                   to = new_latent,
                                   arrows = 2,
                                   values = parameter_table$value[i],
                                   free = TRUE,
                                   lbound = 1e-6,
                                   labels = parameter_table$label[i]))

    # Additionally, we need to specify loadings of 1 / -1 on the covarying items
    # Note: For positive covariances the loadings are all 1. For negative covariances,
    # the loadings on one item are 1 and on the other -1.
    mx_model_int$A[parameter_table$row[i],new_latent]$values[] <- sign(current_value)*1
    mx_model_int$A[parameter_table$row[i],new_latent]$free[] <- FALSE

    mx_model_int$A[parameter_table$col[i],new_latent]$values[] <- 1
    mx_model_int$A[parameter_table$col[i],new_latent]$free[] <- FALSE

    # remove covariance
    mx_model_int$S[parameter_table$row[i],parameter_table$col[i]]$values[] <- 0
    mx_model_int$S[parameter_table$row[i],parameter_table$col[i]]$free[] <- FALSE
    mx_model_int$S[parameter_table$row[i],parameter_table$col[i]]$labels[] <- ""

    mx_model_int$S[parameter_table$col[i],parameter_table$row[i]]$values[] <- 0
    mx_model_int$S[parameter_table$col[i],parameter_table$row[i]]$free[] <- FALSE
    mx_model_int$S[parameter_table$col[i],parameter_table$row[i]]$labels[] <- ""
  }

  if(optimize){
    # in case of phantom variables, we now have to refit the model
    message("Refitting the model with phantom variables. The fit will ",
            "be the same, but the variances of the residuals will change.")

    mx_model_int <- mx_model_int |>
      mxTryHard()

    # check fit
    if(abs(logLik(mx_model_int) - logLik(mx_model)) / logLik(mx_model)  > .01)
      warning("Refactoring the model failed. Observed deviations in likelihood!")
  }else{
    warning("The model with phantom variables has not been optimized")
  }
  return(mx_model_int)
}
