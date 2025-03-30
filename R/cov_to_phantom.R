#' cov_to_phantom
#'
#' In Bayesian networks, all dependencies between variables must be expressed with
#' directed "effects". There is no covariance in the same sense as in OpenMx.
#' Instead, we must replace covariances with effects of a latent phantom variable.
#' For instance x1 <-> x2 is replaced with ph -> y1; ph -> y2; ph <-> ph. The
#' resulting model is identical to the initial OpenMx model in terms of fit, but the residual
#' variance estimates will change. The implementation has been improved by Dr.
#' Christian Gische
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
                           phantom_free,
                           phantom_variance_start,
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

    # get variances; these will be used for scaling the loadings of
    # the phantom variables. The approach used in the
    # following is based on p. 8 in Merkle, E. C., & Rosseel,
    # Y. (2015). blavaan: Bayesian structural equation models via
    # parameter expansion. arXiv preprint arXiv:1511.05604.
    var_ii <- mx_model$S$values[cov_at$row, cov_at$row]
    var_jj <- mx_model$S$values[cov_at$col, cov_at$col]

    # we need to add a latent phantom variable
    current_value <- parameter_table$value[i]
    # Add intercept of 0
    ph_var <- 1
    while(paste0("ph_", parameter_table$label[i], "_", ph_var) %in% mx_model_int$latentVars){
      ph_var <- ph_var + 1
    }
    new_latent <- paste0("ph_", parameter_table$label[i], "_", ph_var)

    mx_model_int <- mxModel(mx_model_int,
                            latentVars = new_latent,
                            # Add intercept of 0
                            mxPath(from = "one",
                                   to = new_latent,
                                   values = 0,
                                   free = FALSE,
                                   labels = paste0("one", mxsem::unicode_directed(),
                                                   new_latent)),
                            # Add variance of 1
                            mxPath(from = new_latent,
                                   to = new_latent,
                                   arrows = 2,
                                   values = phantom_variance_start,
                                   free = TRUE,
                                   lbound = 1e-6,
                                   labels = parameter_table$label[i]))

    # Additionally, we need to specify loadings on the covarying items
    # One loading will be constrained to 1, the other freely estimated
    mx_model_int$A$values[parameter_table$row[i],new_latent] <- sqrt(abs(current_value) * var_ii)
    mx_model_int$A$free[parameter_table$row[i],new_latent] <- FALSE

    mx_model_int$A$values[parameter_table$col[i],new_latent] <- sign(current_value)*sqrt(abs(current_value) * var_ii)
    mx_model_int$A$free[parameter_table$col[i],new_latent] <- FALSE

    # remove covariance
    mx_model_int$S$values[parameter_table$row[i],parameter_table$col[i]] <- 0
    mx_model_int$S$free[parameter_table$row[i],parameter_table$col[i]] <- FALSE
    mx_model_int$S$labels[parameter_table$row[i],parameter_table$col[i]] <- ""

    mx_model_int$S$values[parameter_table$col[i],parameter_table$row[i]] <- 0
    mx_model_int$S$free[parameter_table$col[i],parameter_table$row[i]] <- FALSE
    mx_model_int$S$labels[parameter_table$col[i],parameter_table$row[i]] <- ""
  }

  if(optimize){
    # in case of phantom variables, we now have to refit the model
    message("Refitting the model with phantom variables. The fit will ",
            "be the same, but the variances of the residuals will change.")

    mx_model_int <- mx_model_int |>
      mxTryHard()

    # check fit
    if(abs(logLik(mx_model_int) - logLik(mx_model)) / abs(logLik(mx_model))  > .01)
      warning("Refactoring the model failed. Observed deviations in likelihood!")
  }else{
    warning("The model with phantom variables has not been optimized")
  }
  return(mx_model_int)
}
