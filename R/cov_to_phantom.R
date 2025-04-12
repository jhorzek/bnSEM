#' cov_to_phantom
#'
#' In Bayesian networks, all dependencies between variables must be expressed with
#' directed "effects". There is no covariance in the same sense as in OpenMx.
#' Instead, we must replace covariances with effects of a latent phantom variable.
#' For instance x1 <-> x2 is replaced with ph -> y1; ph -> y2; ph <-> ph. The
#' resulting model is identical to the initial OpenMx model in terms of fit, but the residual
#' variance estimates will change.
#'
#' bnSEM currently supports two different approaches:
#'
#' First, when using phantom_type = "refit", bnSEM will add phantom variables only
#' for covariances and re-estimate the SEM. The advantage of this approach is that
#' each manifest and latent variable of the original model will still have a residual
#' variance in the new model and the Bayesian Network. This is necessary for cpdist
#' to work as expected. The main disadvantages are that (1) refitting may fail and (2)
#' constraints on the variance and covariance parameters cannot be accounted for.
#' The implementation has been improved by Dr. Christian Gische.
#'
#' Second, when using phantom_type = "cholesky", bnSEM will replace the full residual
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
#' @param parameter_table parameter table of an OpenMx model (see ?OpenMx::omxLocateParameters)
#' @param mx_model fitted OpenMx model
#' @param phantom_type type of phantom variable approach to use, There are two approaches:
#' "refit" and "cholesky". See details.
#' @param optimize should the substitute model in case of covariances be optimized?
#' @returns fitted mx_model with phantom variables
#' @importFrom methods is
#' @importFrom stats logLik
#' @keywords internal
cov_to_phantom <- function(parameter_table,
                           mx_model,
                           phantom_type,
                           optimize){

  if(tolower(phantom_type) == "refit-free-loadings"){
    return(cov_to_phantom_refit_loadings(parameter_table = parameter_table,
                                         mx_model = mx_model,
                                         optimize = optimize))
  }else if(tolower(phantom_type) == "refit-free-single-loading"){
    return(cov_to_phantom_refit_single_loading(parameter_table = parameter_table,
                                               mx_model = mx_model,
                                               optimize = optimize))
  }else if(tolower(phantom_type) == "cholesky"){
    return(cov_to_phantom_cholesky(parameter_table = parameter_table,
                                   mx_model = mx_model))
  }else{
    stop("Unkown phantom_type in cov_to_phantom: Please select one of 'refit' or 'cholesky'.")
  }

}


#' cov_to_phantom_refit
#'
#' In Bayesian networks, all dependencies between variables must be expressed with
#' directed "effects". There is no covariance in the same sense as in OpenMx.
#' Instead, we must replace covariances with effects of a latent phantom variable.
#' For instance x1 <-> x2 is replaced with ph -> y1; ph -> y2; ph <-> ph. The
#' resulting model is identical to the initial OpenMx model in terms of fit, but the residual
#' variance estimates will change. The implementation has been improved by Dr.
#' Christian Gische.
#'
#' @param parameter_table parameter table of an OpenMx model (see ?OpenMx::omxLocateParameters)
#' @param mx_model fitted OpenMx model
#' @param optimize should the substitute model in case of covariances be optimized?
#' @returns fitted mx_model with phantom variables
#' @importFrom methods is
#' @importFrom stats logLik
#' @keywords internal
cov_to_phantom_refit_single_loading <- function(parameter_table,
                                                mx_model,
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
                            # Add variance. The main challenge here is that
                            # setting the variance to a value that is too large
                            # often results in non-convergence. Similarly,
                            # if it is too small, convergence also won't happen.
                            # As a simple solution for now, we take the covariance
                            # value as an indicator of how large the variance should
                            # be.
                            mxPath(from = new_latent,
                                   to = new_latent,
                                   arrows = 2,
                                   values = max(c(abs(current_value), .1)),
                                   free = FALSE,
                                   lbound = 1e-6,
                                   labels = paste0(new_latent, mxsem::unicode_undirected(),
                                                   new_latent)))

    # Additionally, we need to specify loadings on the covarying items
    # One loading will be constrained to 1, the other freely estimated
    mx_model_int$A$values[parameter_table$row[i],new_latent] <- 1
    mx_model_int$A$free[parameter_table$row[i],new_latent] <- FALSE

    mx_model_int$A$values[parameter_table$col[i],new_latent] <- 1
    mx_model_int$A$free[parameter_table$col[i],new_latent] <- TRUE

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

#' cov_to_phantom_refit_loadings
#'
#' In Bayesian networks, all dependencies between variables must be expressed with
#' directed "effects". There is no covariance in the same sense as in OpenMx.
#' Instead, we must replace covariances with effects of a latent phantom variable.
#' For instance x1 <-> x2 is replaced with ph -> y1; ph -> y2; ph <-> ph. The
#' resulting model is identical to the initial OpenMx model in terms of fit, but the residual
#' variance estimates will change. For the current implementation, the loadings of the
#' phantom on the observed variables is (1) constraint to equality if the covariance
#' is positive and (2) constraint to have equal values, but opposite signs for covariances
#' that are negative. The variance is constraint to one. This is based on a discussion
#' with Prof. Marcel Paulssen for a bifactor modeling approach.
#'
#' @param parameter_table parameter table of an OpenMx model (see ?OpenMx::omxLocateParameters)
#' @param mx_model fitted OpenMx model
#' @param optimize should the substitute model in case of covariances be optimized?
#' @returns fitted mx_model with phantom variables
#' @importFrom methods is
#' @importFrom stats logLik
#' @keywords internal
cov_to_phantom_refit_loadings <- function(parameter_table,
                                          mx_model,
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

    current_value <- parameter_table$value[i]
    from_label <- colnames(mx_model$S$values)[parameter_table$col[i]]
    to_label   <- rownames(mx_model$S$values)[parameter_table$row[i]]

    # we need to add a latent phantom variable
    ph_var <- 1
    while(paste0("ph_", parameter_table$label[i], "_", ph_var) %in% mx_model_int$latentVars){
      ph_var <- ph_var + 1
    }
    new_latent <- paste0("ph_", parameter_table$label[i], "_", ph_var)
    new_latent_no_unicode <- gsub(x = new_latent,
                                  pattern = mxsem::unicode_undirected(),
                                  replacement = "__")

    mx_model_int <- mxModel(mx_model_int,
                            latentVars = new_latent,
                            # Add intercept of 0
                            mxPath(from = "one",
                                   to = new_latent,
                                   values = 0,
                                   free = FALSE,
                                   labels = paste0("one",
                                                   mxsem::unicode_directed(),
                                                   new_latent)),
                            # Add variance. The main challenge here is that
                            # setting the variance to a value that is too large
                            # often results in non-convergence. Similarly,
                            # if it is too small, convergence also won't happen.
                            # As a simple solution for now, we take the covariance
                            # value as an indicator of how large the variance should
                            # be.
                            mxPath(from = new_latent,
                                   to = new_latent,
                                   arrows = 2,
                                   values = 1,
                                   free = FALSE,
                                   lbound = 1e-6,
                                   labels = paste0(new_latent,
                                                   mxsem::unicode_undirected(),
                                                   new_latent)))

    # Additionally, we need to specify loadings on the covarying items
    mx_model_int$A$values[parameter_table$row[i],new_latent] <- 0
    mx_model_int$A$labels[parameter_table$row[i],new_latent] <- paste0(
      paste0(new_latent_no_unicode, "_", to_label), "[1,1]")
    mx_model_int$A$free[parameter_table$row[i],new_latent] <- FALSE

    mx_model_int$A$values[parameter_table$col[i],new_latent] <- 0
    mx_model_int$A$labels[parameter_table$col[i],new_latent] <- paste0(new_latent_no_unicode,
                                                                       "_",
                                                                       from_label)
    mx_model_int$A$free[parameter_table$col[i],new_latent] <- TRUE

    # Add equality constraint on the loadings
    mx_model_int <- OpenMx::mxModel(mx_model_int,
                                    OpenMx::mxAlgebraFromString(algString = ifelse(
                                      current_value < 0,
                                      paste0("-", mx_model_int$A$labels[parameter_table$col[i],new_latent]),
                                      mx_model_int$A$labels[parameter_table$col[i],new_latent]
                                    ),
                                    name = paste0(new_latent_no_unicode, "_", to_label)))

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


#' cov_to_phantom_cholesky
#'
#' In Bayesian networks, all dependencies between variables must be expressed with
#' directed "effects". There is no covariance in the same sense as in OpenMx.
#' Instead, we must replace covariances with effects of a latent phantom variable.
#'
#' cov_to_phantom_cholesky uses a Cholesky decomposition to represent the covariance
#' matrix of the entire model. This approach has been suggested by Manuel C. Voelkle.
#' See, for example, p. 462 in
#' Klein, A., & Moosbrugger, H. (2000). Maximum likelihood estimation of latent
#' interaction effects with the LMS method. Psychometrika, 65(4), 457-474.
#'
#' The main advantage of using a Cholesky decomposition is that no refitting of the
#' model is required. As a result, many convergence issues are avoided.
#'
#' The main disadvantage of the Cholesky decomposition approach is that we are not
#' just "outsourcing" the covariances to new phantom variables, but also the variances.
#' As a result, all observed and latent variables in the original model are now
#' fully determined by new latent variables. This results in issues with bnlearn,
#' where conditional distributions cannot be computed.
#'
#' @param parameter_table parameter table of an OpenMx model (see ?OpenMx::omxLocateParameters)
#' @param mx_model fitted OpenMx model
#' @returns fitted mx_model with phantom variables
#' @importFrom methods is
#' @importFrom stats logLik
#' @importFrom Matrix nearPD
#' @keywords internal
cov_to_phantom_cholesky <- function(parameter_table,
                                    mx_model){

  mx_model_int <- mx_model

  if(!is(parameter_table, "data.frame"))
    stop("parameter_table must be a data.frame. Use as.data.frame(parameter_table).")

  message("Found covariances in your model. The model will be translated to an",
          " equivalent model with phantom variables. The model with phantom ",
          "variables will be returned in the internal list as 'internal_model'")

  # Cholesky Decomposition - Suggested by Manuel C. Voelkle
  # The transposition provides a new matrix - chol_cov - that
  # reproduces the original covariance matrix with chol_cov %*% t(chol_cov)
  chol_cov <- try(t(chol(mx_model$S$values)), silent = TRUE)
  if(is(chol_cov, "try-error")){
    warning(paste0("The covariance matrix of the SEM is not positive definite, ",
                   "resulting in the Cholesky decomposition to fail. ",
                   "bnSEM will approximate with the nearest positive definite matrix."))
    chol_cov <- t(chol(Matrix::nearPD(x = mx_model$S$values)$mat))
  }
  # When multiplied with a vector of standard-normally distributed variables,
  # we end up with variables that have the same covariance as the initial items
  # items = chol_cov %*% x
  # Based on this, we can compute the variances and covariances of our new
  # chol_cov %*% t(chol_cov)
  # However, we will use another step first: Currently, we will end up with
  # phantom variables that are somewhat arbitrary as they have a variance of
  # 1 and are then connected to the observed variables and latent variables
  # with a directed effect. We will change this such that each phantom variable
  # is connected to one of the manifest / latent variables with an effect of 1
  # (making the phantom item almost identical to the observed one - with the
  # exception of the causal implications of the directed effects of the phantom
  # variables).
  ph_variances <- diag(chol_cov)^2
  names(ph_variances) <- paste0("ph_", rownames(mx_model$S$values))
  ph_effects   <- chol_cov %*% solve(diag(diag(chol_cov)))
  dimnames(ph_effects) <- list(
    rownames(mx_model$S$values),
    paste0("ph_", rownames(mx_model$S$values))
  )

  # First, we remove all variances and covariances from our model
  mx_model_int$S$labels[] <- NA
  mx_model_int$S$values[] <- 0
  mx_model_int$S$free[] <- FALSE

  # We add as many phantom variables as there are columns in our Cholesky decomposition
  for(i in colnames(ph_effects))
    mx_model_int <- mxModel(mx_model_int,
                            latentVars = i,
                            # Add intercept of 0
                            mxPath(from = "one",
                                   to = i,
                                   values = 0,
                                   free = FALSE,
                                   labels = paste0("one", mxsem::unicode_directed(),
                                                   i)),
                            mxPath(from = i,
                                   to = i,
                                   arrows = 2,
                                   values = ph_variances[i],
                                   free = FALSE,
                                   lbound = 1e-6,
                                   labels = paste0("var_", i)))

  # Next, we define directed effects of each of those phantom variables on our
  # model variables based on the Cholesky decomposition
  for(i in rownames(ph_effects)){
    for(j in colnames(ph_effects)){
      value <- ph_effects[i, j]

      mx_model_int <- mxModel(mx_model_int,
                              # Add effect from phantom to observed
                              mxPath(from = j,
                                     to = i,
                                     values = value,
                                     free = TRUE,
                                     labels = paste0(i,
                                                     mxsem::unicode_directed(),
                                                     j)))
    }
  }

  browser()
  mx_model_int <- OpenMx::mxRun(mx_model_int,
                                useOptimizer = FALSE,
                                silent = TRUE)
  logLik(mx_model)
  logLik(mx_model_int)

  if(abs(logLik(mx_model_int) - logLik(mx_model)) / abs(logLik(mx_model))  > .01){
    warning("Refactoring the model failed. Observed deviations in likelihood!")
  }

  return(mx_model_int)
}
