#' pt_to_pt_bn
#'
#' Extracts the regressions underlying the SEM and puts them in a parameter
#' table that bnlearn understands. See https://www.bnlearn.com/examples/custom-fitted/
#' (section "Continuous networks").
#' @param parameter_table parameter table
#' @param mx_model_int mxModel with replaced covariances
#' @returns parameter table for bnlearn
#' @keywords internal
pt_to_pt_bn <- function(parameter_table,
                        mx_model_int){

  variables <- unique(c(mx_model_int$manifestVars, mx_model_int$latentVars))

  # we need to replace all effects with directed effects
  parameter_table_bn <- vector("list", length(variables))
  names(parameter_table_bn) <- variables

  for(v in variables){
    # we want to create a linear model for each of the variables
    coefs <- list()
    # extract intercepts
    coefs[["(Intercept)"]] <- unname(mx_model_int$M$values[,v])

    # get all loadings and regressions for this variable
    vals <- mx_model_int$A$values[v, , drop = FALSE]
    free <- mx_model_int$A$free[v, , drop = FALSE]
    for(eff in which((vals != 0)| free)){
      coefs[[colnames(vals)[eff]]] <- unname(vals[,eff])
    }

    pars <- list(coef = unlist(coefs))

    # extract residual variance
    pars[["sd"]] <- unname(sqrt(mx_model_int$S$values[v,v]))

    parameter_table_bn[[v]] <- pars
  }

  return(parameter_table_bn)
}
