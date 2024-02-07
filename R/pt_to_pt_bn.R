#' pt_to_pt_bn
#'
#' Extracts the regressions underlying the SEM and puts them in a parameter
#' table that bnlearn understands. See https://www.bnlearn.com/examples/custom-fitted/
#' (section "Continuous networks").
#' @param parameter_table parameter table
#' @returns parameter table for bnlearn
#' @keywords internal
pt_to_pt_bn <- function(parameter_table){

  # we need to replace all effects with directed effects
  variables <- unique(unlist(parameter_table[,c("lhs", "rhs")]))
  variables <- variables[variables != ""]

  parameter_table_bn <- vector("list", length(variables))
  names(parameter_table_bn) <- variables

  for(v in variables){
    # we want to create a linear model for each of the variables
    coefs <- list()
    # extract intercepts
    if(any((parameter_table$lhs == v) & (parameter_table$op == "~1"))){
      coefs[["(Intercept)"]] <- parameter_table$est[parameter_table$lhs == v & parameter_table$op == "~1"]
    }else{
      coefs[["(Intercept)"]] <- 0.0
    }

    # get all loadings and regressions for this variable
    for(reg in which((parameter_table$lhs == v) & (parameter_table$op == "~"))){
      coefs[[parameter_table$rhs[reg]]] <- parameter_table$est[reg]
    }
    for(lam in which((parameter_table$rhs == v) & (parameter_table$op == "=~"))){
      coefs[[parameter_table$lhs[lam]]] <- parameter_table$est[lam]
    }
    pars <- list(coef = unlist(coefs))

    # extract residual variance
    if(any((parameter_table$lhs == v) & (parameter_table$op == "~~") & (parameter_table$rhs == v))){
      pars[["sd"]] <- sqrt(parameter_table$est[(parameter_table$lhs == v) &
                                                 (parameter_table$op == "~~") &
                                                 (parameter_table$rhs == v)])
    }else{
      pars[["sd"]] <- 0.0
    }
    parameter_table_bn[[v]] <- pars
  }

  return(parameter_table_bn)
}
