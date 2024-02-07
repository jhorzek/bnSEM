#' create_bn_model
#'
#' Creates a model for bnlearn from the bn parameter table
#' @param parameter_table_bn parameter table for bn
#' @returns bn model string
#' @keywords internal
create_bn_model <- function(parameter_table_bn){
  # create network
  bn_model <- c()
  for(i in 1:length(parameter_table_bn)){

    if(length(parameter_table_bn[[i]][["coef"]]) == 1){
      # it's an exogenous variable!
      bn_model <- c(bn_model,
                    paste0("[", names(parameter_table_bn[i]), "]"))
    }else{
      preds <- names(parameter_table_bn[[i]][["coef"]])
      preds <- preds[preds != "(Intercept)"]
      bn_model <- c(bn_model,
                    paste0("[", names(parameter_table_bn[i]), "|", paste0(preds, collapse = ":"), "]"))
    }
  }
  bn_model <- paste0(bn_model, collapse = "")

  return(bn_model)
}
