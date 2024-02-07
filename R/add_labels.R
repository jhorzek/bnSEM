#' add_labels
#'
#' Adds labels to all unlabeled parameters of the lavaan model
#' @param parameter_table parameter table from lavaan
#' @keywords internal
add_labels <- function(parameter_table){
  # set labels for all parameters
  for(p in 1:nrow(parameter_table)){
    if(parameter_table$label[p] != "")
      next

    parameter_table$label[p] <- paste0(parameter_table$lhs[p],
                                       parameter_table$op[p],
                                       parameter_table$rhs[p])
  }

  return(parameter_table)
}
