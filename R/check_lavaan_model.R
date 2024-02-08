#' check_lavaan_model
#'
#' Check if the lavaan model can be translated to bnlearn
#' @param lavaan_model fitted lavaan model
#' @param phantom_free what to free in the phantom variables. Currently only supports "variance"
#' @keywords internal
check_lavaan_model <- function(lavaan_model,
                               phantom_free){

  # Check for mean structure
  if(!lavaan_model@Options$meanstructure){
    correct_call <- lavaan_model@call
    correct_call$meanstructure = TRUE
    stop("Your lavaan_model must have a meanstructure. Use: \n",
         paste0(deparse(correct_call), collapse = "\n"), "\n")
  }

  # check for multi-group
  if(lavaan_model@Model@ngroups != 1)
    stop("Currently multi-group models are not supported")

  # check for categorical variables
  if(lavaan_model@Model@categorical)
    stop("Currently models with categorical variables are not supported")

  if(!phantom_free %in% c("variance", "none"))
    stop("Currently only phantom_free = 'variance' or 'none' is supported.")
}
