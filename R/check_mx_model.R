#' check_mx_model
#'
#' Check if the OpenMx model can be translated to bnlearn
#' @param mx_model fitted OpenMx model
#' @param optimize should the substitute model in case of covariances be optimized?
#' @importFrom methods is
#' @importFrom stats logLik
#' @keywords internal
check_mx_model <- function(mx_model,
                           optimize){

  if(!is(mx_model, "MxRAMModel"))
    stop("mx_model must be of type MxRAMModel.")

  if(length(mx_model$algebras) != 0)
    stop("bnSEM does not support mxModels with algebras.")

  if(is.na(logLik(mx_model)) & optimize)
    stop("The mx_model has a logLik of NA. Please optimize your model parameters ",
         "(e.g., using mxFit) before passing the model to bnSEM or use optimize = FALSE.")

  if(!all(c("A", "S", "F", "M") %in% names(mx_model$matrices)))
    stop("Expected OpenMx matrices to be called A, S, M, and F")
}
