#' add_parameter
#'
#' Returns data.frame with added parameter
#' @keywords internal
add_parameter <- function(id,
                          lhs,
                          op,
                          rhs,
                          user = 1,
                          block = 1,
                          group = 1,
                          free,
                          ustart = NA,
                          exo= 0,
                          label,
                          plabel,
                          start,
                          est,
                          se = NA){
  return(data.frame(id,
                    lhs,
                    op,
                    rhs,
                    user,
                    block,
                    group,
                    free,
                    ustart,
                    exo,
                    label,
                    plabel,
                    start,
                    est,
                    se))
}
