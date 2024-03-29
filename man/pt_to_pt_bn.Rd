% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pt_to_pt_bn.R
\name{pt_to_pt_bn}
\alias{pt_to_pt_bn}
\title{pt_to_pt_bn}
\usage{
pt_to_pt_bn(parameter_table, mx_model_int)
}
\arguments{
\item{parameter_table}{parameter table}

\item{mx_model_int}{mxModel with replaced covariances}
}
\value{
parameter table for bnlearn
}
\description{
Extracts the regressions underlying the SEM and puts them in a parameter
table that bnlearn understands. See https://www.bnlearn.com/examples/custom-fitted/
(section "Continuous networks").
}
\keyword{internal}
