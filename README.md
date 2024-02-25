
<!-- README.md is generated from README.Rmd. Please edit that file -->

> **Warning**: bnSEM is very experimental. There are probably already
> better packages out there that do the same thing. The main objective
> of bnSEM is to explore similarities between SEM and Bayesian Networks.

# bnSEM

Structural Equation Models (SEM) and Bayesian Networks are closely
related. The objective of bnSEM (**B**ayesian **N**etwork **SEM**) is to
easily transition from SEM to Bayesian Networks. To this end, **bnSEM**
translates SEM fitted with `OpenMx` to Bayesian Networks in `bnlearn`.

## Installation

You can install the development version of bnSEM from
[GitHub](https://github.com/jhorzek/bnSEM/) with:

``` r
# install.packages("devtools")
devtools::install_github("jhorzek/bnSEM")
```

## Example

We first estimate the SEM with **OpenMx**.

``` r
library(mxsem)
model <- '
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'

mx_model <- mxsem(model,
                  data = OpenMx::Bollen) |>
  OpenMx::mxTryHard()
```

Next, we translate the model to a Bayesian Network:

``` r
library(bnSEM)
network <- bnSEM::bnSEM(mx_model = mx_model)
```

To get an impression of the network, you can create a plot:

``` r
library(qgraph)
qgraph::qgraph(network$dag)
```

![](man/figures/pol_dem_plot.png)

## Inference

With our Bayesian Network, we can now investigate the conditional
distribution of variables in our model (see
`vignette("Inference-Example", package = "bnSEM")` for a comparison to
rejection sampling)

``` r
# Check conditional distribution
# probability that dem65 in (1,2) given dem60 > 1:
bnlearn::cpquery(fitted = network$bayes_net,
                 event = (dem65 > 1 & dem65 < 2),
                 evidence = (dem60 > 1))
#> [1] 0.3442592

# Get distribution under this assumption:
dist <- bnlearn::cpdist(fitted = network$bayes_net,
                        node = "dem65",
                        evidence = (dem60 > 1))
hist(dist$dem65)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

``` r

# Now check dem60 < -1
dist <- bnlearn::cpdist(fitted = network$bayes_net,
                        node = "dem65",
                        evidence = (dem60 < -1))
hist(dist$dem65)
```

<img src="man/figures/README-unnamed-chunk-5-2.png" width="100%" />

To check our Bayesian Network, we can also simulate data from the
network and refit our SEM to check if the estimates align:

``` r
sim <- bnlearn::rbn(x = network$bayes_net, n = 100000)

fit_sim <- mxsem(model,
                 data = sim[,mx_model$manifestVars]) |>
  OpenMx::mxTryHard()
```

``` r
round(coef(fit_sim) - coef(mx_model), 3)
#>    ind60→x2    ind60→x3 ind60→dem60 ind60→dem65           a           b 
#>     -50.396     -42.899     -38.483     -14.113      -0.035      -0.026 
#>           c dem60→dem65       y1↔y1       y2↔y2       y3↔y3       y2↔y4 
#>      -0.037       0.025       0.021      -0.547       0.029      -0.317 
#>       y4↔y4       y2↔y6       y6↔y6       x1↔x1       x2↔x2       x3↔x3 
#>      -0.250      -0.016      -0.002       0.480       0.059      -0.065 
#>       y1↔y5       y5↔y5       y3↔y7       y7↔y7       y4↔y8       y6↔y8 
#>       0.000       0.003       0.040      -0.020       0.014      -0.017 
#>       y8↔y8 ind60↔ind60 dem60↔dem60 dem65↔dem65      one→y1      one→y2 
#>      -0.050      -0.448       0.006      -0.009      -0.025       0.005 
#>      one→y3      one→y4      one→y6      one→x1      one→x2      one→x3 
#>      -0.011      -0.003      -0.015      -0.006      -0.009      -0.010 
#>      one→y5      one→y7      one→y8 
#>      -0.001      -0.020      -0.009
```

## Central Challenge

When transitioning from SEM to Bayesian Network, a central challenge is
the use of covariances as parameters in SEM. To translate those
covariances to a model that can be easily used with Bayesian Networks,
bnSEM replaces covariances with direct effects of unobserved phantom
variables. The approach is explained in more detail, for instance, by
Merkle & Rosseel (2015; see p. 8).

Merkle, E. C., & Rosseel, Y. (2015). blavaan: Bayesian structural
equation models via parameter expansion. arXiv preprint
arXiv:1511.05604.
