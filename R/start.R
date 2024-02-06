library(lavaan)
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
'

lavaan_model <- sem(model,
                    data = PoliticalDemocracy,
                    meanstructure = TRUE)

# Extract Parameter Table
if(!lavaan_model@Options$meanstructure){
  correct_call <- lavaan_model@call
  correct_call$meanstructure = TRUE
  stop("Your lavaan_model must have a meanstructure. Use: \n", paste0(deparse(correct_call), collapse = "\n"), "\n")
}

parameter_table <- lavaan_model@ParTable |>
  as.data.frame()

# set labels for all parameters
for(p in 1:nrow(parameter_table)){
  if(parameter_table$label[p] != "")
    next

  parameter_table$label[p] <- paste0(parameter_table$lhs[p],
                                     parameter_table$op[p],
                                     parameter_table$rhs[p])
}

# remove equality constraints
parameter_table <- parameter_table[!parameter_table$op %in% c("<", ">", "=="),]

parameter_table <- parameter_table[,c("lhs", "op", "rhs", "label", "est")]

# we need to replace all effects with directed effects
variables <- unique(unlist(parameter_table[,c("lhs", "rhs")]))
variables <- variables[variables != ""]

parameter_table_bn <- c()

regressions <- vector("list", length(variables))
names(regressions) <- variables

for(v in variables){
  # we want to create a linear model for each of the variables
  pars <- list()
  # extract intercepts
  if(any((parameter_table$lhs == v) & (parameter_table$op == "~1"))){
    pars[["(Intercept)"]] <- parameter_table$est[parameter_table$lhs == v & parameter_table$op == "~1"]
  }else{
    pars[["(Intercept)"]] <- 0.0
  }

  # get all loadings and regressions for this variable
  for(reg in which((parameter_table$lhs == v) & (parameter_table$op == "~"))){
    pars[[parameter_table$rhs[reg]]] <- parameter_table$est[reg]
  }
  for(lam in which((parameter_table$rhs == v) & (parameter_table$op == "=~"))){
    pars[[parameter_table$lhs[lam]]] <- parameter_table$est[lam]
  }
  pars <- list(coef = unlist(pars))

  # extract residual variance
  if(any((parameter_table$lhs == v) & (parameter_table$op == "~~") & (parameter_table$rhs == v))){

    # check if the item also has a covariance
    if(any((parameter_table$lhs == v) & (parameter_table$op == "~~") & (parameter_table$rhs != v)) |
       any((parameter_table$lhs != v) & (parameter_table$op == "~~") & (parameter_table$rhs == v))){
      pars[["sd"]] <- 0
      # add additional latent variable instead
      regressions[[paste0("phantom_lv_", v)]] <- list(
        coef = c("(Intercept)" = 0.0),
        sd = sqrt(parameter_table$est[(parameter_table$lhs == v) &
                                        (parameter_table$op == "~~") &
                                        (parameter_table$rhs == v)])
      )
    }else{
      pars[["sd"]] <- sqrt(parameter_table$est[(parameter_table$lhs == v) &
                                                 (parameter_table$op == "~~") &
                                                 (parameter_table$rhs == v)])
    }
  }else{
    pars[["sd"]] <- 0.0
  }
  regressions[[v]] <- pars
}

# now, we have to implement additional latent variables for all covariances
for(i in which((parameter_table$op == "~~") & (parameter_table$lhs != parameter_table$rhs))){
  new_latent <- parameter_table$label[i]
  regressions[[new_latent]][["coef"]] <- c("(Intercept)" = 0.0)
  regressions[[new_latent]][["sd"]]   <- sqrt(parameter_table$est[i])

  add_coef <- c(1)
  names(add_coef) <- parameter_table$label[i]
  regressions[[paste0("phantom_lv_", parameter_table$lhs[i])]][["coef"]] <- c(regressions[[paste0("phantom_lv_", parameter_table$lhs[i])]][["coef"]], add_coef)
  regressions[[paste0("phantom_lv_", parameter_table$rhs[i])]][["coef"]] <- c(regressions[[paste0("phantom_lv_", parameter_table$rhs[i])]][["coef"]], add_coef)
}

# create network
bn_model <- c()
for(i in 1:length(regressions)){

  if(length(regressions[[i]][["coef"]]) == 1){
    # it's an exogenous variable!
    bn_model <- c(bn_model,
                  paste0("[", names(regressions[i]), "]"))
  }else{
    preds <- names(regressions[[i]][["coef"]])
    preds <- preds[preds != "(Intercept)"]
    bn_model <- c(bn_model,
                  paste0("[", names(regressions[i]), "|", paste0(preds, collapse = ":"), "]"))
  }
}
bn_model <- paste0(bn_model, collapse = "")

# create dag
dag <- bnlearn::model2network(string = bn_model)
plot(dag)

# set up parameters
cfit = bnlearn::custom.fit(dag,
                           dist = regressions)

sim <- bnlearn::rbn(x = cfit, n = 10000)

fit_sim <- sem(model,
               data = sim[,lavaan_model@Data@ov.names[[1]]],
               meanstructure = TRUE)
round(coef(fit_sim) -
        coef(lavaan_model), 3)

bnlearn::cpdist(fitted = cfit,
                nodes = "dem60",
                evidence = list("ind60" = 3),
                method = "lw")

warning("Covariances are implemented incorrectly. They should point to the variances of the residuals!")
