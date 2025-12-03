#' @title Performs F-test on coefficients of a model fit via glmmTMB using one
#'   of several available methods for calculating denomenator degrees of freedom
#'
#' @param model a \code{\link{glmmTMB}} model object
#' @param method Currently supported options are: "containment" (the default),
#'   "pinheiro-bates", "asymptotic" and "satterthwaite".
#' @param type Either 2 (default) or 3 (note that either roman or arabic
#'   numerals are acceptable). Type II tests each term after all other
#'   non-higher-order terms, whereas type III tests each term after all other
#'   terms, including higher-order terms.
#' @param contr_sum Logical. Refit model using sum-to-zero contrasts? Default is
#'   `TRUE`.
#'
#' @return a \code{\link{data.frame}}
#'
#' @exportS3Method
anova.glmmTMB <- function(model, method = "containment", type = 2, contr_sum = TRUE){
  
  #TODO: check to see if this makes a difference
  if (contr_sum == TRUE){
    current_contrast_settings <- options("contrasts")
    options(contrasts = c("contr.sum", "contr.poly"))
    on.exit(options(current_contrast_settings))
  }
  
  if (length(ranef(model)$cond) == 1 & method != 'residual'){
    method = 'residual'
    warning('No random effects are present in the model: switching to residual degrees of freedom')
  }
  
  if (type == "III" || type == 3) {
    type = 3
  } else if (type == "II" || type == 2) {
    type = 2
  } else {
    stop ("type must be either 2 (testing each term after all other non-higher-order terms) or type 3 (testing each term after all other terms)")
  }
  
   if(method == "pinheiro-bates") {
      aov_out <- nlme_aov(model, type)
    } else if (method == "asymptotic") {
      aov_out <- asymptotic_aov(model, type)
    } else if (method == "containment") {
      aov_out <- containment_aov(model, type)
    } else if (method == 'satterthwaite'){
      aov_out <- satterthwaite_aov(model, type)
    } else if (method == 'residual') {
      aov_out <- residual_aov(model, type)
    } else {  
      stop ("Only containment, pinheiro-bates, satterthwaite, residual, and asymptotic methods are supported at this time")
    }
  
  return(aov_out)
}
