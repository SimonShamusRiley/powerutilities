#' Performs F-test on coefficients of a model fit via glmmTMB using one of
#' several available methods for calculating denomenator degrees of freedom
#'
#' @param model a \code{\link{glmmTMB}} model object
#' @param method Currently supported options are: "containment" (the default),
#'   "nlme" and "within-between"
#' @param type Either 2 (default) or 3. Note that although the names are the
#'   same, the meaning of "type 2" and "type 3" are different from what they
#'   mean in SAS. See: `help(car::Anova)` for details.
#' @param test.statistic Either 'F' or 'Chisq'.
#' @param contr_sum Logical. Refit model using sum-to-zero contrasts? Default is
#'   `TRUE`.
#'
#' @return a \code{\link{data.frame}}
#' 
#' @exportS3Method
anova.glmmTMB <- function(model, method = "containment", type = 2, test.statistic = "F", contr_sum = TRUE){
  
  #if (class(model) != "glmmTMB") {
  #  stop ("Only glmmTMB models are supported")
  #}
  
  #TODO: check to see if this makes a difference
  if (contr_sum == TRUE){
    current_contrast_settings <- options("contrasts")
    options(contrasts = c("contr.sum", "contr.poly"))
    on.exit(options(current_contrast_settings))
  }
  
  if (type == "III" || type == 3) {
    type = 3
  } else if (type == "II" || type == 2) {
    type = 2
  } else {
    stop ("Specified type not supported at this time")
  }
  
  if(test.statistic == "F"){
    if(method == "nlme") {
      aov_out <- nlme_aov(model, type)
    } else if (method == "inner-outer") {
      aov_out <- inner_outer_aov(model, type)
    } else if (method == "containment") {
      aov_out <- containment_aov(model, type)
    } else {
      stop ("Only nlme, inner-outter, and containment methods are supported at this time")
    }
  } else if(test.statistic == "Chisq") {
    aov_out <- suppressForeignCheck(glmmTMB:::Anova.glmmTMB(model, type = type))
  } else {
    cat("Only F and Chisq test statistics are supported at this time")
  }
  
  return(aov_out)
}
