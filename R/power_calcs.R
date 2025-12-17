#' @title Statistical Power of F-tests Performed on Models Fit with glmmTMB
#'
#' @description This function calculates the power of F-tests for each fixed
#'   effect term in the model using one of several methods for determining the
#'   denominator degrees of freedom |>
#'
#' @param mod A \code{\link{glmmTMB}} model
#' @param ddf Method for calculating denominator degrees of freedom: containment
#'   (default), pinheiro-bates, satterthwaite, or asymptotic.
#' @param type Either 2 (default) or 3 (note that roman or arabic numerals are
#'   accepted). Type II tests each term after testing all other non-higher-order
#'   terms, whereas type III tests each term after all others, including
#'   higher-order terms.
#' @param alpha The nominal type I error rate. Defaults to 0.05.
#' @export
power_ftest = function(mod, ddf = 'containment', type = 2, alpha = 0.05) {
  
  anova(mod, method = ddf, type = type) |> 
    tibble::rownames_to_column(var = 'factor') |> 
    dplyr::mutate(alpha = alpha, 
                  nc_param = Fvalue*numDF, 
                  F_crit = qf(1-alpha, numDF, denDF, 0),
                  power = 1-pf(F_crit, numDF, denDF, ncp = nc_param)) |> 
    dplyr::select(factor, numDF, denDF, alpha, Fvalue, pvalue, 
                  nc_param, F_crit, power)
  
}  

#' @title Power of Contrasts Performed on Models Fit with glmmTMB
#'
#' @description This function calculates the power of F-tests for each fixed
#'   effect term in the model using one of several methods for determining the
#'   denominator degrees of freedom |>
#'
#' @param emm An emmGrid object associated with a \code{\link{glmmTMB}} model.
#' @param contr_list A named list of contrasts specifications.
#' @param ddf Numeric. The value to use for the denominator degrees of freedom.
#'   Should probably correspond to the value returned from a call to
#'   power_ftest.
#' @param alpha The nominal type I error rate. Defaults to 0.05.
#' @export
power_contrast = function(emm, contr_list, ddf, alpha = 0.05){
  emmeans::contrast(emm, contr_list, df = ddf) |> 
    as.data.frame() |> 
    rename_with(~ if_else(. == "z.ratio", "t.ratio", .), .cols = everything()) |> 
    dplyr::mutate(numDF = 1, 
                  denDF = df, 
                  Fvalue = t.ratio^2, 
                  nc_param = Fvalue*numDF, 
                  alpha = alpha, 
                  F_crit = qf(1-alpha, numDF, denDF, 0),
                  power = 1-pf(F_crit, numDF, denDF, ncp = nc_param)) |> 
    dplyr::rename(pvalue = p.value) |> 
    dplyr::select(contrast:SE, numDF, denDF, alpha, Fvalue, pvalue, 
                  nc_param, F_crit, power)
}
