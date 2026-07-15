#' @title Identify how random effects terms will be ordered in the fitted glmmTMB model
#' @description
#' When performing a "plug-in" power analysis, the values for the random effects 
#' must be supplied in the order in which they are saved internally in the glmmTMB 
#' model. It is not always clear (especially where there are correlated random
#' slopes and intercepts) just what this ordering is. This helper function takes 
#' a model formula and a data set and outputs a table enumerating all the implied random
#' effects and their ordering.
#' 
#' @param formula A model formula (containing random effects).
#' @param data A data frame containing all the columns included in `formula`.
#' @param ... Other arguments passed to `glmmTMB`.
#' 
#' @examples
#' 
#' library(dplyr)
#' library(agridat)
#' library(powerutilities)
#' 
#' # Following the original analysis by Cornelius and Archbold, row x spacing 
#' # is treated as main-plots, and stock within main-plots as sub-plots.
#' apple_des = agridat::archbold.apple |> 
#'   mutate(across(rep:gen, factor), 
#'          mp = interaction(row, spacing, sep = '-'), 
#'          sp = interaction(row, spacing, stock, sep = '-'))
#' 
#' # If we opt to treat blocks as random, then here is the
#' # appropriate model for a split-split-plot RCBD:
#' ssp_mod = ~ spacing*stock*gen + (1|rep/mp/sp)
#' 
#' # The first theta value corresponds to to the standard deviation among 
#' #  sub-plots, the second theta value correspondss to the of the standard deviation
#' #  main-plots and the third theta value corresponds to the standard deviation
#' #  among blocks:
#' theta_finder(formula = ssp_mod, data = apple_des)
#' 
#' # Purely for the sake of illustration, lets add hypothetical random slopes
#' # for soil pH at subplot level. In the data, the values don't need to be 
#' # well simulated or even meaningful:
#' apple_des2 = apple_des |> 
#'   mutate(pH = 7)
#' 
#' # Note that mp and sp are already explicitely defined in terms of their
#' # nesting structure, so defining them here is terms of their interactions
#' # isn't really neccessary
#' ssp_cov_mod = ~ spacing*stock*gen + (1|rep) + (1|rep:mp) + (pH + 1|rep:mp:sp)
#' 
#' # Not one but two new theta values have been added to the model: we also
#' # have the correlation between the pH slope and split-plot intercepts.
#' theta_finder(formula = ssp_cov_mod, data = apple_des2)
#' 
#' # Finally, be aware that redundant random effects terms, which would simply 
#' # be consolidated when fitting a real model, are problematic here. In this example, 
#' # rep is included twice:
#' bad_mod = ~ spacing*stock*gen + (1|rep) + (1|rep/mp)
#' 
#' try(theta_finder(formula = bad_mod, data = apple_des))

#' @importFrom glmmTMB glmmTMB glmmTMBControl 
#' @importFrom nlme VarCorr
#' @importFrom reformulas formatVC findbars
#' @export
theta_finder = function(formula, data, ...){
  dots = list(...)
  if ('doFit' %in% names(dots)){
    dots = dots[!which(names(dots) == 'doFit')]
  }
  if (is.null(findbars(formula))){
    message(simpleMessage('The formula does include random effects'))
    return(invisible())
  }
  formula = update(formula, rep(1, nrow(data)) ~ .)
  
  args = c(list(formula = formula, data = data, 
                dispformula = ~ 0, doFit = F), 
           dots)
  def0 = do.call(what = glmmTMB, args = args)
  
  relist = struc_relist(def0)
  key = paste(names(relist), unlist(relist))
  
  if (any(duplicated(key))){
    dups = unique(key[which(duplicated(key))])
    dup_mssg = paste('The following random effects are duplicated in the formula: ', paste(dups, collapse = ", "))
    stop(simpleError(dup_mssg))
  }
  
  terms = lapply(relist, struc_reterms)
  counter = 1
  
  for (i in 1:length(terms)){
    terms[[i]]$Group[1, 1] = names(terms)[i]
    
    sd_start = counter
    sd_end = counter + (nrow( terms[[i]]$`Std.Dev`)-1)
    
    terms[[i]]$`Std.Dev.`[, 2] = seq(sd_start, sd_end)
    
    cor_start = sd_end+1
    cv = attr(terms[[i]], 'covstruct')
    if (cv == 'us'){
      ncov = sum(lower.tri(terms[[i]]$`Cor.`, diag = F))
      
      if (ncov > 0){
        cor_end = cor_start + (ncov-1)
        terms[[i]]$`Cor.`[lower.tri(terms[[i]]$`Cor.`, diag = F)] = seq(cor_start, cor_end)
      } else {
        cor_end = cor_start-1
      }
      
    } else {
      cor_end = cor_start
      terms[[i]]$`Cor.`[1, 1] = cor_start
    }
    
    counter = cor_end + 1
  }
  
  class(terms) = 'retermslist'
  return(terms)
}

#' @title Fit a glmmTMB model with fixed random effects terms
#' @description
#' This is a wrapper to help with fixing random effects terms in glmmTMB model
#' for use in subsequent power analysis. In particular, it assists and converting
#' random effect standard deviation and correlations into the log standard 
#' deviation and scaled log cholesky factors that glmmTMB expects. 
#' @param formula A two-sided model formula.
#' @param data A data.frame.
#' @param re_terms A vector of random effects standard deviations and, possibly, correlations. 
#' @param REML Logical. Whether to use REML estimation (default) or, alternatively, ML.
#' @param ... Other arguments passed to glmmTMB.
#' @importFrom glmmTMB glmmTMB glmmTMBControl
#' @export
set_glmm = function(formula, data, re_terms = NULL, disp = NULL,
                    REML = TRUE, ...){
  dots = list(...)
  
  if ('doFit' %in% names(dots)){
    dots = dots[!which(names(dots) == 'doFit')]
  }
  
  args0 = c(list(formula = formula, data = data, REML = REML, doFit = F))
  
  def0 = suppressWarnings(do.call(what = glmmTMB, args = c(args0, dots)))
  
  n_disp = length(def0$parameters$betadisp)
  n_disp_terms = length(disp)
  
  if (n_disp != n_disp_terms){
    stop(simpleError(paste0(n_disp, ' dispersion parameters are required but ', n_disp_terms, ' were supplied.')))
  }
  
  if (n_disp > 0){
    starts = list(betadisp = log(disp))
    maps = list(betadisp = factor(rep(NA, n_disp)))
  } else {
    starts = list()
    maps = list()
  }
  
  n_theta = length(def0$parameters$theta)
  n_re_terms = length(re_terms)
  if (n_theta != n_re_terms){
    stop(simpleError(paste0(n_theta, ' re_terms are required but ', n_re_terms, ' were supplied. Use `theta_finder()` for help with specifying random effects terms.')))
  }
  
  if (n_theta == 0){
    if (n_re_terms > 0){
      message(simpleMessage('re_terms is ignored: formula does not include any random effects'))
    }
      if (n_disp == 0){
        args = c(list(formula = formula, data = data), dots) 
      } else {
        args = c(list(formula = formula, data = data, 
                      start = starts, map = maps), 
                 dots) 
      }
  } else {
  
  vc = theta_finder(formula, data)
  
  trans_re_terms = c()
  
  for (i in 1:length(vc)){
    trans_re_terms = c(trans_re_terms, log(re_terms[as.numeric(vc[[i]]$Std.Dev.[,'Std.Dev.'])]))
    
    n = nrow(vc[[i]]$Std.Dev.)
    cvst = attr(vc[[i]], 'covstruct')
    if (cvst == 'diag' | (cvst == 'us' & n == 1)){next}
    
    cor_vals = re_terms[na.omit(as.numeric(vc[[i]]$Cor.))]
    cor_thetas = cor_convert_dispatch[[cvst]](x = cor_vals, n = n)
    trans_re_terms = c(trans_re_terms, cor_thetas)
  }
  
  starts = c(starts, list(theta = trans_re_terms))
  maps = c(maps, list(theta = factor(rep(NA, n_re_terms))))
  args = c(list(formula = formula, data = data, REML = REML,
                start = starts, map = maps), 
           dots)
  }
  do.call(glmmTMB, args)
}

#' @title Statistical Power of F-tests Performed on Models Fit with glmmTMB
#'
#' @description This function calculates the power of F-tests for each fixed
#'   effect term in the model using one of several methods for determining the
#'   denominator degrees of freedom.
#'
#' @param mod A \code{\link{glmmTMB}} model
#' @param ddf Either a method for calculating denominator degrees of freedom (
#'   currently supported options are "df.residual", "asymptotic" and
#'   "kenward-roger"), a numeric vector, or NULL (the default), in which case a
#'   method is selected based on the model type.
#' @param alpha The nominal type I error rate. Defaults to 0.05.
#' @param ... Other values passed to emmeans.
#'                
#' @importFrom glmmTMB glmmTMB glmmTMBControl
#' @importFrom emmeans emmeans joint_tests
#' @importFrom reformulas findbars RHSForm                
#' @export
power_ftest = function(mod, ddf = NULL, alpha = 0.05, ...){
  check_ddf(ddf)
  
  df_final = resolve_ddf(mod, ddf)
  
  numddf = all(inherits(ddf, 'numeric'))
  fe_form = nobars(RHSForm(formula(mod), as.form = T))
  
  dots = list(...)
  args = c(list(mod, fe_form), dots)
  emm = do.call(emmeans, args)
  
  emm@dffun = df_final$dffun
  emm@dfargs = df_final$dfargs
  
  jt = joint_tests(emm) |> 
    as.data.frame()
  
  if (numddf){
    if (!length(ddf) %in% c(1, nrow(jt))) {
      stop(simpleError(sprintf('%s ddf supplied for %s tests', length(ddf), nrow(jt))))
    }
    
    jt = jt |> 
      dplyr::mutate(df2= ddf, p.value = 1-pf(F.ratio, df1, df2))
  }
  
  pow = jt |>
    dplyr::rename(Term = `model term`, NumDF = df1, DenDF = df2,
                  Fval = F.ratio, Pval = p.value) |> 
    dplyr::mutate(NC_param = Fval*NumDF,
                  Fcrit = qf(1-alpha, NumDF, DenDF, 0),
                  Power = 1-pf(Fcrit, NumDF, DenDF, ncp = NC_param)) |>
    dplyr::select(Term, NumDF, DenDF,
                  Fval, Fcrit, Pval, Power)
  attr(pow, 'alpha') = alpha
  attr(pow, 'ddf') = df_final$ddf
  
  class(pow) = c('powertable', 'data.frame')
  return(pow)
}
#' @title Power of Contrasts Performed on Models Fit with glmmTMB
#'
#' @description This function calculates the power of contrasts.
#'
#' @param emm An emmGrid object associated with a `glmmTMB` model.
#' @param contr_list A (named) list of contrast specifications.
#' @param ddf Either a method for calculating denominator degrees of freedom (
#'   currently supported options are "df.residual", "asymptotic" and
#'   "kenward-roger"), a numeric value, or NULL (the default), in which case a
#'   method is selected based on the model type.
#' @param alpha  Numeric. The nominal type I error rate. Defaults to 0.05.
#' @param n_sims Numeric. The number of simulations to use for calculating
#'               type M error rate. If set to zero, the closed form, asymptotic
#'               calculations are used.
#' @param ... Other arguments passed to emmeans::contrast.         
#' 
#' @importFrom emmeans emmeans contrast    
#' @importFrom glmmTMB glmmTMB
#' @importFrom retrodesign retrodesign retro_design_closed_form               
#' @importFrom reformulas findbars
#' @export
power_contrast = function(emm, contr_list, ddf = NULL, 
                          alpha = 0.05, n_sims = 1e4, ...){
  if (!inherits(emm, 'emmGrid')){
    stop(simpleError('"emm" must be the result of a call to emmeans()'))
  } 
  if (!inherits(alpha, 'numeric') | alpha > 1 | alpha < 0){
    stop(simpleError('alpha must be a numeric value between 0 and 1'))
  }
  if (!inherits(n_sims, 'numeric') | n_sims < 0){
    stop(simpleError('n_sims must be a non-negative numeric value'))
  }
  
  check_ddf(ddf)
  
  dots = list(...)
  
  if ('ratios' %in% names(dots)){
    if (isFALSE(dots$ratios)){
      message(simpleMessage('Setting `ratios = FALSE`: power calculations must be performed on the link scale.'))
      dots = dots[!names(dots)=='ratios']
    }
  }
  
  if ('null' %in% names(dots)){
    if (null != 0){
      message(simpleMessage('Setting `null = 0`: powerutilities does not yet support non-zero null hypotheses.'))
      dots$null = 0 
    }
  }
  
  if ('predict.type' %in% names(emm@misc)){
    if (emm@misc$predict.type != 'emmeans'){
      message(simpleMessage('Setting `type = "emmeans"`: power calculations must be performed on the link scale.'))
      emm = update(emm, type = 'emmeans')
    }
  }
  
  fixed = is.null(findbars(formula(emm@model.info$call)))
  numddf = inherits(ddf, 'numeric')
  if (numddf & length(ddf) > 1){
    warning(simpleWarning('multiple ddf values supplied, only the first will be used'))
  }
  
  gen = FALSE
  
  call_list = as.list(emm@model.info$call)
  if ('family' %in% names(call_list)){
    fam = eval(call_list$family)
    if (fam$family != 'gaussian' | fam$link != 'identity'){
      gen = TRUE
    }
  }
  
  if (identical(ddf, 'df.residual') | (is.null(ddf) & (fixed & !gen))){
    model = eval(emm@model.info$call)
    df_final = resolve_ddf(emm, request = ddf)
  }
  
  df_final = resolve_ddf(emm, request = ddf)
  
  if (numddf){
    emm = update(emm, df = ddf)
  } else {
    emm@dffun = df_final$dffun
    emm@dfargs = df_final$dfargs
  }
  
  con = do.call(emmeans::contrast, c(list(emm, contr_list, ratios = FALSE), dots)) |> 
    as.data.frame() |> 
    dplyr::rename_with(
      .fn = ~ dplyr::case_when(
        . == "z.ratio" ~ "t.ratio",
        . %in% c("lower.CL", "asymp.LCL") ~ "LCL",
        . %in% c("upper.CL", "asymp.UCL") ~ "UCL",
        TRUE ~ .
      ),
      .cols = everything()
    ) |> 
    dplyr::mutate(NumDF = 1, 
                  DenDF = df, 
                  Fval = t.ratio^2, 
                  NC_param = Fval*NumDF, 
                  Fcrit = qf(1-alpha, NumDF, DenDF, 0),
                  Power = 1-pf(Fcrit, NumDF, DenDF, ncp = NC_param)) |> 
    dplyr::rename_with(capwords, .cols = 1:SE) |> 
    dplyr::rename(Pval = p.value) |> 
    dplyr::select(Contrast:SE, 
                  any_of(c('LCL', 'UCL')), NumDF, DenDF,  
                  Fval, Fcrit, Pval, Power) 
  
  retro_fun = ifelse(n_sims == 0,
                     retro_design_closed_form, 
                     retrodesign)
  if (n_sims == 0){
    more_errs = with(con, mapply(retro_design_closed_form,
                                 A = Estimate, s = SE, 
                                 MoreArgs = list(alpha = alpha))) 
  } else {
    more_errs = with(con, mapply(retrodesign,
                                 A = Estimate, s = SE, df = DenDF,
                                 MoreArgs = list(alpha = alpha, 
                                                 n.sims = n_sims))) 
  }
  
  con$TypeS  = unlist(more_errs['type_s', ])
  con$TypeM = unlist(more_errs['type_m', ])
  
  out = con |> 
    dplyr::mutate(TypeM = ifelse(abs(Estimate) < 1.5e-8, Inf, TypeM))
  
  attr(out, 'alpha') = alpha
  attr(out, 'ddf') = df_final$ddf
  
  class(out) = c('powertable', 'data.frame')
  return(out)
}



