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
#' @return A 'retermslist' object numbering each of the random effect variances
#' and correlations implied by the model formula.
#' 
#' @examplesIf requireNamespace("agridat", quietly = TRUE)
#' 
#' library(dplyr)
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
#' #  sub-plots, the second theta value corresponds to the of the standard deviation
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
#' # Note that mp and sp are already explicitly defined in terms of their
#' # nesting structure, so defining them here is terms of their interactions
#' # isn't really necessary
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
#' 
#' @importFrom stats update
#' @importFrom glmmTMB glmmTMB glmmTMBControl 
#' @importFrom reformulas findbars
#' @export
theta_finder = function(formula, data, ...){
  dots = list(...)
  if ('doFit' %in% names(dots)){
    dots = dots[names(dots) != 'doFit']
  }
  if (is.null(findbars(formula))){
    message(simpleMessage('The formula does not include random effects'))
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
#' 
#' @param formula A two-sided model formula.
#' @param data A data.frame.
#' @param re_terms A vector of random effects standard deviations and, possibly, correlations.
#' @param disp The value at which to fix the residual dispersion. 
#' @param REML Logical. Whether to use REML estimation (default) or, alternatively, ML.
#' @param ... Other arguments passed to glmmTMB.
#' 
#' @return A `glmmTMB` object.
#' 
#' @examples
#' 
#' # Create synthetic data set:
#' nrep = 8 # Number of replicates
#' nfac = 2 # Number of treatments
#' dat = expand.grid(Rep = factor(1:nrep),
#'                   Trt = factor(LETTERS[1:nfac])) 
#' dat$Y = ifelse(dat$Trt == 'A', 0, 3)
#' 
#' # "Fit" the model while setting residual standard deviation to 2:
#' mod = set_glmm(Y ~ Trt, data = dat, disp = 2)
#' 
#' # Confirm it worked:
#' sigma(mod)
#' 
#' @importFrom glmmTMB glmmTMB glmmTMBControl
#' @importFrom stats na.omit
#' @export
set_glmm = function(formula, data, re_terms = NULL, disp = NULL,
                    REML = TRUE, ...){
  dots = list(...)
  
  if ('doFit' %in% names(dots)){
    dots = dots[names(dots) != 'doFit']
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
  out = do.call(glmmTMB, args)
  return(out)
}

#' @title Statistical Power of F-tests Performed on Models Fit with glmmTMB
#'
#' @description This function calculates the power of F-tests for each fixed
#'   effect term in the model using one of several methods for determining the
#'   denominator degrees of freedom.
#'
#' @param mod A `glmmTMB` model
#' @param ddf Either a method for calculating denominator degrees of freedom (
#'   currently supported options are "df.residual", "asymptotic" and
#'   "kenward-roger"), a numeric vector, or NULL (the default), in which case a
#'   method is selected based on the model type.
#' @param alpha The nominal type I error rate. Defaults to 0.05.
#' @param ... Other values passed to emmeans.
#'
#' @return A `powertable` object.
#' 
#' @examples
#' 
#' # Create synthetic data set:
#' nrep = 8 # Number of replicates
#' nfac = 2 # Number of treatments
#' dat = expand.grid(Rep = factor(1:nrep),
#'                   Trt = factor(LETTERS[1:nfac])) 
#' dat$Y = ifelse(dat$Trt == 'A', 0, 3)
#' 
#' # "Fit" the model while setting residual standard deviation to 2:
#' mod = set_glmm(Y ~ Trt, data = dat, disp = 2)
#' 
#' # Calculate power
#' power_ftest(mod)
#' 
#' @importFrom glmmTMB glmmTMB glmmTMBControl
#' @importFrom emmeans emmeans joint_tests
#' @importFrom reformulas findbars RHSForm  
#' @importFrom stats qf pf    
#' @importFrom dplyr mutate rename select
#' @importFrom rlang .data
#'       
#' @export
power_ftest = function(mod, ddf = NULL, alpha = 0.05, ...){
  check_ddf(ddf)
  
  df_final = resolve_ddf(mod, ddf)
  
  numddf = all(inherits(ddf, 'numeric'))
  fe_form = nobars(RHSForm(formula(mod), as.form = T))
  
  dots = list(...)
  args = c(list(object = mod, specs = fe_form), dots)
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
      mutate(df2 = ddf,  p.value = 1-pf(.data$F.ratio, .data$df1, .data$df2))
  }
  
  pow = jt |>
    rename(Term = "model term", NumDF = "df1", DenDF = "df2",
                  Fval = "F.ratio", Pval = "p.value") |> 
    dplyr::mutate(NC_param = .data$Fval*.data$NumDF,
                  Fcrit = qf(1-alpha, .data$NumDF, .data$DenDF, 0),
                  Power = 1-pf(.data$Fcrit, .data$NumDF, .data$DenDF, ncp = .data$NC_param)) |>
    dplyr::select('Term', 'NumDF', 'DenDF',
                  'Fval', 'Fcrit', 'Pval', 'Power')
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
#' @return A `powertable` object.       
#' 
#' @examples
#' 
#' library(emmeans)
#' 
#' # Create synthetic data set:
#' nrep = 8 # Number of replicates
#' nfac = 2 # Number of treatments
#' dat = expand.grid(Rep = factor(1:nrep),
#'                   Trt = factor(LETTERS[1:nfac])) 
#' dat$Y = ifelse(dat$Trt == 'A', 0, 3)
#' 
#' # "Fit" the model while setting residual standard deviation to 2:
#' mod = set_glmm(Y ~ Trt, data = dat, disp = 2)
#' 
#' # Create emmeans object:
#' emm = emmeans(mod, ~ Trt)
#' 
#' # Define contrasts (here there is only one):
#' contr = list('A-B' = c(1, -1))
#' 
#' # Calculate power of the contrast
#' power_contrast(emm, contr)
#' 
#' @importFrom emmeans emmeans contrast    
#' @importFrom glmmTMB glmmTMB
#' @importFrom retrodesign retrodesign retro_design_closed_form               
#' @importFrom reformulas findbars
#' @importFrom stats qf pf update formula
#' @importFrom dplyr mutate rename rename_with select case_when any_of everything 
#' @importFrom rlang .data
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
    if (dots$null != 0){
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
  
  ## I think this is to remove
  # if (identical(ddf, 'df.residual') | (is.null(ddf) & (fixed & !gen))){
  #   model = eval(emm@model.info$call)
  #   df_final = resolve_ddf(emm, request = ddf)
  # }
  
  df_final = resolve_ddf(emm, request = ddf)
  
  if (numddf){
    emm = update(emm, df = ddf)
  } else {
    emm@dffun = df_final$dffun
    emm@dfargs = df_final$dfargs
  }
  
  con = do.call(contrast, c(list(emm, contr_list, ratios = FALSE), dots)) |> 
    as.data.frame() |> 
    rename_with(
      .fn = ~ case_when(
        . == "z.ratio" ~ "t.ratio",
        . %in% c("lower.CL", "asymp.LCL") ~ "LCL",
        . %in% c("upper.CL", "asymp.UCL") ~ "UCL",
        TRUE ~ .
      ),
      .cols = everything()
    ) |> 
    mutate(NumDF = 1, 
           DenDF = .data$df, 
           Fval = .data$t.ratio^2, 
           NC_param = .data$Fval*.data$NumDF, 
           Fcrit = qf(1-alpha, .data$NumDF, .data$DenDF, 0),
           Power = 1-pf(.data$Fcrit, .data$NumDF, .data$DenDF, ncp = .data$NC_param)) |> 
    rename_with(capwords, .cols = 1:SE) |> 
    rename(Pval = "p.value") |> 
    select(Contrast:SE, 
           any_of(c('LCL', 'UCL')), 'NumDF', 'DenDF',  
           'Fval', 'Fcrit', 'Pval', 'Power') 
  
  if (n_sims == 0){
    more_errs = mapply(retro_design_closed_form,
                       A = con$Estimate, s = con$SE, 
                       MoreArgs = list(alpha = alpha))
  } else {
    more_errs = mapply(retrodesign,
                       A = con$Estimate, s = con$SE, df = con$DenDF,
                       MoreArgs = list(alpha = alpha, n.sims = n_sims))
  }
  
  con$TypeS  = unlist(more_errs['type_s', ])
  con$TypeM = unlist(more_errs['type_m', ])
  
  out = con |> 
    mutate(TypeM = ifelse(abs(.data$Estimate) < 1.5e-8, Inf, .data$TypeM))
  
  attr(out, 'alpha') = alpha
  attr(out, 'ddf') = df_final$ddf
  
  class(out) = c('powertable', 'data.frame')
  return(out)
}



