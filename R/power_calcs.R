#' @title Statistical Power of F-tests Performed on Models Fit with glmmTMB
#'
#' @description This function calculates the power of F-tests for each fixed
#'   effect term in the model using one of several methods for determining the
#'   denominator degrees of freedom.
#'
#' @param mod A \code{\link{glmmTMB}} model
#' @param ddf Method for calculating denominator degrees of freedom 
#'  ('kenward-roger' (default), 'df.residual', or 'asymptotic') or else
#'  a numeric vector of DF values. 
#' @param alpha The nominal type I error rate. Defaults to 0.05.
#' 
#' @examples#' 
#' # Power analysis in a split-plot RCBD experiment on oat yield
#' oats = agridat::yates.oats |> 
#'        dplyr::mutate(dplyr::across(nitro, factor))
#' 
#' # All random effects are arbitrarily given standard deviations of 10
#' oat_mod = glmmTMB(yield ~ gen*nitro + (1|block/gen), 
#'                   data = oats, 
#'                   REML = T, 
#'                   dispformula = ~0, 
#'                   start = list(theta = c(log(c(10, 10)))),
#'                   map = list(theta = factor(c(NA, NA))), 
#'                   control = glmmTMBControl(zerodisp_val = log(10))) 
#'                   
#' # In a perfectly balanced case, Kenward-Roger produces the same denominator
#' # degrees of freedom as the classical containment method
#' power_ftest(oat_mod, ddf = 'kr') 
#' 
#' # Residual degrees of freedom, appropriate for non-nested designs, is 
#' # anti-conservative
#' power_ftest(oat_mod, ddf = 'df.residual')  
#' 
#' # Wald-type asympotic tests are the most anti-conservative
#' power_ftest(oat_mod, ddf = 'asymptotic')                     
#'                   
#' @export
power_ftest = function(mod, ddf = 'kenward-roger', alpha = 0.05){
  require(glmmTMB)
  require(emmeans)
  
  if (!inherits(ddf, c('character', 'numeric'))){
    stop(simpleError('ddf must be either a numeric vector or else one of: "kenward-roger", "asymptotic", "df.residual"'))
  }
  
  if (inherits(ddf, 'character')) {
    ddf = match.arg(ddf, choices = c('kenward-roger', 'kr', 
                                     'asymptotic', 'df.residual')) |> 
      switch('kenward-roger' = 'kenward-roger', 
             'kr' = 'kenward-roger', 
             'asymptotic' = 'asymptotic', 
             'df.residual' = 'df.residual')
    
    if (ddf %in% c('kenward-roger') & is.null(reformulas::findbars(formula(mod)))){
      message('No random effects are present in the model: switching to "df.residual"')
      ddf = 'df.residual'
    }
    
    jt = emmeans::joint_tests(mod, ddf = ddf) 
    
  } else {
    jt = emmeans::joint_tests(mod, df = 0) |> 
      dplyr::mutate(df2= ddf, 
                    p.value = 1-pf(F.ratio, df1, df2)) 
    
  }
  
  pow = jt |>
    as.data.frame() |>
    dplyr::rename(Term = `model term`, NumDF = df1, DenDF = df2,
                  Fval = F.ratio, Pval = p.value) |> 
    dplyr::mutate(NC_param = Fval*NumDF,
                  Fcrit = qf(1-alpha, NumDF, DenDF, 0),
                  Power = 1-pf(Fcrit, NumDF, DenDF, ncp = NC_param)) |>
    dplyr::select(Term, NumDF, DenDF,
                  Fval, Fcrit, Pval, Power)
  attr(pow, 'alpha') = alpha
  attr(pow, 'ddf') = ifelse(inherits(ddf, 'character'), ddf,
                            'user-specified')
  
  class(pow) = c('powertable', 'data.frame')
  return(pow)
}

#' @title Power of Contrasts Performed on Models Fit with glmmTMB
#'
#' @description This function calculates the power of contrasts.
#'
#' @param emm An emmGrid object associated with a \code{\link{glmmTMB}} model.
#' @param contr_list A (named) list of contrast specifications.
#' @param alpha  Numeric. The nominal type I error rate. Defaults to 0.05.
#' @param n_sims Numeric. The number of simulations to use for calculating
#'               type M error rate. If set to zero, the closed form, asymptotic
#'               calculations are used.
#' @param ... Other arguments passed to emmeans::contrast.         
#' 
#' @examples
#' # Power analysis in a split-plot RCBD experiment on oat yield
#' oats = agridat::yates.oats |> 
#'        dplyr::mutate(dplyr::across(nitro, factor))
#' 
#' # All random effects are arbitrarily given standard deviations of 10
#' oat_mod = glmmTMB(yield ~ gen*nitro + (1|block/gen), 
#'                   data = oats, 
#'                   REML = T, 
#'                   dispformula = ~0, 
#'                   start = list(theta = c(log(c(10, 10)))),
#'                   map = list(theta = factor(c(NA, NA))), 
#'                   control = glmmTMBControl(zerodisp_val = log(10))) 
#' 
#' # Here, the denominator DF are determined when EMMs are calculated:
#' (oat_emm = emmeans(oat_mod, ~ gen, ddf = 'kenward-roger'))
#' 
#' contr = list('GoldenRain vs Others' = c(1, -0.5, -0.5))
#' 
#' power_contrast(oat_emm, contr)             
#' 
#' # With manually specified denominator DFs (note that the argument
#' # is not `ddf =`, but `df =`)
#' (oat_emm2 = emmeans(oat_mod, ~ gen, df = 12))
#' 
#' power_contrast(oat_emm2, contr)  
#'                   
#' @export
power_contrast = function(emm, contr_list, alpha = 0.05, n_sims = 1e4, ...){
  require(glmmTMB)
  require(emmeans)
  require(retrodesign)
  
  if (!inherits(emm, 'emmGrid')){
    stop(simpleError('"emm" must be the result of a call to emmeans()'))
  } 
  if (!inherits(alpha, 'numeric') | alpha > 1 | alpha < 0){
    stop(simpleError('alpha must be a numeric value between 0 and 1'))
  }
  if (!inherits(n_sims, 'numeric') | n_sims < 0){
    stop(simpleError('n_sims must be a non-negative numeric value'))
  }
  
  dots = list(...)
  
  if ('ratios' %in% names(dots)){
    message(simpleMessage('Setting `ratios = FALSE`: power calculations must be performed on the link scale.'))
    dots = dots[which(names(dots) != 'ratios')]
  }
  
  if ('null' %in% names(dots)){
    message(simpleMessage('Setting `null = 0`: powerutilities does not yet support non-zero null hypotheses.'))
    dots$null = 0
  }
  
  if ('predict.type' %in% names(emm@misc)){
    if (emm@misc$predict.type != 'emmeans'){
    message(simpleMessage('Setting `type = "emmeans"`: power calculations must be performed on the link scale.'))
    emm = update(emm, type = 'emmeans')
    }
  }
  
  capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s, 1, 1)),
                             {s <- substring(s, 2); if(strict) tolower(s) else s},
                             sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
  }

  con = do.call(emmeans::contrast, c(list(emm, contr_list, ratios = FALSE), dots)) |> 
    as.data.frame() |> 
    dplyr::rename_with(~ dplyr::if_else(. == "z.ratio", "t.ratio", .),
                       .cols = everything()) |> 
    dplyr::mutate(NumDF = 1, 
                  DenDF = df, 
                  Fval = t.ratio^2, 
                  NC_param = Fval*NumDF, 
                  Fcrit = qf(1-alpha, NumDF, DenDF, 0),
                  Power = 1-pf(Fcrit, NumDF, DenDF, ncp = NC_param)) |> 
    dplyr::rename_with(capwords, .cols = 1:SE) |> 
    dplyr::rename(Pval = p.value) |> 
    dplyr::select(Contrast:SE, NumDF, DenDF,  
                  Fval, Fcrit, Pval, Power) 
  
  if (n_sims == 0){
    more_errs = with(con, mapply(retrodesign::retro_design_closed_form,
                                 A = Estimate, s = SE,
                                 MoreArgs = list(alpha = alpha))) 
  } else {
    more_errs = with(con, mapply(retrodesign::retrodesign, 
                                 A = Estimate, s = SE,
                                 df = DenDF, n.sims = n_sims,
                                 MoreArgs = list(alpha = alpha))) 
  }
   
  con$TypeS  = unlist(more_errs['type_s', ])
  con$TypeM = unlist(more_errs['type_m', ])
  
  out = con |> 
    mutate(TypeM = ifelse(abs(Estimate) < 1.5e-8, Inf, TypeM))
  
  attr(out, 'alpha') = alpha
  ddf = switch(deparse(emm@dffun)[2], 
               "pbkrtest::Lb_ddf(k, dfargs$unadjV, dfargs$adjV)" = 'kenward-roger', 
               'Inf' = 'asymptotic', 
               "stats::df.residual(dfargs$object)" = 'df.residual')
  
  if (ddf == 'asymptotic' & !is.null(emm@misc$df)){
    ddf = 'user-specified'
  }
  
  attr(out, 'ddf') = ddf
  
  class(out) = c('powertable', 'data.frame')
  return(out)
}

#' @title Print the results of a power analysis
#'
#' @description This function controls the formatting when printing `powertable`
#' objects.
#'
#' @param x A `powertable` object
#' @param digits Integer. The number of digits to print for all columns other 
#' than the p-value and power. 
#' @param pdigits Integer. The number of digits to print for the p-value and 
#' power columns
#' @export
print.powertable = function(x, digits = 1, pdigits = getOption('pdigits', default = 4), ...){
  
  fmt_like_pval = function(p) {
    ifelse(p < 10^(-pdigits),
           paste0('<.', paste(rep('0', pdigits - 1), collapse = ''), '1'),
           sprintf(paste0('%.', pdigits, 'f'), p))
  }
  
  out = x |>
    dplyr::mutate(dplyr::across(dplyr::any_of(c('Pval', 'Power', 'TypeS')), fmt_like_pval),
                  across(where(is.numeric), ~sprintf(paste0('%.', digits, 'f'), .)))
  
  print.data.frame(out)
  
  cat(paste0('\nDegrees-of-freedom method: ', attr(x, 'ddf'), '\n\u03B1 = ', attr(x, 'alpha')))
}

#' @exportS3Method knitr::knit_print
knit_print.powertable = function(x, ...) {
  out <- paste(capture.output(print.powertable(x, ...)), collapse = "\n")
  if (knitr::is_html_output()) {
    knitr::asis_output(paste0("<pre>", out, "</pre>"))
  } else {
    knitr::asis_output(paste0("\n```\n", out, "\n```\n"))
  }
}
