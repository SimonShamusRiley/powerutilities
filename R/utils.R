#' @noRd
struc_relist = function(x){
  reterms = x$condList$reTrms$cnms
  types = sapply(x$condReStruc, \(re){names(re$blockCode)})
  groups = names(reterms)
  out = mapply(`attr<-`,  x = reterms, value = types, 
               MoreArgs = list(which = 'covstruct'), SIMPLIFY = F)
  out = mapply(`attr<-`,  x = out, value = groups, 
               MoreArgs = list(which = 'group'), SIMPLIFY = F)
  return(out)
}

#' @noRd
covstruct_dispatch = list(
  us   = function(x, ...) {
    out = list(Group = array(NA, dim = c(length(x), 1), 
                       dimnames = list(NULL, 'Group')),
         `Std.Dev.` = array(rep(x, 2), dim = c(length(x), 2), 
                            dimnames = list(NULL, c('Name', 'Std.Dev.'))),
         `Cor.` = array(as.character(NA), dim = c(length(x), length(x)), 
                        dimnames = list(NULL, c('Cor.', rep('', length(x)-1)))))
    attr(out, 'covstruct') = 'us'
    return(out)
  },
  cs   = function(x, ...) {
    out = list(Group = array(NA, dim = c(length(x), 1), 
                       dimnames = list(NULL, 'Group')),
         `Std.Dev.` = array(rep(x, 2), dim = c(length(x), 2), 
                            dimnames = list(NULL, c('Name', 'Std.Dev.'))),
         `Cor.` = array(as.character(NA), dim = c(length(x), 1), 
                        dimnames = list(NULL, 'Cor.')))
    attr(out, 'covstruct') = 'cs'
    return(out)
  },
  # homcs = function(x, ...) {
  #   out = list(Group = array(NA, dim = c(1, 1), 
  #                            dimnames = list(NULL, 'Group')),
  #              `Std.Dev.` = array(rep(common_prefix(x), 2), dim = c(1, 2), 
  #                                 dimnames = list(NULL, c('Name', 'Std.Dev.'))),
  #              `Cor.` = array(as.character(NA), dim = c(1, 1), 
  #                             dimnames = list(NULL, 'Cor.')))
  #   attr(out, 'covstruct') = 'homcs'
  #   return(out)
  # },
  ar1  = function(x, ...) {
    out = list(Group = array(NA, dim = c(1, 1), 
                       dimnames = list(NULL, 'Group')),
         `Std.Dev.` = array(rep(common_prefix(x), 2), dim = c(1, 2), 
                            dimnames = list(NULL, c('Name', 'Std.Dev.'))),
         `Cor.` = array(as.character(NA), dim = c(1, 1), 
                        dimnames = list(NULL, 'Cor.')))
    attr(out, 'covstruct') = 'ar1'
    return(out)
  },
  hetar1 = function(x, ...) {
    out = list(Group = array(NA, dim = c(length(x), 1), 
                             dimnames = list(NULL, 'Group')),
               `Std.Dev.` = array(rep(x, 2), dim = c(length(x), 2), 
                                  dimnames = list(NULL, c('Name', 'Std.Dev.'))),
               `Cor.` = array(as.character(NA), dim = c(length(x), 1), 
                              dimnames = list(NULL, 'Cor.')))
    attr(out, 'covstruct') = 'hetar1'
    return(out)
  },
  gau = function(x, ...) {
    out = list(Group = array(NA, dim = c(1, 1), 
                             dimnames = list(NULL, 'Group')),
               `Std.Dev.` = array(rep(common_prefix(x), 2), dim = c(1, 2), 
                                  dimnames = list(NULL, c('Name', 'Std.Dev.'))),
               `Cor.` = array(as.character(NA), dim = c(1, 1), 
                              dimnames = list(NULL, 'Cor.')))
    attr(out, 'covstruct') = 'gau'
    return(out)
  },
  exp = function(x, ...) {
    out = list(Group = array(NA, dim = c(1, 1), 
                             dimnames = list(NULL, 'Group')),
               `Std.Dev.` = array(rep(common_prefix(x), 2), dim = c(1, 2), 
                                  dimnames = list(NULL, c('Name', 'Std.Dev.'))),
               `Cor.` = array(as.character(NA), dim = c(1, 1), 
                              dimnames = list(NULL, 'Cor.')))
    attr(out, 'covstruct') = 'exp'
    return(out)
  },
  ou = function(x, ...) {
    out = list(Group = array(NA, dim = c(1, 1), 
                             dimnames = list(NULL, 'Group')),
               `Std.Dev.` = array(rep(common_prefix(x), 2), dim = c(1, 2), 
                                  dimnames = list(NULL, c('Name', 'Std.Dev.'))),
               `Cor.` = array(as.character(NA), dim = c(1, 1), 
                              dimnames = list(NULL, 'Cor.')))
    attr(out, 'covstruct') = 'ou'
    return(out)
  },
  diag = function(x, ...) {
    out = list(Group = array(NA, dim = c(length(x), 1), 
                             dimnames = list(NULL, 'Group')),
               `Std.Dev.` = array(rep(x, 2), dim = c(length(x), 2), 
                                  dimnames = list(NULL, c('Name', 'Std.Dev.'))),
               `Cor.` = array(as.character(NA), dim = c(length(x), length(x)), 
                              dimnames = list(NULL, c('Cor.', rep('', length(x)-1)))))
    attr(out, 'covstruct') = 'diag'
    return(out)
  }
)

#' @noRd
#' @importFrom stats qlogis 
cor_convert_dispatch = list(
  us = function(x, ...){
    put_cor(C = x,  input_val = 'vec')
  }, 
  cs = function(x, ...){
    a = 1/(list(...)$n-1)
    qlogis((x+a)/(1+a))
  }, 
  # homcs = function(x, ...){
  #   a = 1/(list(...)$n - 1)
  #   qlogis((x+a)/(1+a))
  # }, 
  ar1 = function(x, ...){
    x/sqrt(1-x^2)
  }, 
  hetar1 = function(x, ...){
    x/sqrt(1-x^2)
  }, 
  gau = function(x, ...){
    -log(-log(x))/2
  }, 
  exp = function(x, ...){
    -log(-log(x))
  }, 
  ou = function(x, ...){
    log(-log(x))
  })

#' @noRd
struc_reterms = function(x, ...) {
  cs = attr(x, 'covstruct')
  if (is.null(cs)) stop("x has no 'covstruct' attribute")
  if (!(cs %in% names(cor_convert_dispatch))) stop(paste0(x, 'covstruct is not (yet) supported'))
  
  cs = match.arg(cs, choices = names(covstruct_dispatch))
  covstruct_dispatch[[cs]](x, ...)
}

#' @noRd
common_prefix = function(x) {
  x = x[!is.na(x)]
  if (length(x) <= 1) return(if (length(x) == 1) x else NA_character_)
  
  x = sort(x, method = 'radix')  # locale-independent byte order
  first = x[1]
  last  = x[length(x)]
  
  min_len = min(nchar(first), nchar(last))
  if (min_len == 0) return('')
  
  matched = substring(first, 1:min_len, 1:min_len) == substring(last, 1:min_len, 1:min_len)
  n_common = if (all(matched)) min_len else which.min(matched) - 1
  
  substr(first, 1, n_common)
}

#' @title Extract the possibly fixed residual dispersion parameter
#' 
#' @description
#' This is a helper function used for getting `sigma()` when it is not estimated
#' but fixed.
#' 
#' @param mod A glmmTMB model.
#' @param ... Ignored. 
#' 
#' @return Numeric. The value of the (estimated or specified) residual dispersion in the model.
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
#' # Extract residual SD:
#' extract_disp(mod)
#' 
#' @export
extract_disp = function(mod, ...){
  pars = mod$obj$env$parList()
  if ('betadisp' %in% names(pars)){
    disp = exp(pars$betadisp)
  } else if ('betad' %in% names(pars)) {
    disp = exp(pars$betad)
  } else {
    stop(simpleError('Models on the residual dispersion are not currently supported'))
  }
  return(disp)
}

#' @noRd
#' @importFrom stats family 
is_gaus_mod = function(emm){
  call = as.list(emm@model.info$call)
  family_arg = ifelse('family' %in% names(call), TRUE, FALSE)
  
  if (family_arg) {
    if (call$family$family == 'gaussian' & call$family$link == 'identity') {
      gaus = TRUE
    } else {
      gaus = FALSE
    }
  } else {
    gaus = TRUE
  } 
  return(gaus)
}

#' @noRd
is_reml = function(emm) {
  call = as.list(emm@model.info$call)
  reml_arg = ifelse('REML' %in% names(call), TRUE, FALSE)
  
  if (reml_arg) {
    if (call$REML == TRUE) {
      reml = TRUE
    } else {
      reml = FALSE
    }
  } else {
    reml = FALSE
  } 
  return(reml)
}

#' @noRd
#' @importFrom stats formula
#' @importFrom reformulas findbars
is_fe_mod = function(emm){
  is.null(findbars(formula(emm@model.info$call)))
}

#' @importFrom reformulas RHSForm nobars
#' @importFrom emmeans emmeans
#' @importFrom pbkrtest Lb_ddf
#' @importFrom glmmTMB glmmTMB getME dof_KR
#' @importFrom stats df.residual vcov formula
resolve_ddf = function(object, request = NULL){
   if (inherits(object, 'glmmTMB')){
    model = object
    fe_form = nobars(RHSForm(formula(object), as.form = T))
    emm = emmeans(object, specs = fe_form)
  } else if (inherits(object, 'emmGrid')){
    emm = object
  } else {
    stop(simpleError('object should be a glmmTMB model or emmGrid object constructed from one'))
  }
  
  fixed = is_fe_mod(emm)
  gaus = is_gaus_mod(emm)
  
  numddf = inherits(request, 'numeric')
  
  # Define default DDF method when none is specified
  if (is.null(request)){
    if (gaus){
      if (fixed){
          ddf = 'df.residual'
        } else {
          ddf = 'kenward-roger'
        }
      } else {
        ddf = 'asymptotic'
      }
  } else if (numddf){
    ddf = 'user-specified'
  } else {
    ddf = request
  }
  
  # Provide warnings and/or corrections for potentially inappropriate ddf specifications
  if (ddf == 'kenward-roger'){
    if (fixed & gaus){
        ddf = 'df.residual'
        message(simpleMessage('kenward-roger is only appropriate for models fit with REML, switching to ddf = "df.residual"'))
      } else if (!is_reml(emm)){
        ddf = 'asymptotic'
        message(simpleMessage('kenward-roger is only appropriate for models fit with REML, switching to ddf = "asymptotic"'))
      }
    }
  
  if (gaus & !fixed & ddf != 'kenward-roger'){
    message(simpleMessage('For gaussian mixed models, it is recomended to use ddf = "kenward-roger"'))
  }
  
  if (gaus & fixed & ddf == 'asymptotic'){
    message(simpleMessage('For gaussian fixed-effects models, it is recommended to use ddf = "df.residual"'))
  }
 
  if (inherits(object, 'emmGrid')){
    if (!numddf){
      if (ddf %in% c('df.residual', 'kenward-roger')){
        model = eval(emm@model.info$call)
      }
    }
  }
  
  dfargs = switch(ddf, 
                  asymptotic = list(), 
                  df.residual = list(object = eval(emm@model.info$call)), 
                  `kenward-roger` = list(V = vcov(model)$cond, 
                                       adjV = attr(dof_KR(model), 'vcov')), 
                  `user-specified` = list(df = request))
  dffun = switch(ddf, 
                 asymptotic = \(k, dfargs){Inf}, 
                 df.residual = \(k, dfargs){df.residual(dfargs$object)}, 
                 `kenward-roger` = \(k, dfargs){Lb_ddf(k, dfargs$V, dfargs$adjV)}, 
                 `user-specified` = \(k, dfargs){dfargs$df})
  out = list(dffun = dffun, dfargs = dfargs, ddf = ddf)
  return(out)
}

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

check_ddf = function(ddf, methods = c('df.residual', 'asymptotic', 'kenward-roger')){
  if (identical(ddf, NULL)){return(TRUE)}
  if (all(inherits(ddf, 'numeric'))){return(TRUE)}
  ifelse(ddf %in% methods, TRUE, stop(simpleError(paste('If non-NULL, ddf should be a numeric vector or one of:', paste(methods, collapse = ', ')))))
}
