
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
#' @export
extract_disp = function(mod, ...){
  pars = mod$obj$env$parList()
  if ('betadisp' %in% names(pars)){
    exp(pars$betadisp)
  } else if ('betad' %in% names(pars)) {
    exp(pars$betad)
  } else {
    stop(simpleError('Models on the residual dispersion are not currently supported'))
  }
}


