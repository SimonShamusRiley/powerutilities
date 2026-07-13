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

#'@export
print.retermslist = function(x) {
  pad_columns <- function(mat_list, cols) {
    widths <- sapply(cols, function(col) {
      max(sapply(mat_list, function(m) max(nchar(m[, col]), na.rm = TRUE)))
    })
    names(widths) <- cols
    
    lapply(mat_list, function(m) {
      for (col in cols) {
        w <- widths[[col]]
        x <- m[, col]
        is_na <- is.na(x)
        x[!is_na] <- formatC(x[!is_na], width = -w)
        x[is_na] <- strrep(" ", w)
        m[, col] <- x
      }
      m
    })
  }
  
  comb_terms = lapply(x, \(x) cbind(x$Group, x$`Std.Dev.`, x$`Cor.`)) |> 
    pad_columns(cols = c('Group', 'Name'))
  
  for (n in 1:length(comb_terms)){
    colnames(comb_terms[[n]])[4] = paste0('Cor. (', attr(x[[n]], 'covstruct'), ')')
    prmatrix(comb_terms[[n]], quote = F, na.print = '', rowlab = rep('', nrow(comb_terms[[n]])))
    cat('\n')
  }
}
