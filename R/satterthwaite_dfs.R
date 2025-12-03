satterthwaite_aov = function(model = model, type = type){
  y_name<- names(model$modelInfo$respCol)
  dc <- dataClasses(model)
  TMBaov <- suppressPackageStartupMessages(car::Anova(model, type=type))
  
  # Pull the DFs associated with each term
  basic_aov_dfs <- base_aov_dfs(model)
  
  # Setup the final output
  if(row.names(TMBaov)[1] == "(Intercept)") {
    df_output <- basic_aov_dfs[basic_aov_dfs$vartype=="fixed", ]
  } else {
    df_output <- basic_aov_dfs[basic_aov_dfs$vartype=="fixed", ][-1, ]
  }
  
  df_output$vartype <- NULL
  df_output$denDf <- NA
  
  X = model.matrix(mod, 'X')
  batch = factor(attr(X, 'assign'))
  M = t(model.matrix(~ batch))
  M = M[rownames(M) != '(Intercept)', ]
  M = M/rowSums(M) 
  
  chisq <- as.vector(TMBaov$Chisq)
  nDF <- as.vector(TMBaov$Df)
  Fval <- chisq/nDF
  dDF <- glmmTMB::dof_satt(mod, M)
  Pval <- pf(Fval, nDF, dDF, lower.tail = FALSE)
  
  aod <- data.frame(numDF = nDF, denDF = dDF, Fvalue = round(Fval, 2), pvalue = round(Pval, 4))
  row.names(aod) <- df_output$terms
  class(aod) <- c("bdf", "containment", "data.frame")
  
  if (type == 3) {
    attr(aod, "heading") <-  paste("Analysis of Deviance Table (Type III F-tests)", "\n\nResponse: ", y_name)
  } else if (type == 2){
    attr(aod, "heading") <-  paste("Analysis of Deviance Table (Type II F-tests)", "\n\nResponse: ", y_name)
  }
  
  return(aod)
}