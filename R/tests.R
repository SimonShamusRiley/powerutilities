# # Test
# library(tidyverse)
# library(glmmTMB)
# library(emmeans)
# library(reformulas)
# library(pbkrtest)
# library(retrodesign)
# 
# dat = expand.grid(A = factor(1:2), 
#                   B = factor(1:2), 
#                   Rep = factor(1:8)) |> 
#   mutate(Y = ifelse(A == '1', 1, 3), 
#          Y = ifelse(B == '1', Y+0, Y+1))
# 
# #### Linear Models ####
# lm = glmmTMB(Y~A*B, data = dat, family = gaussian(link = 'identity'), 
#               start = list(betadisp = log(5)), 
#               map = list(betadisp = factor(NA)), 
#               REML = T)
# theta_finder(~ A*B, data = dat)
# lm = set_glmm(Y~A*B, data = dat, disp = 5)
# extract_disp(lm)
# 
# #### Linear Mixed Models ####
# theta_finder(~ A*B + (1|Rep/A), data = dat)
# lmm = set_glmm(Y~A*B+ (1|Rep/A), data = dat, re_terms = c(1, 2), disp = 5)
# VarCorr(lmm, sigma = extract_disp(lmm))
# 
# #### Generalized Linear Models ####
# glm = glmmTMB(Y~A*B, data = dat, family = poisson(link = 'log'), 
#               #start = list(betadisp = log(2)), 
#               #map = list(betadisp = factor(NA)), 
#               REML = F)
# theta_finder(~A*B, data = dat)
# glm = set_glmm(Y~A*B, data = dat, family = poisson(link = 'log'), 
#                REML = F)
# glm
# #### Generalized Linear Models ####
# theta_finder(~A*B + (1|Rep), data = dat)
# glmm = set_glmm(Y~A*B + (1|Rep), data = dat, family = poisson(link = 'log'), 
#                 re_terms = 2, REML = F)
# VarCorr(glmm)
# 
# power_ftest(lm, ddf = NULL)
# power_ftest(lm, ddf = 4)
# power_ftest(lm, ddf = 'df.residual')
# power_ftest(lm, ddf = 'asymptotic')
# power_ftest(lm, ddf = 'kenward-roger')
# 
# power_ftest(lmm, ddf = NULL)            
# power_ftest(lmm, ddf = 4)
# power_ftest(lmm, ddf = 'df.residual')
# power_ftest(lmm, ddf = 'asymptotic')
# power_ftest(lmm, ddf = 'kenward-roger') 
# 
# power_ftest(glm, ddf = NULL)            
# power_ftest(glm, ddf = 4)
# power_ftest(glm, ddf = 'df.residual')
# power_ftest(glm, ddf = 'asymptotic')
# power_ftest(glm, ddf = 'kenward-roger') 
# 
# 
# power_ftest(glmm, ddf = NULL)            
# power_ftest(glmm, ddf = 4)
# power_ftest(glmm, ddf = 'df.residual')
# power_ftest(glmm, ddf = 'asymptotic')
# power_ftest(glmm, ddf = 'kenward-roger') 
# 
# contr_list = list('A2 - A1' = c(-0.5, 0.5, -0.5, 0.5))
# lm_emm = emmeans(lm, ~ A:B)
# glm_emm = emmeans(glm, ~A:B)
# lmm_emm = emmeans(lmm, ~A:B)
# glmm_emm = emmeans(glmm, ~A:B)
# 
# power_contrast(lm_emm, contr_list, ddf = NULL)
# power_contrast(lm_emm, contr_list, ddf = 2)
# power_contrast(lm_emm, contr_list, ddf = 'df.residual')
# power_contrast(lm_emm, contr_list, ddf = 'asymptotic')
# power_contrast(lm_emm, contr_list, ddf = 'kenward-roger')
# 
# power_contrast(lmm_emm, contr_list, ddf = NULL)            
# power_contrast(lmm_emm, contr_list, ddf = 4)
# power_contrast(lmm_emm, contr_list, ddf = 'df.residual')
# power_contrast(lmm_emm, contr_list, ddf = 'asymptotic')
# power_contrast(lmm_emm, contr_list, ddf = 'kenward-roger') 
# 
# power_contrast(glm_emm, contr_list, ddf = NULL)            
# power_contrast(glm_emm, contr_list, ddf = 4)
# power_contrast(glm_emm, contr_list, ddf = 'df.residual')
# power_contrast(glm_emm, contr_list, ddf = 'asymptotic')
# power_contrast(glm_emm, contr_list, ddf = 'kenward-roger') 
# 
# power_contrast(glmm_emm, contr_list, ddf = NULL)            
# power_contrast(glmm_emm, contr_list, ddf = 4)
# power_contrast(glmm_emm, contr_list, ddf = 'df.residual')
# power_contrast(glmm_emm, contr_list, ddf = 'asymptotic')
# power_contrast(glmm_emm, contr_list, ddf = 'kenward-roger') 
