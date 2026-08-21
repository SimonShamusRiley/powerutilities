test_that("set_glmm produces glmmTMB object", {
  dat = expand.grid(A = factor(1:2), 
                    B = factor(1:3),
                    Rep = factor(1:8)) 
  dat$Y = with(dat, as.numeric(A) + as.numeric(B))

  lm_mod = set_glmm(Y~A*B, data = dat, disp = 2)
  
  expect_equal(class(lm_mod), 'glmmTMB')
})

test_that("LM residual std. dev. is correctly set", {
  dat = expand.grid(A = factor(1:2), 
                    B = factor(1:3),
                    Rep = factor(1:8)) 
  dat$Y = with(dat, as.numeric(A) + as.numeric(B))
  disp_val = 2
  lm_mod = set_glmm(Y~A*B, data = dat, disp = disp_val)
  disp_set = exp(getME(lm_mod, 'beta')[['betadisp']])
  expect_equal(disp_val, disp_set)
})

test_that("LMM RE std. dev. correctly set", {
  dat = expand.grid(A = factor(1:2), 
                    B = factor(1:3),
                    Rep = factor(1:8)) 
  dat$Y = with(dat, as.numeric(A) + as.numeric(B))
  re_vals = c(1, 2)
  lmm_mod = set_glmm(Y~A*B + (1|Rep/A), data = dat, re_terms = re_vals, disp = 3)
  re_set = exp(getME(lmm_mod, 'theta'))
  expect_equal(re_vals, re_set)
})

test_that('LMM RE US cov correctly set', {
  dat = expand.grid(A = factor(1:2),
                    Rep = factor(1:8)) 
  dat$X = runif(nrow(dat))
  dat$Y = with(dat, as.numeric(A) + 2*X)
  re_vals = c(.1, .2, .3, .4)
  
  lmm_mod = set_glmm(Y~A*X + (1|Rep) + (X + 1|Rep:A), data = dat,
                     re_terms = re_vals, disp = 3)
  
  vc = VarCorr(lmm_mod)$cond
  re_set = unlist(sapply(vc, attr, 'stddev'), use.names = F)
  cor_set = as.vector(attr(vc$`Rep:A`, 'correlation'))
  
  expect_equal(re_vals[1:3], re_set)
  expect_equal(re_vals[4], unique(cor_set[2:3]))
})

test_that('LMM RE AR1 cov correctly set', {
  dat = expand.grid(A = factor(1:2),
                    Rep = factor(1:8)) 
  dat$X = runif(nrow(dat))
  dat$Y = with(dat, as.numeric(A) + 2*X)
  re_vals = c(.1, .2, .3, .4)
  
  lmm_mod = set_glmm(Y~A*X + (1|Rep) + (X + 1|Rep:A), data = dat,
                     re_terms = re_vals, disp = 3)
  
  vc = VarCorr(lmm_mod)$cond
  re_set = unlist(sapply(vc, attr, 'stddev'), use.names = F)
  cor_set = as.vector(attr(vc$`Rep:A`, 'correlation'))
  
  expect_equal(re_vals[1:3], re_set)
  expect_equal(re_vals[4], unique(cor_set[2:3]))
})


