#==============================================================================#
# Designing Powerful Experiments: Practical Approach to Power Analysis in R    #
#==============================================================================#

#### Housekeeping ####
# Install required packages
install.packages('tidyverse')
install.packages('devtools')
install.packages('readxl')
install.packages('writexl')
install.packages('glmmTMB')
install.packages('emmeans')
install.packages('edibble')

devtools::install_github('SimonShamusRiley/powerutilities', 
                         dependencies = TRUE, build_vignettes = TRUE)

# Load required packages
library(tidyverse)
library(glmmTMB)
library(emmeans)
library(powerutilities)
library(edibble)
library(readxl)
library(writexl)

# To learn more about the powerutilities package
vignette(package = 'powerutilities')
vignette(topic = 'intro-powerutilities')

# Set working directory
setwd("C:/Users/au802896/OneDrive - Aarhus universitet/Desktop/stat power demo")

#### Example 1: A two-sample t-test ####
# First, we'll demonstrate the "plug-in" approach and show that for simple
# cases it produces the same results as the tools in base R:  
?power.t.test()
power.t.test(n = 8, delta = 3, sd = 1.5, sig.level = 0.05, power = NULL)
power.t.test(n = NULL, delta = 3, sd = 1.5, sig.level = 0.05, power = 0.96)

# In the plug-in approach, we encode the *same* information in a fake data
# set and model, then let the model fitting software do (most) of the 
# calculations for us

# Create data set
ex1_dat = expand.grid(Rep = factor(1:8), 
                  Trt = c('A', 'B')) |> 
  mutate(Y = ifelse(Trt == 'A', 0, 3))
print(ex1_dat)

# "fit" model while fixing residual SD at 1.5:
ex1_mod = glmmTMB(Y ~ Trt, data = ex1_dat,
                  # set residual sd to zero...
                  dispformula = ~ 0,
                  # ...then define what value to use for "zero"
                  control = glmmTMBControl(zerodisp_val = log(1.5)))
# It worked!
sigma(ex1_mod)

# Now we can use tools from the powerutilties package
# First, the power of the F-test. Set the denominator degrees of freedom
# method to "residual" - we'll return to this later:
power_ftest(ex1_mod, ddf = 'residual')

# Then, the power for specific contrasts (to learn more about specifying
# contrasts in R, see the vignettes for powerutilities and the emmeans 
# packages)
# estimate marginal means:
ex1_emm = emmeans(ex1_mod, ~ Trt)
print(ex1_emm)

# define contrast(s):
ex1_contr_list = list('A - B' = c(1, -1))

# calculate power:
power_contrast(ex1_emm, ex1_contr_list, ddf = 14)

#### Example 2: A Split-plot RCBD Experiment ####
# There was a study oats performed in 1931 at Rothamsted Station, the data from
# which was used by Frank Yates in an article entitled "Complex Experiments".
# We are going to perform a power analysis on a proposed recreation of that study:
# - Factorial treatment structure, with 3 cultivars (GoldenRain, Marvellous and 
#   Victory) and 4 nitrogen application rates (0, 20, 40, 60 kg/ha)
# - Cultivars are randomly allocated to main-plots within each of six blocks, and 
#   nitrogen levels randomized to subplots within each mainplot (i.e., cultivar)
#   within each block
ex2_dat = expand.grid(Block = 1:6, 
                      Cultivar = c('GoldenRain', 'Marvellous', 'Victory'), 
                      Nitrogen = c(0, 20, 40, 60)) |> 
  arrange(Cultivar, Nitrogen)
print(ex2_dat)

# Note that unless we have row and column values for each plot, we cannot see
# the experimental design in the data, meaning there is no way for R to know
# how to analyse the data if we don't tell it via the model formula, just like
# no other person would be able to know how to analyse the data correctly if 
# we don't give them details about the experimental setup
print(ex2_dat)

# Now, we need to define the treatment means for each group. This can be quite
# tedious in R, so we're actually going to do it in Excel and then simply 
# read the values in. 
# First, create a table with each unique combination of treatments and save it:
ex2_trt_lvls = ex2_dat |> 
  distinct(Cultivar, Nitrogen) |> 
  arrange(Cultivar, Nitrogen)
ex2_trt_lvls

write_xlsx(ex2_trt_lvls, path = 'ex2_trt_dict.xlsx')

# Now open the ex2_trt_dict.xlsx sheet in excel, and enter the following
# information yield information:
# Victory and GoldenRain have a VarEff value of 1400, while Marvelleous has 1500,
# Victory and Marvellous have NEff of (0, 150, 300, 450), while GoldenRain has
# NEff values of (0, 200, 400, 600). Yield for each treatment combination is just
# the sum of VarEff and NEff

# Read the treatment dictionary back in and merge it
# to the fake data set
ex2_trt_dict = read_xlsx(path = 'ex2_trt_dict.xlsx') 
ex2_trt_dict

ex2_dat = left_join(ex2_dat, ex2_trt_dict)|> 
  mutate(across(Block:Nitrogen, factor))
print(ex2_dat)

# "Fit model". Here, in addition to setting the residual variance, we 
# also need to set the random effect variances, which we do using 
# a combination of the "start" and "map" arguments
ex2_mod = glmmTMB(Yield ~ Cultivar*Nitrogen + (1|Block/Cultivar), 
                  data = ex2_dat, REML = T, 
                  family = gaussian(link = 'identity'),
                  dispformula = ~0, 
                  control = glmmTMBControl(zerodisp_val = log(200)), 
                  start = list(theta = c(log(200), log(250))), 
                  map = list(theta = factor(c(NA, NA))))
VarCorr(ex2_mod)
sigma(ex2_mod)

# Now we can assess power. To start, let's use "containment" degrees of freedom,
# which is appropriate for this kind of design:
power_ftest(ex2_mod, ddf = 'containment')

# Here, we can play with different designs

# Now let's consider different DenDF calculations. "Residual" uses the same
# DFs for all tests, basically ignoring the experimental design
power_ftest(ex2_mod, ddf = 'residual')

# "Asymptotic" both ignores the design and exaggerates and p-values and power
# should be taken as lower/upper limits (respectively)
power_ftest(ex2_mod, ddf = 'asymptotic')

# Now let's look at some specific contrasts. Imagine that Marvellous and 
# Golden Rain are new varieties which we want to compare against the current
# standard, Victory, at each N rate. First, we estimate marginal means using
# Nitrogen as a grouping variable:
ex2_emm = emmeans(ex2_mod, ~ Cultivar | Nitrogen)
print(ex2_emm)

# Note that glmmTMB will use infinite denDFs (i.e., asymptotic tests) by 
# default

# Now define contrasts among cultivars (can ignore nitrogen: since its a 
# grouping variable contrasts will be applied at each N rate)
ex2_contr_list = list('GoldenRain - Victory' = c(1, 0, -1), 
                      'Marvellous - Victory' = c(0, 1, -1))

# Calculate power for each contrast
power_contrast(ex2_emm, ex2_contr_list, ddf = 45)

#### Example 3: RCBD for Binomial Data ####
# This example is inspired by a study I worked on comparing methods for 
# controlling thrips. In this simplified version, we'll have three treatments 
# applied in each of two ways. And since that's only 6 treatment combinations,
# I'm just going to construct my treatement dictionary in R:
ex3_trt_dict = expand.grid(Insecticide = c('Water', 'Spinosad', 'Sesame Oil'), 
                           Application = c('Curative', 'Prophylactic')) |> 
  mutate(Mortality = c(0.05, 0.20, 0.15, 
                       0.05, 0.35, 0.25))

print(ex3_trt_dict)

# Now our study is going to be arranged in blocks, but imagine we are limited
# in our number of thrips, and we want to see what happens if we have more blocks
# with fewer thrips per plant or fewer blocks with more thrips per plant
ex3a_des = data.frame(Block = 1:8, 
                      Ttl_Thrips = 5)

ex3b_des = data.frame(Block = 1:4, 
                      Ttl_Thrips = 8)

ex3a_dat = merge(ex3_trt_dict, ex3a_des) |> 
  mutate(across(c(Block, Insecticide, Application), factor))
print(ex3a_dat)

ex3b_dat = merge(ex3_trt_dict, ex3b_des) |> 
  mutate(across(c(Block, Insecticide, Application), factor))
print(ex3b_dat)

# Now we "fit" the model. Note that since the binomial distribution
# has only one parameter, we do not need to use the "dispformula = ~0",
# as there is no dispersion parameter to estimate/set. However, we
# do have the added complication that the random effects variances are
# on the link scale, so it is not intuitive to understand how much variability
# is actually implied by a particular SD.
ex3a1_mod = glmmTMB(Mortality ~ Insecticide*Application + (1|Block), data = ex3a_dat, 
                   family = binomial(link = 'logit'), weights = Ttl_Thrips,
                   start = list(theta = log(.1)), 
                   map = list(theta = factor(NA)))

ex3b1_mod = glmmTMB(Mortality ~ Insecticide*Application + (1|Block), data = ex3b_dat, 
                   family = binomial(link = 'logit'), weights = Ttl_Thrips,
                   start = list(theta = log(.1)), 
                   map = list(theta = factor(NA)))

ex3a2_mod = glmmTMB(Mortality ~ Insecticide*Application + (1|Block), data = ex3a_dat, 
                   family = binomial(link = 'logit'), weights = Ttl_Thrips,
                   start = list(theta = log(1)), 
                   map = list(theta = factor(NA)))

ex3b2_mod = glmmTMB(Mortality ~ Insecticide*Application + (1|Block), data = ex3b_dat, 
                   family = binomial(link = 'logit'), weights = Ttl_Thrips,
                   start = list(theta = log(1)), 
                   map = list(theta = factor(NA)))

# First, let's compare designs using asymptotic test, which is the default for
# GLMs and GLMMs in R. We can see that the design with more blocks and fewer
# thrips per plant has higher power across the board
power_ftest(ex3a1_mod, ddf = 'asymptotic')
power_ftest(ex3b1_mod, ddf = 'asymptotic')

# We can also see the effect of using different DF calculations can differently
# impact the power: for insecticide, design "B" has about 9 percentage points
# less power when using Inf DF, while design "B" has 14 percentage points less
# power when using residual DF:
power_ftest(ex3a1_mod, ddf = 'residual')
power_ftest(ex3b1_mod, ddf = 'residual')

# We can also examine the effect of block variance
power_ftest(ex3a1_mod, ddf = 'residual')
power_ftest(ex3a2_mod, ddf = 'residual')

# Now let's consider the power of some contrasts. But first, note something
# odd about the emmeans? The estimated probabilities are not what we set them to
# be. This is because what we specified was a marginal probability and what we're
# estimating is a conditional probability. For the normal distribution, 
# the marginal and conditional means are the same, but not so for other 
# distributions, and the magnitude of the difference depends in part on the
# random effects variances
(ex3a1_emm = emmeans(ex3a1_mod, ~ Insecticide*Application, type = 'response'))

# Now lets consider two contrasts, where the true difference is 0.1
ex3_contr_list = list('Cur Ses - Cur Water' = c(-1, 0, 1, 0, 0, 0), 
                      'Pro Spin - Pro Ses' = c(0, 0, 0, 0, 1, -1) )

# We can see that even though the *difference* is the same in these two contrasts
# what is being tested is actually the *odds ratio* (because of the logit link
# function, for log link functions it is *ratio* that is tested). 
power_contrast(ex3a1_emm, ex3_contr_list, ddf = Inf)

#### Example 4: Integration with edibble ####
# The edibble package has lots of nice tools for experimental design, 
# which can be incorporated into the power analysis. So instead of using
# expand.grid to create the base for our fake data set, we could use edibble:
exp_des = design(name = 'Oat Study') |> 
  set_units(Block = 6, 
            MainPlot = nested_in(Block, 3), 
            SubPlot = nested_in(MainPlot, 4)) |> 
  set_trts(Cultivar = c('GoldenRain', 'Marvellous', 'Victory'),
           Nitrogen = c(0, 20, 40, 60)) |> 
  allot_trts(Cultivar ~ MainPlot, 
             Nitrogen ~ SubPlot) |> 
  assign_trts('random')

plot(exp_des)
exp_des_tab = serve_table(exp_des)
print(exp_des_tab)







