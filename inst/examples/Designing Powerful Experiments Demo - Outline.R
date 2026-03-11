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

#### Example 1: A two-sample t-test ####
# First, we'll demonstrate the "plug-in" approach and show that for simple
# cases it produces the same results as the tools in base R:  



#### Example 2: A Split-plot RCBD Experiment ####
# There was a study oats performed in 1931 at Rothamsted Station, the data from
# which was used by Frank Yates in an article entitled "Complex Experiments".
# We are going to perform a power analysis on a proposed recreation of that study:
# - Factorial treatment structure, with 3 cultivars (GoldenRain, Marvellous and 
#   Victory) and 4 nitrogen application rates (0, 20, 40, 60 kg/ha)
# - Cultivars are randomly allocated to main-plots within each of six blocks, and 
#   nitrogen levels randomized to subplots within each mainplot (i.e., cultivar)
#   within each block



#### Example 3: RCBD for Binomial Data ####
# This example is inspired by a study I worked on comparing methods for 
# controlling thrips. In this simplified version, we'll have three treatments 
# applied in each of two ways. And since that's only 6 treatment combinations,
# I'm just going to construct my treatement dictionary in R:



#### Example 4: Integration with edibble ####
# The edibble package has lots of nice tools for experimental design, 
# which can be incorporated into the power analysis. So instead of using
# expand.grid to create the base for our fake data set, we could use edibble:







