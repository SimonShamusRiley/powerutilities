# powerutilities
<!-- badges: start --> 
<!-- badges: end -->

This package provides some tools for power analysis when using what [Piepho et. al (2022)](http://www.doi.org/10.1017/S0021859622000466) 
called the "plug-in approach" to power analysis, as described in chapter 
*Precision, Power, Sample Size and Planning* in 
[Generalized Linear Mixed Models: Modern Concepts, Methods and Applications](https://www.taylorfrancis.com/chapters/mono/10.1201/9780429092060-24/precision-power-sample-size-planning-walter-stroup-marina-ptukhina-julie-garai?context=ubx&refId=7081f25c-0d97-41eb-b09c-a81ef861e2ef).

The basic idea is to encode information about sample size and expected
treatments means (and differences) for the proposed experiment in a fake data
set, to then specify a model encoding information about the design of the
proposed experiment, and finally to "fit" the model with plug-in values for the
random effects and dispersion parameters (based on pilot studies, previously
published research or estimated guesswork). The result is that the mixed
modelling software does all the requisite work of calculating standard errors
upon which power can be derived. By manipulating the fake data set, model, or
variance parameter values, it becomes possible to explore how differences in
sample size, effect size, experimental design, or random effects will impact
power.

The vignettes provided with the package serve a number of aims:

- To provide a fuller, more detailed description of this approach to power 
analysis

- Present a variety of practical examples demonstrating the various ways this 
approach can serve the concrete needs of practicing researchers (especially 
those in the agricultural sciences)

- Provide a vehicle for introducing and discussing a number of complexities related
to the use of mixed models, and generalized linear mixed models in particular, which I
think are not as fully or widely appreciated as they should be. 

**Please note that this effort remains very much a work in progress**. Any issues,
concerns, questions, bugs, feature requests, etc., can be posted on the 
[Issues page](https://github.com/SimonShamusRiley/powerutilities/issues).

## Installation

You can install the development version of powerutilities from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github('SimonShamusRiley/powerutilities', 
                         build_vignettes = T)
```



