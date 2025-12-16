# powerutilities
<!-- badges: start --> 
<!-- badges: end -->

The powerutilities package provides some tools for power analysis when using
what [Piepho et. al (2022)](http://www.doi.org/10.1017/S0021859622000466) called
the "plug-in approach" to power analysis, or what might also be called the
Littel-Stroup method, as described in the chapter *Precision, Power, Sample Size
and Planning* in [Generalized Linear Mixed Models: Modern Concepts, Methods and
Applications](https://www.taylorfrancis.com/chapters/mono/10.1201/9780429092060-24/precision-power-sample-size-planning-walter-stroup-marina-ptukhina-julie-garai?context=ubx&refId=7081f25c-0d97-41eb-b09c-a81ef861e2ef).

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

In the course of such investigations, it's important that an appropriate choice
of (approximate) denominator degrees of freedom be used, but unfortunately the
`glmmTMB` package - among the only R packages which enables fitting of GLMMs will
fixing parameter values - defaults to the use of infinite degrees of freedom 
(Wald tests; asymptotic results). 

## Installation

You can install the development version of powerutilities from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("SimonShamusRiley/powerutilities")
```

Additional descriptions and worked examples will be forthcoming.