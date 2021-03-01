# CoMoMo

Welcome to CoMoMo: CoMoMo combines multiple mortality forecasts using different model combinations.

CoMoMo is an R package for combining the mortality rate forecasts from the Generalized Age Period cohort models. 

## Prerequisites

The codes use the cross-validation approach implemented in the `StMoMo` package. It can be installed using the devtools package:

``` r
devtools::install_github("amvillegas/StMoMo", ref = "GroupLasso")
```

## Installation

The easiest way to install CoMoMo is using the devtools package:

``` r
devtools::install_github("kessysalvatory/CoMoMo")
```

## Example

This is a basic example which shows you how to combine forecasts from multiple mortality models:

``` r
library(demography)
library(StMoMo)
library(CoMoMo)
library(tibble)
library(ggplot2)
library(patchwork)

# Getting some data

MorData <- hmd.mx(country = 'GBRTENW',username = "your email", password = "your password")
DataStMoMo <- StMoMoData(MorData, "male")

agesFit <- 50:89
yearsFit <- 1960:1990
nAg <- length(agesFit); nYr <- length(yearsFit)

# Choosing which models to fit
# Define the models

LC <- lc()
APC <- apc()
CBD <- cbd(link = "log")
M7 <- m7(link = "log")
RH <- rh(approxConst = TRUE)
f2 <- function(x, ages) mean(ages) - x
f3 <- function(x, ages) {ifelse((mean(ages) - x)<0, 0, (mean(ages) - x))}
constPlat <- function(ax, bx, kt, b0x, gc, wxt, ages){
  nYears <- dim(wxt)[2]
  x <- ages
  t <- 1:nYears
  c <- (1 - tail(ages, 1)):(nYears - ages[1])
  xbar <- mean(x)
  phiReg <- lm(gc ~ 1 + c + I(c ^ 2), na.action = na.omit)
  phi <- coef(phiReg)
  gc <- gc - phi[1] - phi[2]*c-phi[3]*c^2
  kt[2,] <- kt[2,] + 2*phi[3]*t
  kt[1,] <- kt[1,]+phi[2]*t+phi[3]*(t^2 - 2*xbar*t)
  ax <- ax+phi[1] - phi[2]*x+phi[3]*x^2
  ci <- rowMeans(kt, na.rm = TRUE)
  ax <- ax+ci[1]+ci[2]*(xbar-x)
  kt[1, ] <- kt[1, ] - ci[1]
  kt[2, ] <- kt[2, ] - ci[2]
  list(ax = ax,bx = bx,kt = kt,b0x = b0x,gc = gc)
}
PLAT <- StMoMo(link = "log", staticAgeFun = TRUE, periodAgeFun = c("1", f2),
               cohortAgeFun = "1", constFun = constPlat)
               
# model list 

models <- list("LC" = LC, "RH" = RH, "APC" = APC, "CBD" = CBD, "M7" = M7, "PLAT" = PLAT)

# Generate the metadata for stacked regression ensemble
# h is the forecasting horizon
# First, we train the models via cv to produce the cross-validated forecasts
# Combine the forecasts and observed rates to form a metadata.

metaData <- stackMetadata(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 15)

# Estimating the weights using different learners

# When normalize = TRUE all the weights sum to a unit

stack_ridge_weight <- stack(metaData, metalearner = "Ridge", normalize = TRUE)
stack_lasso_weight <- stack(metaData, metalearner = "Lasso", normalize = TRUE)
stack_nnls_weight <- stack(metaData, metalearner = "nnls", normalize = TRUE)
stack_linear_weight <- stack(metaData, metalearner = "Linear", normalize = TRUE)
stack_elastic_weight <- stack(metaData, metalearner = "Elastic", normalize = TRUE)

# When normalize = False all the weights may not sum to a unit

stack_ridge_weight0 <- stack(metaData, metalearner = "Ridge", normalize = FALSE)
stack_lasso_weight0 <- stack(metaData, metalearner = "Lasso", normalize = FALSE)
stack_nnls_weight0 <- stack(metaData, metalearner = "nnls", normalize = FALSE)
stack_linear_weight0 <- stack(metaData, metalearner = "Linear", normalize = FALSE)
stack_elastic_weight0 <- stack(metaData, metalearner = "Elastic", normalize = FALSE)

# plot the weights 

stack_ridge_plot <- plot(stack_ridge_weight)
stack_lasso_plot <- plot(stack_lasso_weight)
stack_nnls_plot <- plot(stack_nnls_weight)
stack_linear_plot <- plot(stack_linear_weight)
stack_elastic_plot <- plot(stack_elastic_weight)

# use the Bayesian model averaging (bma) to estimate the weights
# use single validation set appproach to estimate the weights 

bma_weight_val <- bma(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 15, method = "sv")

# Cross-validation is used to calculate the cross-validation mean squared errors 

bma_weight_cv <- bma(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 15, method = "cv")


# use the Model Confidence Set (mcs) to choose the models to combine
# use single validation set appproach to estimate the weights 

mcs_weight_val <- mcs(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 15, method = "sv")

# Models are chosen via cross-validation
# Models are chosen for each horizon

mcs_weight_cv <- mcs(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 15,  method = "cv")

# Fitting the models

modelFits <- fitCoMoMo(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit)

# combining the fitted mortality models and different weights 

modcom <- CoMoMo(modelFits, weight = stack_ridge_weight)

# forecast the mortality rates

final <- forecast(modcom, h = 15)
```

## Questions 

Please feel free to open an issue with any questions you may have. You can contact me at s.kessy@unsw.edu.au
