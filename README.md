# CoMoMo

Welcome to CoMoMo: CoMoMo combines multiple mortality forecasts using different model combinations.

CoMoMo is an R package for combining the mortality rate forecasts from the Generalized Age Period cohort models. 

## Prerequisites

The codes use the cross-validation approach implemented in the `StMoMo` package. It can be installed using the devtools package:

``` r
devtools::install_github("amvillegas/StMoMo", ref = "GroupLasso")
```

## Installation

If you have R 4.0.2 or later and use Windows or Mac, the easiest way to install CoMoMo is using the devtools package:

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
library(gganimate)

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
PLAT <- plat()
            
# model list 

models <- list("LC" = LC, "RH" = RH, "APC" = APC, "CBD" = CBD, "M7" = M7, "PLAT" = PLAT)

# Generate the metadata for stacked regression ensemble
# h is the forecasting horizon
# First, we train the models via cv to produce the cross-validated forecasts
# Combine the forecasts and observed rates to form a metadata.

metaData <- metadata(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 15)

# Estimating the weights using different learners

# When normalize = TRUE all the weights sum to a unit and dynamic = TRUE to get horizon-specific weights 

stack_ridge_weight <- stack(metaData, metalearner = "Ridge", normalize = TRUE, dynamic = TRUE)
stack_lasso_weight <- stack(metaData, metalearner = "Lasso", normalize = TRUE, dynamic = TRUE)
stack_nnls_weight <- stack(metaData, metalearner = "nnls", normalize = TRUE, dynamic = TRUE)
stack_linear_weight <- stack(metaData, metalearner = "Linear", normalize = TRUE, dynamic = TRUE)
stack_elastic_weight <- stack(metaData, metalearner = "Elastic", normalize = TRUE, dynamic = TRUE)


# When normalize = False all the weights may not sum to a unit

stack_ridge_weight0 <- stack(metaData, metalearner = "Ridge", normalize = FALSE)
stack_lasso_weight0 <- stack(metaData, metalearner = "Lasso", normalize = FALSE)
stack_nnls_weight0 <- stack(metaData, metalearner = "nnls", normalize = FALSE)
stack_linear_weight0 <- stack(metaData, metalearner = "Linear", normalize = FALSE)
stack_elastic_weight0 <- stack(metaData, metalearner = "Elastic", normalize = FALSE)

# When dynamic = FALSE, we get a single set of constant weights for all the forecasting horizons. 

dstack_ridge_weight <- stack(metaData, metalearner = "Ridge", normalize = TRUE, dynamic = FALSE)
dstack_lasso_weight <- stack(metaData, metalearner = "Lasso", normalize = TRUE, dynamic = FALSE)
dstack_nnls_weight <- stack(metaData, metalearner = "nnls", normalize = TRUE, dynamic = FALSE)
dstack_linear_weight <- stack(metaData, metalearner = "Linear", normalize = TRUE, dynamic = FALSE)
dstack_elastic_weight <- stack(metaData, metalearner = "Elastic", normalize = TRUE, dynamic = FALSE)


# plot the weights 

stack_ridge_plot <- plot(stack_ridge_weight)
stack_lasso_plot <- plot(stack_lasso_weight)
stack_nnls_plot <- plot(stack_nnls_weight)
stack_linear_plot <- plot(stack_linear_weight)
stack_elastic_plot <- plot(stack_elastic_weight)

# use the Bayesian model averaging (bma) to estimate the weights
# use single validation set appproach to estimate the weights 

bma_weight_val <- bma(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 15, method = "sv")

# plot the weights 

bma_weight_val_plot <- plot(bma_weight_val)

# Cross-validation is used to calculate the cross-validation mean squared errors 

bma_weight_cv <- bma(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 15, method = "cv")

# plot the weights 

bma_weight_cv_plot <- plot(bma_weight_cv)

# use the Model Confidence Set (mcs) to choose the models to combine
# use single validation set appproach to estimate the weights 

mcs_weight_val <- mcs(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 15, method = "sv")

# plot the weights 

mcs_weight_val_plot <- plot(mcs_weight_val)

# Models are chosen via cross-validation
# Models are chosen for each horizon

mcs_weight_cv <- mcs(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 15,  method = "cv")

# plot the weights 

mcs_weight_cv_plot <- plot(mcs_weight_cv)

# Estimates the weights using the frequentist approach 
# Estimates the weight at horizon h as the ratio of cvmse(h)/sum(cvmse(h))

fweights <- frequentist(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 15)

# Fitting the models

modelFits <- fitCoMoMo(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit)

# combining the fitted mortality models and different weights 

modcom <- CoMoMo(modelFits, weight = stack_ridge_weight)

# forecast the mortality rates

final <- forecast(modcom, h = 15)
```

## Questions 

Please feel free to open an issue with any questions you may have. You can contact me at s.kessy@unsw.edu.au
