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
devtools::install_github(""kessysalvatory/CoMoMo"")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(demography)
library(StMoMo)
library(CoMoMo)
library(tibble)
library(ggplot2)
library(patchwork)

MorData <- hmd.mx(country = 'GBRTENW',username = "salvatory@aims.ac.tz", password = "Salva=0606")
DataStMoMo <- StMoMoData(MorData, "male")

agesFit <- 50:89
yearsFit <- 1960:1990
nAg <- length(agesFit); nYr <- length(yearsFit)

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

models <- list(LC, RH, APC, M7, PLAT, CBD)

# model names 

modelNames <- c("LC","RH", "APC", "CBD", "M7", "PLAT")

# use the function model() to assign names to the model list

modlist <- model(models, modelNames)

# Compute the weights
# use the stacked regression ensemble (stack) with nnls as a metalearner 

stacked <- stack(modlist, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 5, metalearner = "nnls")
plot(stacked)

# use the Bayesian model averaging (bma) to estimate the weights

bmw <- bma(modlist, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 5, method = "cv")
plot(bmw)

# use the Model Confidence Set (mcs) to choose the models to combine

mcss <- mcs(modlist, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 5, B = 5000, l=3, alpha = 0.1,  method = "cv")


# combined forecasts
# stack

final <- CoMoMo(modlist, weight = stacked, data = DataStMoMo, ages.fit = agesFit,years.fit = yearsFit, h = 5)

# Bayesian

final <- CoMoMo(modlist, weight = bmw, data = DataStMoMo, ages.fit = agesFit,years.fit = yearsFit, h = 5)

# mcs

final <- CoMoMo(modlist, weight = mcss, data = DataStMoMo, ages.fit = agesFit,years.fit = yearsFit, h = 5)
```

## Questions 


Please feel free to open an issue with any questions you may have. You can contact me at s.kessy@unsw.edu.au
