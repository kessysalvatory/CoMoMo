#' @title Model combination using different methods.
#'
#' @export

fitCoMoMo <- function(models, data = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL, ages = NULL, years = NULL, h = NULL)

{

  # Check the forecast horizon

  if (h>0 && h!= as.integer(h)) stop("The forecast horizon h must be a positive integer.")

  year.eval <- max(years.fit)

  year.max <- year.eval + h

  # Check if more than one model is supplied

  if(length(models)<2) stop("Argument models needs to contain more than one model.")

  # Determine fitting ages and years are specified

  if(!is.null(data)) {
    if (class(data) != "StMoMoData") {
      stop("Argument data needs to be of class StMoMoData.")
    }

    ages <- data$ages
    years <- data$years

  } else {
    if (is.null(Dxt) || is.null(Ext)) {
      stop("Either argument data or arguments Dxt and Ext need to be provided.")
    }
    if (is.null(ages)) ages <- 1:nrow(Dxt)
    if (is.null(years)) years <- 1:ncol(Dxt)
  }

  if (is.null(ages.fit)) ages.fit <- ages
  if (is.null(years.fit)) years.fit <- years

# names of specified models

  Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

  if (is.null(h) && year.max <= max(data$years))

{

    fit <- lapply(year.eval:year.max, function(yearval) lapply(models, function(x)
      StMoMo::fit(x, data = data, ages.fit = ages.fit, years.fit = years.fit[1]:yearval)))

    return(fit)
  }

# produce predictions for multiple evaluation years.

 else if (!is.null(h) && year.max <= max(data$years))

{

  fitModels <-  lapply(year.eval:year.max, function(yearval) lapply(models, function(x) StMoMo::fit(x, data = data, ages.fit = ages.fit, years.fit = years.fit[1]:yearval)))

  forModels <-  lapply(1:h, function(s) lapply(1:length(models), function(x) forecast(fitModels[[s]][[x]], h = length(h:s))))

  forRates <-  lapply(1:h, function(i) lapply(forModels[[i]], function(x) {tidyr::pivot_longer(dplyr::bind_cols(ages = ages.fit, as.data.frame(x$rates)),
                                                                                               cols = 2:(length(h:i) + 1), names_to = "year", values_to = "rate")%>%dplyr::mutate(year = as.numeric(year), h = year - min(year) + 1)}))

  output <- dplyr::bind_rows(lapply(1:h, function(k) lapply(1:length(forRates[[k]]), function(x) dplyr::mutate(forRates[[k]][[x]], model = Specified_Models[[x]]))))

  return(output)

 }

# produce predictions for multiple evaluation years when year.max > max(data$years)

 else if (!is.null(h) && year.max > max(data$years))

 {

   fitModels <-  lapply(year.eval:max(data$years), function(yearval) lapply(models, function(x) StMoMo::fit(x, data = data, ages.fit = ages.fit, years.fit = years.fit[1]:yearval)))

   finalmodel <- c(fitModels, rep(tail(fitModels, n = 1),  year.max - max(data$years)))

   forModels <-  lapply(1:h, function(s) lapply(1:length(models), function(x) forecast(finalmodel[[s]][[x]], h = length(h:s))))

   forRates <-  lapply(1:h, function(i) lapply(forModels[[i]], function(x) {tidyr::pivot_longer(dplyr::bind_cols(ages = ages.fit, as.data.frame(x$rates)),
                                                                                               cols = 2:(length(h:i) + 1), names_to = "year", values_to = "rate")%>%dplyr::mutate(year = as.numeric(year), h = year - min(year) + 1)}))
   output <- dplyr::bind_rows(lapply(1:h, function(k) lapply(1:length(forRates[[k]]), function(x) dplyr::mutate(forRates[[k]][[x]], model = Specified_Models[[x]]))))

   return(output)

 }

}

#' @export

CoMoMo  <- function(weight, data = NULL,...) {
 UseMethod("CoMoMo")
}

#' @usage combines the mortality rate forecasts using simple model averaging

#' @return Returns an object of class \code{CoMoMo.default} with the following components:
#'
#' \item{Comb}{Returns the combinated forecasts for different horizon.}
#'
#' \item{Pred}{Returns the predictions of individual mortality models for different horizon.}
#'
#' \item{comb.method }{Returns the combination approach.}
#'
#' \item{mse}{Returns the mean squared errors.}
#'
#'  @examples comb <- CoMoMo(modlist, data = DataStMoMo, ages.fit = agesFit,years.fit = yearsFit, h = 5)
#'

#' @export

CoMoMo.default <- function(models, data = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL, ages = NULL, years = NULL, h = NULL) {

  # Check the forecast horizon

  if (h>0 && h!= as.integer(h)) stop("The forecast horizon h must be a positive integer.")

  # Determine fitting ages and years are specified

  if(!is.null(data)) {
    if (class(data) != "StMoMoData") {
      stop("Argument data needs to be of class StMoMoData.")
    }

    ages <- data$ages
    years <- data$years

  } else {
    if (is.null(Dxt) || is.null(Ext)) {
      stop("Either argument data or arguments Dxt and Ext need to be provided.")
    }
    if (is.null(ages)) ages <- 1:nrow(Dxt)
    if (is.null(years)) years <- 1:ncol(Dxt)
  }

  prediction <- fitCoMoMo(models = models,  data = data, ages.fit = ages.fit, years.fit = years.fit, h = h)

  simple <- prediction%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(average = exp(mean(log(rate))))%>%dplyr::ungroup()%>%
   tidyr::pivot_longer(cols = average, values_to = "rate", names_to = "model")

   allForecast_one <- dplyr::bind_rows(prediction,  simple)

  # observed rates

   obsRatesDF <- tidyr::pivot_longer(dplyr::bind_cols(ages = data$ages, as.data.frame(data$Dxt/data$Ext)),
                                   cols = 2:(length(data$years) + 1), names_to = "year", values_to = "obsrate")%>%dplyr::mutate(year = as.numeric(year))

   allForecast <- dplyr::left_join(allForecast_one, obsRatesDF)

  # compute errors

  mseDF <- allForecast%>%dplyr::mutate(err2 = (log(rate) - log(obsrate))^2)%>%dplyr::group_by(h, model)%>%dplyr::summarise(mse = mean(err2, na.rm = TRUE))%>%dplyr::ungroup()

  output <-list(comb =na.omit(simple), pred = prediction, comb.method = "simple", mse =  na.omit(mseDF))

  class(output) <- "CoMoMo.default"

  return(output)
}

#' @usage combines the mortality rate forecasts using bayesian model averaging

#' @return Returns an object of class \code{CoMoMo.bma} with the following components:
#'
#' \item{Comb}{Returns the combinated forecasts for different horizon.}
#'
#' \item{Pred}{Returns the predictions of individual mortality models for different horizon.}
#'
#' \item{comb.method }{Returns the combination approach.}
#'
#' \item{mse}{Returns the mean squared errors.}
#'
#'  @examples
#'
#'  weights <- bma(modlist, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 5, method = "cv")
#'  comb <- CoMoMo(modlist, weight =  weights, data = DataStMoMo, ages.fit = agesFit,years.fit = yearsFit, h = 5)
#'

#' @export

CoMoMo.bma <- function(models, data = NULL, weight = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL, ages = NULL, years = NULL, h = NULL) {


  # Check inputs

  if (class(weight) != "bma") {
    stop("Argument data needs to be of class bma.")
  }

  # Check the forecast horizon

  if (h>0 && h!= as.integer(h)) stop("The forecast horizon h must be a positive integer.")

   prediction <- fitCoMoMo(models = models, data = data, ages.fit = ages.fit, years.fit = years.fit, h = h)

   if (max(weight$weights$h)==h)

   {

    pbma <- prediction%>%dplyr::left_join(weight$weights)%>%dplyr::mutate(forbma = log(rate)*weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(bayesian = exp(sum(forbma)))%>%dplyr::ungroup()%>%
    tidyr::pivot_longer(cols = bayesian, values_to = "rate", names_to = "model")

    allForecast_one <- dplyr::bind_rows(prediction,  pbma)

    # observed rates

    obsRatesDF <- tidyr::pivot_longer(dplyr::bind_cols(ages = data$ages, as.data.frame(data$Dxt/data$Ext)),cols = 2:(length(data$years)+1),
                                      names_to = "year", values_to = "obsrate")%>%dplyr::mutate(year = as.numeric(year))

    allForecast <- dplyr::left_join(allForecast_one, obsRatesDF)

    # compute errors

    mseDF <- allForecast%>%dplyr::mutate(err2 = (log(rate) - log(obsrate))^2)%>%dplyr::group_by(h, model)%>%dplyr::summarise(mse = mean(err2, na.rm = TRUE))%>%dplyr::ungroup()

    output <-list(comb = na.omit(pbma), pred = prediction, comb_method = "bma", weights = weight$weights, mse = na.omit(mseDF))

    class(output) <- "CoMoMo.bma"

    return(output)

   }

 else if (max(weight$weights$h) < h)

 {

   # updates the weights here

   start <- max(weight$weights$h)

   interval<- h - start

   weightsDF <- (dplyr::bind_rows(lapply(rep(list(as.data.frame(tail(weight$weights, length(models)))), interval), function(x) x%>%dplyr::mutate(model = modelNames)))%>%dplyr::mutate(h=rep(1:interval + start, each=length(models))))[,c("h","weights","model")]

   weightsDFall <- rbind(weight$weights, weightsDF)

   weight$weights <- weightsDFall

   pbma <- prediction%>%dplyr::left_join(weight$weights)%>%dplyr::mutate(forbma = log(rate)*weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(bayesian = exp(sum(forbma)))%>%dplyr::ungroup()%>%
     tidyr::pivot_longer(cols = Bayesian, values_to = "rate", names_to = "model")

   allForecast_one <- dplyr::bind_rows(prediction,  pbma)

   # observed rates

   obsRatesDF <- tidyr::pivot_longer(dplyr::bind_cols(ages = data$ages, as.data.frame(data$Dxt/data$Ext)),cols = 2:(length(data$years)+1),
                                     names_to = "year", values_to = "obsrate")%>%dplyr::mutate(year = as.numeric(year))

   allForecast <- dplyr::left_join(allForecast_one, obsRatesDF)

   # compute errors

   mseDF <- allForecast%>%dplyr::mutate(err2 = (log(rate) - log(obsrate))^2)%>%dplyr::group_by(h, model)%>%dplyr::summarise(mse = mean(err2, na.rm = TRUE))%>%dplyr::ungroup()

   output <-list(comb = na.omit(pbma), pred = prediction, comb_method = "bma", mse = na.omit(mseDF))

   class(output) <- "CoMoMo.bma"

   return(output)

 }

   else if (max(weight$weights$h) > h)

   {

     # updates the weights here

     weight$weights <- head(weight$weights, length(models)*h)

     pbma <- prediction%>%dplyr::left_join(weight$weights)%>%dplyr::mutate(forbma = log(rate)*weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(bayesian = exp(sum(forbma)))%>%dplyr::ungroup()%>%
       tidyr::pivot_longer(cols = Bayesian, values_to = "rate", names_to = "model")

     allForecast_one <- dplyr::bind_rows(prediction,  pbma)

     # observed rates

     obsRatesDF <- tidyr::pivot_longer(dplyr::bind_cols(ages = data$ages, as.data.frame(data$Dxt/data$Ext)),cols = 2:(length(data$years)+1),
                                       names_to = "year", values_to = "obsrate")%>%dplyr::mutate(year = as.numeric(year))

     allForecast <- dplyr::left_join(allForecast_one, obsRatesDF)

     # compute errors

     mseDF <- allForecast%>%dplyr::mutate(err2 = (log(rate) - log(obsrate))^2)%>%dplyr::group_by(h, model)%>%dplyr::summarise(mse = mean(err2, na.rm = TRUE))%>%dplyr::ungroup()

     output <-list(comb = na.omit(pbma), pred = prediction, comb_method = "bma", mse = na.omit(mseDF))

     class(output) <- "CoMoMo.bma"

     return(output)

   }
}

#' @export
#'
#' @usage combines the mortality rate forecasts using stacked regression ensemble

#' @return Returns an object of class \code{CoMoMo.stack} with the following components:
#'
#' \item{Comb}{Returns the combinated forecasts for different horizon.}
#'
#' \item{Pred}{Returns the predictions of individual mortality models for different horizon.}
#'
#' \item{comb.method }{Returns the combination approach.}
#'
#' \item{mse}{Returns the mean squared errors.}
#'
#' \item{metalearner}{Returns the used metalearner.}
#'
#'  @examples
#'
#'  weights <- stack(modlist, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 5, metalearner = "nnls")
#'  comb <- CoMoMo(modlist, weight =  weights, data = DataStMoMo, ages.fit = agesFit,years.fit = yearsFit, h = 5)
#'


CoMoMo.stack <- function(models, data = NULL, weight = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL, ages = NULL, years = NULL, h = NULL) {

  # Check inputs

  if (class(weight) != "stack") {
    stop("Argument data needs to be of class stack.")
  }

  # Check the forecast horizon

  if (h>0 && h!= as.integer(h)) stop("The forecast horizon h must be a positive integer.")

  prediction <- fitCoMoMo(models = models,  data = data, ages.fit = ages.fit, years.fit = years.fit, h = h)

 if (max(weight$weights$h)==h)

 {

  pstack <- prediction%>%dplyr::filter (h <= max(weight$weights$h))%>%dplyr::left_join(weight$weights)%>%dplyr::mutate(forstack= log(rate)*weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(stack = exp(sum(forstack)))%>%dplyr::ungroup()%>%
    tidyr::pivot_longer(cols = stack, values_to = "rate", names_to = "model")

   allForecast_one <- dplyr::bind_rows(prediction, pstack)

   # observed rates

   obsRatesDF <- tidyr::pivot_longer(dplyr::bind_cols(ages = data$ages, as.data.frame(data$Dxt/data$Ext)),
                             cols = 2:(length(data$years) + 1), names_to = "year", values_to = "obsrate")%>%dplyr::mutate(year = as.numeric(year))

   allForecast <- dplyr::left_join(allForecast_one, obsRatesDF)

   # compute errors

   mseDF <- allForecast%>%dplyr::mutate(err2 = (log(rate) - log(obsrate))^2)%>%dplyr::group_by(h, model)%>%dplyr::summarise(mse = mean(err2, na.rm = TRUE))%>%dplyr::ungroup()

   output <-list(comb = na.omit(pstack), pred = prediction, comb_method = "stack", metalearner = weight$metalearner, mse = na.omit(mseDF))

   class(output) <- "CoMoMo.stack"

   return(output)

 }


else if (max(weight$weights$h) < h)

{

  # updates the weights here

  start <- max(weight$weights$h)

  interval<- h-start

  weightsDF <- (dplyr::bind_rows(lapply(rep(list(as.data.frame(tail(weight$weights, length(models)))), interval), function(x) x%>%dplyr::mutate(model = modelNames)))%>%dplyr::mutate(h=rep(1:interval + start, each=length(models))))[,c("h","weights","model")]

  weightsDFall <- rbind(weight$weights, weightsDF)

  weight$weights <- weightsDFall

  pstack <- prediction%>%dplyr::left_join(weight$weights)%>%dplyr::mutate(forstack= log(rate)*weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(stack = exp(sum(forstack)))%>%dplyr::ungroup()%>%
    tidyr::pivot_longer(cols = stack, values_to = "rate", names_to = "model")

  allForecast_one <- dplyr::bind_rows(prediction, pstack)

  # observed rates

  obsRatesDF <- tidyr::pivot_longer(dplyr::bind_cols(ages = data$ages, as.data.frame(data$Dxt/data$Ext)),
                                    cols = 2:(length(data$years) + 1), names_to = "year", values_to = "obsrate")%>%dplyr::mutate(year = as.numeric(year))

  allForecast <- dplyr::left_join(allForecast_one, obsRatesDF)

  # compute errors

  mseDF <- allForecast%>%dplyr::mutate(err2 = (log(rate) - log(obsrate))^2)%>%dplyr::group_by(h, model)%>%dplyr::summarise(mse = mean(err2, na.rm = TRUE))%>%dplyr::ungroup()

  output <-list(comb = na.omit(pstack), pred = prediction, comb_method = "stack", metalearner = weight$metalearner, mse = na.omit(mseDF))

   class(output) <- "CoMoMo.stack"

   return(output)
}

  else if (max(weight$weights$h) > h)

  {

    # updates the weights here

    weight$weights <- head(weight$weights, length(models)*h)

    pstack <- prediction%>%dplyr::left_join(weight$weights)%>%dplyr::mutate(forstack= log(rate)*weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(stack = exp(sum(forstack)))%>%dplyr::ungroup()%>%
      tidyr::pivot_longer(cols = stack, values_to = "rate", names_to = "model")

    allForecast_one <- dplyr::bind_rows(prediction, pstack)

    # observed rates

    obsRatesDF <- tidyr::pivot_longer(dplyr::bind_cols(ages = data$ages, as.data.frame(data$Dxt/data$Ext)),
                                      cols = 2:(length(data$years) + 1), names_to = "year", values_to = "obsrate")%>%dplyr::mutate(year = as.numeric(year))

    allForecast <- dplyr::left_join(allForecast_one, obsRatesDF)

    # compute errors

    mseDF <- allForecast%>%dplyr::mutate(err2 = (log(rate) - log(obsrate))^2)%>%dplyr::group_by(h, model)%>%dplyr::summarise(mse = mean(err2, na.rm = TRUE))%>%dplyr::ungroup()

    output <-list(comb = na.omit(pstack), pred = prediction, comb_method = "stack", metalearner = weight$metalearner, mse = na.omit(mseDF))

    class(output) <- "CoMoMo.stack"

    return(output)
  }

}

#' @export
#'
#' @usage combines the mortality rate forecasts using model confidence set

#' @return Returns an object of class \code{CoMoMo.mcs} with the following components:
#'
#' \item{Comb}{Returns the combinated forecasts for different horizon.}
#'
#' \item{Pred}{Returns the predictions of individual mortality models for different horizon.}
#'
#' \item{comb.method }{Returns the combination approach.}
#'
#' \item{mse}{Returns the mean squared errors.}
#'
#'  @examples
#'
#'  weights <- weights <- mcs(modlist, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 5, B = 5000, l=3, alpha = 0.1,  method = "cv")
#'  comb <- CoMoMo(modlist, weight =  weights, data = DataStMoMo, ages.fit = agesFit,years.fit = yearsFit, h = 5)
#'

CoMoMo.mcs <- function(models, data = NULL, weight = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL, ages = NULL, years = NULL, h = NULL) {


  # Check inputs

  if (class(weight) != "mcs") {
    stop("Argument data needs to be of class mcs.")
  }

  # Check the forecast horizon

  if (h>0 && h!= as.integer(h)) stop("The forecast horizon h must be a positive integer.")

  prediction <- fitCoMoMo(models = models,  data = data, ages.fit = ages.fit, years.fit = years.fit, h = h)


  if (max(weight$weights$h)==h)

  {

  pmcs <- prediction%>%dplyr::left_join(weight$weights)%>%dplyr::mutate(formcs = log(rate)*weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(mcs = exp(sum(formcs)))%>%dplyr::ungroup()%>%
    tidyr::pivot_longer(cols = mcs, values_to = "rate", names_to = "model")

  allForecast_one <- dplyr::bind_rows(prediction, pmcs)

  # observed rates

  obsRatesDF <- tidyr::pivot_longer(dplyr::bind_cols(ages = data$ages, as.data.frame(data$Dxt/data$Ext)),
                                    cols = 2:(length(data$years)+1), names_to = "year", values_to = "obsrate")%>%dplyr::mutate(year = as.numeric(year))

  allForecast <- dplyr::left_join(allForecast_one, obsRatesDF)

  # compute errors

  mseDF <- allForecast%>%dplyr::mutate(err2 = (log(rate) - log(obsrate))^2)%>%dplyr::group_by(h, model)%>%dplyr::summarise(mse = mean(err2, na.rm = TRUE))%>%dplyr::ungroup()

  output <-list(comb = na.omit(pmcs), pred = prediction, comb_method = "mcs",  mse = na.omit(mseDF))

  class(output) <- "CoMoMo.mcs"

  return(output)

}


  else if (max(weight$weights$h) < h)

  {

    # updates the weights here

    start <- max(weight$weights$h)

    interval<- h-start

    weightsDF <- (dplyr::bind_rows(lapply(rep(list(as.data.frame(tail(weight$weights, length(models)))), interval), function(x) x%>%dplyr::mutate(model = modelNames)))%>%dplyr::mutate(h=rep(1:interval + start, each=length(models))))[,c("h","weights","model")]

    weightsDFall <- rbind(weight$weights, weightsDF)

    weight$weights <- weightsDFall

    pmcs <- prediction%>%dplyr::left_join(weight$weights)%>%dplyr::mutate(formcs = log(rate)*weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(mcs = exp(sum(formcs)))%>%dplyr::ungroup()%>%
      tidyr::pivot_longer(cols = mcs, values_to = "rate", names_to = "model")

    allForecast_one <- dplyr::bind_rows(prediction, pmcs)

    # observed rates

    obsRatesDF <- tidyr::pivot_longer(dplyr::bind_cols(ages = data$ages, as.data.frame(data$Dxt/data$Ext)),
                                      cols = 2:(length(data$years)+1), names_to = "year", values_to = "obsrate")%>%dplyr::mutate(year = as.numeric(year))

    allForecast <- dplyr::left_join(allForecast_one, obsRatesDF)

    # compute errors

    mseDF <- allForecast%>%dplyr::mutate(err2 = (log(rate) - log(obsrate))^2)%>%dplyr::group_by(h, model)%>%dplyr::summarise(mse = mean(err2, na.rm = TRUE))%>%dplyr::ungroup()

    output <-list(comb = na.omit(pmcs), pred = prediction, comb_method = "mcs",  mse = na.omit(mseDF))

    class(output) <- "CoMoMo.mcs"

    return(output)

  }

  else if (max(weight$weights$h) > h)

  {

    # updates the weights here

    weight$weights <- head(weight$weights, length(models)*h)

    pmcs <- prediction%>%dplyr::left_join(weight$weights)%>%dplyr::mutate(formcs = log(rate)*weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(mcs = exp(sum(formcs)))%>%dplyr::ungroup()%>%
      tidyr::pivot_longer(cols = mcs, values_to = "rate", names_to = "model")

    allForecast_one <- dplyr::bind_rows(prediction, pmcs)

    # observed rates

    obsRatesDF <- tidyr::pivot_longer(dplyr::bind_cols(ages = data$ages, as.data.frame(data$Dxt/data$Ext)),
                                      cols = 2:(length(data$years)+1), names_to = "year", values_to = "obsrate")%>%dplyr::mutate(year = as.numeric(year))

    allForecast <- dplyr::left_join(allForecast_one, obsRatesDF)

    # compute errors

    mseDF <- allForecast%>%dplyr::mutate(err2 = (log(rate) - log(obsrate))^2)%>%dplyr::group_by(h, model)%>%dplyr::summarise(mse = mean(err2, na.rm = TRUE))%>%dplyr::ungroup()

    output <-list(comb = na.omit(pmcs), pred = prediction, comb_method = "mcs",  mse = na.omit(mseDF))

    class(output) <- "CoMoMo.mcs"

    return(output)

  }

}


















