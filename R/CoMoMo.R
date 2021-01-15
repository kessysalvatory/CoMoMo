#' @title Model combination using different methods.
#'
#' @export

fitCoMoMo <- function(models, data = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL, ages = NULL, years = NULL, h = NULL)

{
  # Check if more than one model is supplied

  if(length(models)<2) stop("Argument models needs to contain more than one model.")

  # check if the supplied models are different

  if(length(unique(unname(models))) < length(models)) stop("Models must be different.")

  # Check the forecast horizon

  if (h<0 || h!= as.integer(h)) stop("The forecast horizon h must be a positive integer.")

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

year.eval <- max(years.fit)

year.max <- year.eval + h

# produce predictions for multiple evaluation years.

 if (!is.null(h) && year.max <= max(data$years))

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

#' @usage combines the mortality rate forecasts using different model combinations methods

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
#'  @examples
#'
#'  comb <- CoMoMo(modlist, data = DataStMoMo, ages.fit = agesFit,years.fit = yearsFit, h = 5)
#'
#'  weights <- bma(modlist, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 5, method = "cv")
#'
#'  comb <- CoMoMo(modlist, weight =  weights, data = DataStMoMo, ages.fit = agesFit,years.fit = yearsFit, h = 5)
#'
#'
#' @export

CoMoMo  <- function(weight, data = NULL,...) {
  UseMethod("CoMoMo")
}


#' @export

CoMoMo.default <- function(models, data = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL, ages = NULL, years = NULL, h = NULL) {

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

  output <-list(comb.rates = na.omit(simple), model.rates = prediction, comb.method = "simple", mse = na.omit(mseDF))

  class(output) <- "CoMoMo.default"

  return(output)
}

#' @export

CoMoMo.weight <- function(models, data = NULL, weight = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL, ages = NULL, years = NULL, h = NULL) {

  prediction <- fitCoMoMo(models = models, data = data, ages.fit = ages.fit, years.fit = years.fit, h = h)

   if (max(weight$weights$h)==h)

   {
     mortfor <- prediction%>%dplyr::left_join(weight$weights)%>%dplyr::mutate(forrates = log(rate)*model.weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(comb = exp(sum(forrates)))%>%dplyr::ungroup()%>%
       tidyr::pivot_longer(cols = comb, values_to = "rate", names_to = "model")

     allForecast_one <- dplyr::bind_rows(prediction,  mortfor)

     # observed rates

     obsRatesDF <- tidyr::pivot_longer(dplyr::bind_cols(ages = data$ages, as.data.frame(data$Dxt/data$Ext)),cols = 2:(length(data$years)+1),
                                       names_to = "year", values_to = "obsrate")%>%dplyr::mutate(year = as.numeric(year))

     allForecast <- dplyr::left_join(allForecast_one, obsRatesDF)

     # compute errors

     mseDF <- allForecast%>%dplyr::mutate(err2 = (log(rate) - log(obsrate))^2)%>%dplyr::group_by(h, model)%>%dplyr::summarise(mse = mean(err2, na.rm = TRUE))%>%dplyr::ungroup()

     output <-list(comb.rates = na.omit(mortfor), model.rates = prediction, comb.method = weight$comb.method, mse = na.omit(mseDF))


     class(output) <- "CoMoMo.weight"

     return(output)

   }

 else if (max(weight$weights$h) < h)

 {
   # updates the weights here when having weights for fewer horizons

   start.horizon <- max(weight$weights$h)

   interval<- h - start.horizon

   weightsDF <- (dplyr::bind_rows(lapply(rep(list(as.data.frame(tail(weight$weights, length(models)))), interval), function(x) x%>%dplyr::mutate(model = names(models))))%>%dplyr::mutate(h = rep(1:interval + start.horizon, each=length(models))))[, c("h","model.weights","model")]

   weightsDFall <- rbind(weight$weights, weightsDF)

   weight$weights <- weightsDFall

   mortfor <- prediction%>%dplyr::left_join(weight$weights)%>%dplyr::mutate(forrates = log(rate)*model.weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(comb = exp(sum(forrates)))%>%dplyr::ungroup()%>%
     tidyr::pivot_longer(cols = comb, values_to = "rate", names_to = "model")

   allForecast_one <- dplyr::bind_rows(prediction, mortfor)

   # observed rates

   obsRatesDF <- tidyr::pivot_longer(dplyr::bind_cols(ages = data$ages, as.data.frame(data$Dxt/data$Ext)),cols = 2:(length(data$years)+1),
                                     names_to = "year", values_to = "obsrate")%>%dplyr::mutate(year = as.numeric(year))

   allForecast <- dplyr::left_join(allForecast_one, obsRatesDF)

   # compute errors

   mseDF <- allForecast%>%dplyr::mutate(err2 = (log(rate) - log(obsrate))^2)%>%dplyr::group_by(h, model)%>%dplyr::summarise(mse = mean(err2, na.rm = TRUE))%>%dplyr::ungroup()

   output <-list(comb.rates = na.omit(mortfor), model.rates = prediction, comb.method = weight$comb.method, mse = na.omit(mseDF))

   class(output) <- "CoMoMo.weight"

   return(output)

 }

   else if (max(weight$weights$h) > h)

   {
     # updates the weights when predicting shorter horizon than the horizon

     # for the weights

     weight$weights <- head(weight$weights, length(models)*h)

     mortfor <- prediction%>%dplyr::left_join(weight$weights)%>%dplyr::mutate(forrates = log(rate)*model.weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(comb = exp(sum(forrates)))%>%dplyr::ungroup()%>%
       tidyr::pivot_longer(cols = comb, values_to = "rate", names_to = "model")

     allForecast_one <- dplyr::bind_rows(prediction,  mortfor)

     # observed rates

     obsRatesDF <- tidyr::pivot_longer(dplyr::bind_cols(ages = data$ages, as.data.frame(data$Dxt/data$Ext)),cols = 2:(length(data$years)+1),
                                       names_to = "year", values_to = "obsrate")%>%dplyr::mutate(year = as.numeric(year))

     allForecast <- dplyr::left_join(allForecast_one, obsRatesDF)

     # compute errors

     mseDF <- allForecast%>%dplyr::mutate(err2 = (log(rate) - log(obsrate))^2)%>%dplyr::group_by(h, model)%>%dplyr::summarise(mse = mean(err2, na.rm = TRUE))%>%dplyr::ungroup()

     output <-list(comb.rates = na.omit(mortfor), model.rates = prediction, comb.method = weight$comb.method, mse = na.omit(mseDF))

     class(output) <- "CoMoMo.weight"

     return(output)

   }
}

















