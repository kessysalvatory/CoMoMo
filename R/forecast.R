#' @export
#'
#' @importFrom forecast forecast
#'

forecast.CoMoMo <- function(object, h = h, ...)

{
  Specified_Models <- lapply(1:length(object$model.fits), function(x) names(object$model.fits[x]))

  forModels <-  lapply(1:length(object$model.fits), function(x) forecast(object$model.fits[[x]], h = h))

  if (h==1)

  {
    forData <- lapply(forModels, function(x) dplyr::bind_cols(ages = object$model.fits[[1]]$ages, as.data.frame(x$rates)))

    for (i in 1:length(object$model.fits))

    {
      names(forData[[i]])[2] <- as.character(forModels[[1]]$years)
    }

    forRates <- lapply(forData, function(x) {tidyr::pivot_longer(x, cols = 2:(h + 1), names_to = "year", values_to = "rate")%>%dplyr::mutate(year = as.numeric(year), h = year - min(year) + 1)})

  }

  else {

    forRates <- lapply(forModels, function(x) {tidyr::pivot_longer(dplyr::bind_cols(ages = object$model.fits[[1]]$ages, as.data.frame(x$rates)), cols = 2:(h + 1), names_to = "year",
                                                                   values_to = "rate")%>%dplyr::mutate(year = as.numeric(year), h = year - min(year) + 1)})
  }

  prediction <- dplyr::bind_rows(lapply(1:length(forRates), function(x) dplyr::mutate(forRates[[x]], model = Specified_Models[[x]])))

  if (max(object$weight$weights$h)==h)

  {
    mortfor <- prediction%>%dplyr::left_join(object$weight$weights)%>%dplyr::mutate(forrates = log(rate)*model.weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(comb = exp(sum(forrates)))%>%dplyr::ungroup()%>%
      tidyr::pivot_longer(cols = comb, values_to = "rate", names_to = "model")

    allForecast_one <- dplyr::bind_rows(prediction, mortfor)

    output <-list(comb.rates = na.omit(mortfor), model.rates = prediction, comb.method = object$weight$comb.method)

    return(output)

  }

  else if (max(object$weight$weights$h) < h)

  {
    # updates the weights here when having weights for fewer horizons

    start.horizon <- max(object$weight$weights$h)

    interval<- h - start.horizon

    weightsDF <- (dplyr::bind_rows(lapply(rep(list(as.data.frame(tail(object$weight$weights, length(models)))), interval), function(x) x%>%dplyr::mutate(model = names(models))))%>%dplyr::mutate(h = rep(1:interval + start.horizon, each=length(models))))[, c("h","model.weights","model")]

    weightsDFall <- rbind(object$weight$weights, weightsDF)

    object$weight$weights <- weightsDFall

    mortfor <- prediction%>%dplyr::left_join(object$weight$weights)%>%dplyr::mutate(forrates = log(rate)*model.weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(comb = exp(sum(forrates)))%>%dplyr::ungroup()%>%
      tidyr::pivot_longer(cols = comb, values_to = "rate", names_to = "model")

    output <- structure(list(comb.rates = na.omit(mortfor), model.rates = prediction, comb.method = object$weight$comb.method))

    return(output)

  }

  else if (max(object$weight$weights$h) > h)

  {
    # updates the weights when predicting shorter horizon than the horizon

    # for the weights

    object$weight$weights <- head(object$weight$weights, length(models)*h)

    mortfor <- prediction%>%dplyr::left_join(object$weight$weights)%>%dplyr::mutate(forrates = log(rate)*model.weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(comb = exp(sum(forrates)))%>%dplyr::ungroup()%>%
      tidyr::pivot_longer(cols = comb, values_to = "rate", names_to = "model")

    output <-list(comb.rates = na.omit(mortfor), model.rates = prediction, comb.method = object$weight$comb.method)

    return(output)

  }
}

#' @export
#'
forecast.CoMoMo.simple <- function(object, h = h, ...)

{
  # names of specified models

  Specified_Models <- lapply(1:length(object$model.fits), function(x) names(object$model.fits[x]))

  forModels <-  lapply(1:length(object$model.fits), function(x) forecast(object$model.fits[[x]], h = h))

  if (h==1)

  {
    forData <- lapply(forModels, function(x) dplyr::bind_cols(ages = object$model.fits[[1]]$ages, as.data.frame(x$rates)))

    for (i in 1:length(object$model.fits))

    {
      names(forData[[i]])[2] <- as.character(forModels[[1]]$years)
    }

    forRates <- lapply(forData, function(x) {tidyr::pivot_longer(x, cols = 2:(h + 1), names_to = "year", values_to = "rate")%>%dplyr::mutate(year = as.numeric(year), h = year - min(year) + 1)})

  }

  else {

    forRates <- lapply(forModels, function(x) {tidyr::pivot_longer(dplyr::bind_cols(ages = object$model.fits[[1]]$ages, as.data.frame(x$rates)), cols = 2:(h + 1), names_to = "year",
                                                                   values_to = "rate")%>%dplyr::mutate(year = as.numeric(year), h = year - min(year) + 1)})
  }

  prediction <- dplyr::bind_rows(lapply(1:length(forRates), function(x) dplyr::mutate(forRates[[x]], model = Specified_Models[[x]])))

  simple <- prediction%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(average = exp(mean(log(rate))))%>%dplyr::ungroup()%>%
    tidyr::pivot_longer(cols = average, values_to = "rate", names_to = "model")

  output <-list(comb.rates = na.omit(simple), model.rates = prediction, comb.method = "simple")

  return(output)

}




