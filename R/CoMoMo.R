#' @title Model combination using different methods.
#'
#' @export

fitCoMoMo <- function(models, data = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL, ages = NULL, years = NULL)

{
  # Check if more than one model is supplied

  if(length(models)<2) stop("Argument models needs to contain more than one model.")

  # check if the supplied models are different

  if(length(unique(unname(models))) < length(models)) stop("Models must be different.")

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

   fitModels <-  lapply(models, function(x) StMoMo::fit(x,  Dxt = Dxt, Ext = Ext, ages = ages, years = years,
                      data = data, ages.fit = ages.fit, years.fit = years.fit))

   return(fitModels)

}

#' @usage combines the mortality rate forecasts using different model combinations methods
#'
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
#'

CoMoMo  <- function(weight,...) {
  UseMethod("CoMoMo")
}

#' @export

CoMoMo.default <- function(object, h = h) {

  # names of specified models

  Specified_Models <- lapply(1:length(object), function(x) names(object[x]))

  forModels <-  lapply(1:length(object), function(x) forecast(object[[x]], h = h))

   if (h ==1)

  {

    tempdata <- lapply(forModels, function(x) dplyr::bind_cols(ages = object[[1]]$ages, as.data.frame(x$rates)))

    for (i in 1:length(object))

    {
      names(tempdata[[i]])[2] <- as.character(forModels[[1]]$years)
    }

    forRates <- lapply(tempdata, function(x) {tidyr::pivot_longer(x, cols = 2:(h + 1), names_to = "year", values_to = "rate") %>%dplyr::mutate(year = as.numeric(year), h = year - min(year) + 1)})

   }

  else {

  forRates <- lapply(forModels, function(x) {tidyr::pivot_longer(dplyr::bind_cols(ages = object[[1]]$ages, as.data.frame(x$rates)), cols = 2:(h + 1), names_to = "year",
                                                                   values_to = "rate")%>%dplyr::mutate(year = as.numeric(year), h = year - min(year) + 1)})

  }

  prediction <- dplyr::bind_rows(lapply(1:length(forRates), function(x) dplyr::mutate(forRates[[x]], model = Specified_Models[[x]])))

  simple <- prediction%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(average = exp(mean(log(rate))))%>%dplyr::ungroup()%>%
   tidyr::pivot_longer(cols = average, values_to = "rate", names_to = "model")

   output <-list(comb.rates = na.omit(simple), model.rates = prediction, comb.method = "simple")

  return(output)
}

#' @export

CoMoMo.weight <- function(object, h = h, weight = weight) {

  # names of specified models

  Specified_Models <- lapply(1:length(object), function(x) names(object[x]))

  forModels <-  lapply(1:length(object), function(x) forecast(object[[x]], h = h))

  if (h ==1)

  {

    tempdata <- lapply(forModels, function(x) dplyr::bind_cols(ages = object[[1]]$ages, as.data.frame(x$rates)))

    for (i in 1:length(object))

    {
      names(tempdata[[i]])[2] <- as.character(forModels[[1]]$years)
    }

    forRates <- lapply(tempdata, function(x) {tidyr::pivot_longer(x, cols = 2:(h + 1), names_to = "year", values_to = "rate") %>%dplyr::mutate(year = as.numeric(year), h = year - min(year) + 1)})

  }

  else {

  forRates <- lapply(forModels, function(x) {tidyr::pivot_longer(dplyr::bind_cols(ages = object[[1]]$ages, as.data.frame(x$rates)), cols = 2:(h + 1), names_to = "year",
                                                                 values_to = "rate")%>%dplyr::mutate(year = as.numeric(year), h = year - min(year) + 1)})
  }

  prediction <- dplyr::bind_rows(lapply(1:length(forRates), function(x) dplyr::mutate(forRates[[x]], model = Specified_Models[[x]])))

  if (max(weight$weights$h)==h)

  {
    mortfor <- prediction%>%dplyr::left_join(weight$weights)%>%dplyr::mutate(forrates = log(rate)*model.weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(comb = exp(sum(forrates)))%>%dplyr::ungroup()%>%
      tidyr::pivot_longer(cols = comb, values_to = "rate", names_to = "model")

    allForecast_one <- dplyr::bind_rows(prediction,  mortfor)

    output <-list(comb.rates = na.omit(mortfor), model.rates = prediction, comb.method = weight$comb.method)

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

    output <-list(comb.rates = na.omit(mortfor), model.rates = prediction, comb.method = weight$comb.method)

    return(output)

  }

  else if (max(weight$weights$h) > h)

  {
    # updates the weights when predicting shorter horizon than the horizon

    # for the weights

    weight$weights <- head(weight$weights, length(models)*h)

    mortfor <- prediction%>%dplyr::left_join(weight$weights)%>%dplyr::mutate(forrates = log(rate)*model.weights)%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(comb = exp(sum(forrates)))%>%dplyr::ungroup()%>%
      tidyr::pivot_longer(cols = comb, values_to = "rate", names_to = "model")

    output <-list(comb.rates = na.omit(mortfor), model.rates = prediction, comb.method = weight$comb.method)

    return(output)

  }
}

















