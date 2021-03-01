
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
#'  comb <- CoMoMo(modlist, data = DataStMoMo, ages.fit = ages.fit, years.fit = years.fit, h = 5)
#'
#'  weights <- bma(modlist, data = DataStMoMo, ages.fit = ages.fit, years.fit = years.fit, h = 5, method = "cv")
#'
#'  comb <- CoMoMo(modlist, weight =  weights, data = DataStMoMo,  ages.fit = ages.fit, years.fit = years.fit, h = 5)
#'
#'
#' @export

CoMoMo  <- function(object, weight)

{
  output <- structure(list(model.fits = object, weight = weight))

  class(output) <- "CoMoMo"

  return(output)

}

# forecast.CoMoMo.default <- function(object, h = h) {
#
#   # names of specified models
#
#   Specified_Models <- lapply(1:length(object$model.fits), function(x) names(object$model.fits[x]))
#
#   forModels <-  lapply(1:length(object$model.fits), function(x) forecast(object$model.fits[[x]], h = h))
#
#    if (h ==1)
#
#   {
#
#    tempdata <- lapply(forModels, function(x) dplyr::bind_cols(ages = object[[1]]$ages, as.data.frame(x$rates)))
#
#     for (i in 1:length(object$model.fits))
#
#    {
#       names(tempdata[[i]])[2] <- as.character(forModels[[1]]$years)
#
#     }
#
#    forRates <- lapply(tempdata, function(x) {tidyr::pivot_longer(x, cols = 2:(h + 1), names_to = "year", values_to = "rate") %>%dplyr::mutate(year = as.numeric(year), h = year - min(year) + 1)})
#
#     }
#
#   else {
#
#  forRates <- lapply(forModels, function(x) {tidyr::pivot_longer(dplyr::bind_cols(ages = object[[1]]$ages, as.data.frame(x$rates)), cols = 2:(h + 1), names_to = "year",
#                                                                    values_to = "rate")%>%dplyr::mutate(year = as.numeric(year), h = year - min(year) + 1)})
#
#  }
#
#    prediction <- dplyr::bind_rows(lapply(1:length(forRates), function(x) dplyr::mutate(forRates[[x]], model = Specified_Models[[x]])))
#
#   simple <- prediction%>%dplyr::group_by(ages, year, h)%>%dplyr::summarise(average = exp(mean(log(rate))))%>%dplyr::ungroup()%>%
#    tidyr::pivot_longer(cols = average, values_to = "rate", names_to = "model")
#
#    output <-list(comb.rates = na.omit(simple), model.rates = prediction, comb.method = "simple")
#
#    return(output)
#
#  }

#'

















