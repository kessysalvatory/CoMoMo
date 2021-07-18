#' @title Generating the metadata. 
#'
#' @description Stacked regression ensemble often proceeds in two steps in generating the final predictions. 
#' The first step consists of multiple base models which separately generate cross-validated predictions from the training data. 
#' The predictions from various models and the observed response variable constitute the metadata.
 #' @references
#'
#' Kessy, Salvatory, Michael Sherris, Andrés Villegas, and Jonathan Ziveyi. 2021. 
#' “Mortality Forecasting Using Stacked Regression Ensembles.” SSRN Electronic Journal. https://doi.org/10.2139/ssrn.3823511.
#' 
#'  @param h number of years for forecasting horizon. 
#'
#'  @return Returns an object of the class \code{stackmeta} with the following components: 
#'  \item{metadata}{A list of the metadata for different forecasting horizons.}
#'  \item{cvmse}{Returns the cross-validation mean squared error for each all the indiviual models for different horizons.}
#'  \item{models}{A list of the models used to generate the metadata.}

#' @export
#'
metadata <- function(models, data = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL, ages = NULL, years = NULL, h = NULL)

{ # the function compute the metadata using block-cross-validation
  
  # check if the supplied models are different

  if(length(unique(unname(models))) < length(models)) stop("Models must be different.")

  # # Check the forecast horizon

  if (h<0 || h!= as.integer(h)) stop("The forecast horizon h must be a positive integer.")

  # checking the horizons that weights can be generated

  if (h > (max(years.fit) - years.fit[1]-1)) {
    stop("The time series must be longer than h.")
  }

  # Check if more than one model is supplied

  if(length(models)<2) stop("Argument models needs to contain more than one model.")

  # Assign models names

  if (is.null(names(models))) stop("Assign names to the models")

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
    if (is.null(ages))  ages <- 1:nrow(Dxt)
    if (is.null(years)) years <- 1:ncol(Dxt)
  }

  if (is.null(ages.fit)) ages.fit <- ages
  if (is.null(years.fit)) years.fit <- years

  # names of specified models

  Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

  output0 <- cvloss(models = models, data = data, ages.fit = ages.fit, years.fit = years.fit, h = h, Dxt = Dxt, Ext = Ext, ages = ages, years = years)

  cvModels1 <- output0$cvModels; cvRates1 <- output0$cvRates

  xtrain0 <- output0$xtrain;  naRows <- output0$naRows;  CVerror <-  output0$CVE

  xtrain <- lapply(1:h, function(x)  xtrain0[[x]][-naRows[[x]],])

  ytrain0 <- lapply(1:h, function(x) lapply(1:length(models), function(m) log(as.matrix(tidyr::pivot_longer(dplyr::bind_cols(ages = ages.fit, as.data.frame(cvModels1[[x]][[m]]$rates)),
                                           cols = 2:(length(years.fit) + 1), names_to = "year", values_to = "rate")%>% dplyr::select(rate)))))

  ytrain <- lapply(1:h,  function(x) ytrain0[[x]][[1]][-naRows[[x]],])

  data <- lapply(1:h, function(x) as.matrix(cbind(as.data.frame(xtrain[[x]]), rate =  ytrain[[x]])))

  result <- structure(list(metadata = data, cvmse = CVerror, models = names(models)))

  class(result) <- "stackmeta"

  return(result)

}
