#' @title fitting different methods.
#'
#' @export

fitCoMoMo <- function(models, data = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL, ages = NULL, years = NULL)

{
  # Check if more than one model is supplied

  if(length(models)<2) stop("Argument models needs to contain more than one model.")

  # check if the supplied models are different

  if(length(unique(unname(models))) < length(models)) stop("Models must be different.")

  # Determine model.fitsting ages and years are specified

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
