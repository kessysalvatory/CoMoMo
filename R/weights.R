 #' @title Weights estimation using different methods.
#'
#' @description We consider four different model combination approaches:
#'  Simple Model Averaging, Bayesian Model Averaging, Model Confidence Set,
#'  and Stacked Regression Ensemble. These methods vary from each other depending on
#'  how they use the historical data to choose the combining weights or the models to be combined.
#'
#' @details Simple Model Averaging assigns equal weights to all the models.
#'
#' @details Bayesian Model Averaging estimates the weights using the posterior model probabilities.
#'
#' @details  Model Confidence Set chooses the subset of superior mortality models to combine
#' where each model is assigned equal weight
#'
#' @details  Stacked Regression Ensemble combines point forecasts from
#' multiple base learners using the weights that optimize a
#' cross-validation criterion.
#'
#' @param models is the specified models to be combined.
#'
#' @param h number of years for forecasting horizon.
#'
#' @param data an optional object of type StMoMoData containing
#'
#' information on deaths and exposures to be used for training the model.
#'
#' This is typically created with function \code{\link{StMoMoData}}.
#'
#' If this is not provided then the training data is taken from
#'
#' arguments, \code{Dxt}, \code{Ext}, \code{ages}, \code{years}.
#'
#' @param Dxt optional matrix of deaths data.
#'
#' @param Ext optional matrix of observed exposures of the same
#' dimension of \code{Dxt}.
#'
#' @param ages optional vector of ages corresponding to rows of
#' \code{Dxt} and \code{Ext}.
#'
#' @param years optional vector of years corresponding to rows of
#' \code{Dxt} and \code{Ext}.
#'
#' @param ages.fit optional vector of ages to include in the
#' training. Must be a subset of \code{ages}.
#'
#' @param years.fit optional vector of years to include in the
#' training. Must be a subset of \code{years}.
#'
#' @param holdout optional scalar of number of years withheld to calculate the projection bias.
#' Must be a subset of \code{years}.
#'
#' @param method optional paramater specifying how the models are trained. Cross-validation (cv) is default option
#' Single-validation (sv) option can also be specified.
#'
#' @return Returns an object of class \code{bma} with the following components:
#'
#' \item{Weights}{Returns the combination weights for different horizon.}
#'
#' \item{Method}{Returns the trainining method either cv or sv.}
#'
#'  @examples weights <- bma(modlist, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 5, method = "cv")

#' @export

valloss <- function(models, data = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL, ages = NULL, years = NULL, holdout = NULL, h = NULL)

{
  # Check the forecast horizon

  if (h>0 && h!= as.integer(h)) stop("The forecast horizon h must be a positive integer.")

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

  fit <- lapply(models, function(x) forecast(StMoMo::fit(x, data = data, ages.fit = ages.fit, years.fit = years.fit[1]:max(years.fit) - holdout), h = holdout, Dxt = Dxt, Ext = Ext, ages = ages, years = years))

  rates <- lapply(fit, function(x){tidyr::pivot_longer(dplyr::bind_cols(ages = ages.fit, as.data.frame(x$rates)), cols = 2:(holdout + 1), names_to = "year",
                                                       values_to = "rate")%>%dplyr::mutate(year = as.numeric(year), h = (holdout + 1))})

  rates_observed <- tidyr::pivot_longer(dplyr::bind_cols(ages = data$ages, as.data.frame(data$Dxt/data$Ext)),cols = 2:(length(data$years) + 1),
                                        names_to = "year", values_to = "obsrate")%>%dplyr::mutate(year = as.numeric(year))

  output <-list(rates = rates, obsrate = rates_observed)

}

#' @export

cvloss <- function(models, data = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL,
                   ages = NULL, years = NULL, h = NULL)

{

  # Check the forecast horizon

  if (h>0 && h!= as.integer(h)) stop("The forecast horizon h must be a positive integer.")

  # Check if more than one model is supplied

  if(length(models)<2) stop("Argument models needs to contain more than one model.")

  # names of specified models

  Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

  if (h>0 && h!= as.integer(h)) stop("This is invalid horizon")

  cvModels1 <- lapply(1:h, function(h) lapply(models, function(x) StMoMo::cv.StMoMo(x, h = h, data = data, ages.train = ages.fit, years.train = years.fit, returnY = TRUE, Dxt = Dxt, Ext = Ext, ages = ages, years = years)))

  cvRates1 <- lapply(cvModels1, function(m) lapply(m, function (x) {tidyr::pivot_longer(dplyr::bind_cols(ages = ages.fit, as.data.frame(x$cv.rates)), cols = 2:(length(years.fit) + 1),
                                                                                        names_to = "year", values_to = "rate")%>% dplyr::select(rate)}))

  xtrain <- lapply(1:length(1:h), function(x) log(as.matrix(dplyr::bind_cols(cvRates1[[x]]))))

  for (i in 1:length(1:h))

  { colnames(xtrain[[i]]) <- unlist(Specified_Models)

  }

  naRows <- lapply(1:length(1:h), function(x) which(is.na(rowSums(xtrain[[x]]))))

  cvRates_mcs <- lapply(1:length(1:h), function(x) dplyr::bind_cols(cvRates1[[x]]))

  # Cvmse

  cvmse <- lapply(1:length(1:h), function(h) lapply(1:length(models), function(x) cvModels1[[h]][[x]]["cv.mse"]))

  cverr<-  lapply(1:length(1:h), function(x) data.frame(h = x, cv = dplyr::bind_rows(cvmse[[x]])))

  CVerror <- dplyr::bind_rows(lapply(cverr, function(x) x %>%dplyr::mutate(model = unlist(Specified_Models))))

  for (i in 1:length(1:h))

  { colnames(cvRates_mcs[[i]]) <- unlist(Specified_Models)

  }

  Forecasts_mcs <- lapply(1:length(1:h), function(x)  cvRates_mcs[[x]][-naRows[[x]],])

  ytrain <- lapply(1:length(1:h), function(x) lapply(1:length(models), function(m) tidyr::pivot_longer(dplyr::bind_cols(ages = ages.fit, as.data.frame(cvModels1[[x]][[m]]$rates)), cols = 2:(length(years.fit) + 1),
                                                                                                       names_to = "year", values_to = "rate")%>% dplyr::select(rate)))

  observedrates_mcs <- lapply(1:length(1:h), function(x) lapply(1:length(models), function(m) ytrain[[x]][[m]][-naRows[[x]],]))

  obsRatesDF4 <- lapply(1:length(1:h), function(x)  dplyr::bind_cols(rep(observedrates_mcs[[x]][[1]], ncol(Forecasts_mcs[[1]]))))

  for (i in 1:length(1:h))

  { colnames(obsRatesDF4[[i]]) <- unlist(Specified_Models)

  }

  cvloss <- lapply(1:length(1:h), function(x)  log(Forecasts_mcs[[x]]) - log(obsRatesDF4[[x]]))

  out <- list(cvloss = cvloss, CVE =  CVerror, cvModels = cvModels1, cvmse = cvmse)

}


#' @export

bma <- function(models, method = "cv", data = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL, ages = NULL, years = NULL, holdout = round(length(years.fit)/3,0), h = NULL)

{

  # Check if more than one model is supplied

  if(length(models)<2) stop("Argument models needs to contain more than one model.")

  # check if the supplied models are different

  if( length(unique(unname(models)))==1) stop("Models must be different.")

  if (h>0 && h!= as.integer(h)) stop("The forecast horizon h must be a positive integer.")

  Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

  if ( holdout < round(length(years.fit)/3,0))
  {
    cat(paste("Warning: holdout is small \n"))
  }


  if (method!="cv" && method!="sv") stop("This is undefined method")

  if (method == "cv")

  {

    output0 <- cvloss(models = models, data = data, ages.fit = ages.fit, years.fit = years.fit, h = h, Dxt = Dxt, Ext = Ext, ages = ages, years = years)

    weights_bma<-  lapply(1:length(1:h), function(x) data.frame(h = x, weights = exp(-0.5*dplyr::bind_rows(output0$cvmse[[x]]))/sum(exp(-0.5*dplyr::bind_rows(output0$cvmse[[x]])))))

    out <- dplyr::bind_rows(lapply(weights_bma, function(x) x %>%dplyr::mutate(model = unlist(Specified_Models))))

    output <- out%>%dplyr::rename(weights = cv.mse)

    result <- list(weights = output, method = "cv", cvmse =  output0$CVE)

    class(result) <- "bma"

    return(result)
  }


  else if (method == "sv")

  {

    # names of specified models

    Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

    rates <- valloss(models = models, data = data, ages.fit = ages.fit, years.fit = years.fit, holdout = holdout, h = h,Dxt = Dxt, Ext = Ext, ages = ages, years = years)$rates

    rates_bma <- dplyr::bind_rows(lapply(1:length(models), function(x) dplyr::mutate(rates[[x]], model = Specified_Models[[x]])))

    rates_observed <- valloss(models = models, data = data, ages.fit = ages.fit, years.fit = years.fit, holdout = holdout, h = h, Dxt = Dxt, Ext = Ext, ages = ages, years = years)$obsrate

    rates_all <-dplyr::left_join(rates_bma, rates_observed)

    projection_errors <- rates_all%>%dplyr::mutate(bias = (log(rate) - log(obsrate)))%>%dplyr::group_by(model)%>%dplyr::summarise(mae = mean(bias, na.rm = TRUE))%>%dplyr::ungroup()%>%dplyr::select(mae)

    weights_bma <- exp(-abs(projection_errors))/sum(exp(-abs(projection_errors)))

    output <- (dplyr::bind_rows(lapply(rep(list(weights_bma), h), function(x) x %>%dplyr::mutate(model = Specified_Models)))%>%dplyr::mutate(h = rep(1:h, each=length(Specified_Models))))[, c("h","mae","model")]

    colnames(output) <- c("h","weights", "model")

    output$model <- factor(output$model, levels = Specified_Models, labels = Specified_Models)

    res <- list(weights = output, method = "sv")

    class(res) <- "bma"

    return(res)

  }
}

#' @export
#'
#' @param metalearner a supervised machine learning algorithm for learning the weights in the stacked regression ensemble.
#'
#' Default option is non-negative least square regression (nnls)
#'
#' Other options are standard linear regression (Linear), lasso regression (Lasso)
#'
#' Ridge regression (Ridge) and Elastic net (Elastic)
#'
#' @return Returns an object of class \code{stack} with the following components:
#'
#' \item{Weights}{Returns the combination weights for different horizon.}
#'
#' \item{Method}{Returns the meta-learner used to learn the weights.}
#'
#' @examples weights <- stack(modlist, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 5, metalearner = "nnls")


stack <- function(models, data = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL, ages = NULL, years = NULL, h = NULL, metalearner = "nnls", normalize = TRUE, saveMetadata =  TRUE)

{

  # check if xtrain and ytrain exist

  h.prev <- 0

  horizon.exist <-  FALSE

  if ( file.exists("horizon.rda"))

  {

    horizon.exist <- TRUE

    load("horizon.rda")

  }

  ageLength <- 0

  ageLength.exist <-  FALSE

  if ( file.exists("ageLength.rda"))

  {
    ageLength.exist <-  TRUE

    load("ageLength.rda")

  }

  yearLength <- 0

  yearLength.exist <-  FALSE

  if ( file.exists("yearLength.rda"))

  {
    yearLength.exist <-  TRUE

    load("yearLength.rda")

  }

  Dxtlength <- 0

  Dxtlength.exist <-  FALSE

  if ( file.exists("Dxtlength.rda"))

  {
    Dxtlength.exist <-  TRUE

    load("Dxtlength.rda")

  }

  Extlength <- 0

  Extlength.exist <-  FALSE

  if ( file.exists("Extlength.rda"))

  {
    Extlength.exist <-  TRUE

    load("Extlength.rda")

  }


  agesLength <- 0

  agesLength.exist <-  FALSE

  if ( file.exists("agesLength.rda"))

  {
    agesLength.exist <-  TRUE

    load("agesLength.rda")

  }

  yearsLength <- 0

  yearsLength.exist <-  FALSE

  if ( file.exists("yearsLength.rda"))

  {
    yearsLength.exist <-  TRUE

    load("yearsLength.rda")

  }

  Datalabel <- 0

  Datalabel.exist <-  FALSE

  if ( file.exists("Datalabel.rda"))

  {
    Datalabel.exist <-  TRUE

    load("Datalabel.rda")

  }

  Dataseries <- 0

  Dataseries.exist <-  FALSE

  if ( file.exists("Dataseries.rda"))

  {
    Dataseries.exist <-  TRUE

    load("Dataseries.rda")

  }

  prevModels <- 0

  prevModels.exist <-  FALSE

  if ( file.exists("prevModels.rda"))

  {
    prevModels.exist <-  TRUE

    load("prevModels.rda")

  }

  int <- intersect(prevModels, models)

  if( file.exists("data.rda") && file.exists("data0.rda") && file.exists("naRows0.rda") && file.exists("CVerror.rda") && horizon.exist && h.prev == h && ageLength.exist && yearLength.exist

      && yearLength == length(years.fit) && ageLength == length(ages.fit) && Datalabel.exist && Datalabel == data$label && Dataseries.exist && Dataseries == data$series && Dxtlength.exist && Dxtlength == length(data$Dxt) && Extlength.exist&&
      Extlength == length(data$Ext) &&  yearsLength.exist &&  agesLength == length(data$ages) && yearsLength == length(data$years) && agesLength.exist && prevModels.exist && length(int) == length(models) &&  length(int) == length(prevModels)) {

    # load data

    load("data.rda")

    load("data0.rda")

    load("naRows0.rda")

    load("CVerror.rda")

  } else {


    # Check the forecast horizon

    if (h>0 && h!= as.integer(h)) stop("The forecast horizon h must be a positive integer.")

    h.prev <- h

    prevModels <- models

    ageLength <- length(ages.fit)

    yearLength <- length(years.fit)

    Datalabel <- data$label

    Dataseries <- data$series

    Dxtlength <- length(data$Dxt)

    Extlength <- length(data$Ext)

    agesLength <- length(data$ages)

    yearsLength <- length(data$years)

    # checking the horizons that weights can be generated

    if (h > (max(years.fit) - years.fit[1]-1)) {
      stop("The time series must be longer than h.")
    }

    # Check if more than one model is supplied

    if(length(models)<2) stop("Argument models needs to contain more than one model.")

    # check if the supplied models are different

    if( length(unique(unname(models)))==1) stop("Models must be different.")

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

    if (metalearner!="Lasso" && metalearner!="Ridge" && metalearner!="Elastic" && metalearner!="nnls" && metalearner!="Linear") stop("unknown metalearner.")

    cvModels1 <- lapply(1:h, function(h) lapply(models, function(x) StMoMo::cv.StMoMo(x, h = h, data = data, ages.train = ages.fit, years.train = years.fit, returnY = TRUE, Dxt = Dxt, Ext = Ext, ages = ages, years = years)))

    cvRates1 <- lapply(cvModels1, function(m) lapply(m, function (x) {tidyr::pivot_longer(dplyr::bind_cols(ages = ages.fit, as.data.frame(x$cv.rates)), cols = 2:(length(years.fit) + 1),
                                                                                          names_to = "year", values_to = "rate")%>% dplyr::select(rate)}))

    # for stacked regression

    xtrain <- lapply(1:length(1:h), function(x) log(as.matrix(dplyr::bind_cols(cvRates1[[x]]))))

    for (i in 1:length(1:h))

    { colnames(xtrain[[i]]) <- unlist(Specified_Models)

    }


    naRows <- lapply(1:length(1:h), function(x) which(is.na(rowSums(xtrain[[x]]))))


    ytrain <- lapply(1:length(1:h), function(x) lapply(1:length(models), function(m) log(as.matrix(tidyr::pivot_longer(dplyr::bind_cols(ages = ages.fit, as.data.frame(cvModels1[[x]][[m]]$rates)), cols = 2:(length(years.fit) + 1),
                                                                                                                       names_to = "year", values_to = "rate")%>% dplyr::select(rate)))))

    # Cvmse

    CV <-lapply(1:length(1:h), function(h) lapply(1:length(models), function(x) cvModels1[[h]][[x]]["cv.mse"]))

    cverr<-  lapply(1:length(1:h), function(x) data.frame(h = x, cv = dplyr::bind_rows(CV[[x]])))

    CVerror <- dplyr::bind_rows(lapply(cverr, function(x) x %>%dplyr::mutate(model = unlist(Specified_Models))))


    if (saveMetadata)

    {

      save(naRows, file = "naRows0.rda")

      save(CVerror, file = "CVerror.rda")

      save(ytrain, file = "data0.rda")

      save(h.prev, file = "horizon.rda")

      save(xtrain, file = "data.rda")

      save(ageLength, file = "ageLength.rda")

      save(yearLength, file = "yearLength.rda")

      save(Datalabel, file = "Datalabel.rda")

      save(Dataseries, file = "Dataseries.rda")

      save(Extlength, file = "Extlength.rda")

      save(Dxtlength, file = "Dxtlength.rda")

      save(agesLength, file = "agesLength.rda")

      save(yearsLength, file = "yearsLength.rda")

      save(prevModels, file = "prevModels.rda")

    }

    else if (!saveMetadata)

    {

      if (file.exists("prevModels.rda")) file.remove("prevModels.rda")

      if (file.exists("data.rda")) file.remove("data.rda")

      if (file.exists("data0.rda")) file.remove("data0.rda")

      if (file.exists("naRows0.rda")) file.remove("naRows0.rda")

      if (file.exists("CVerror.rda")) file.remove("CVerror.rda")

      if (file.exists("horizon.rda")) file.remove("horizon.rda")

      if (file.exists("Dxtlength.rda")) file.remove("Dxtlength.rda")

      if (file.exists("Extlength.rda")) file.remove("Extlength.rda")

      if (file.exists("ageLength.rda")) file.remove("ageLength.rda")

      if (file.exists("yearLength.rda")) file.remove("yearLength.rda")

      if (file.exists("Datalabel.rda")) file.remove("Datalabel.rda")

      if (file.exists("Dataseries.rda")) file.remove("Dataseries.rda")

      if (file.exists("yearsLength.rda")) file.remove("yearsLength.rda")

      if (file.exists("agesLength.rda")) file.remove("agesLength.rda")

    }

  }

  # Train meta model

  if (normalize)

  {  # learning weights using ridge regression

    if (metalearner=="Ridge"){

      Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

      ridge.model <- lapply(1:length(1:h), function(x) glmnet::glmnet(xtrain[[x]][-naRows[[x]],], ytrain[[x]][[1]][-naRows[[x]],], alpha = 0))

      cv_ridge <- lapply(1:length(1:h), function(x) glmnet::cv.glmnet(xtrain[[x]][-naRows[[x]],], ytrain[[x]][[1]][-naRows[[x]],], alpha = 0))

      coeffients <- lapply(1:length(1:h), function(x) coef(ridge.model[[x]], s = cv_ridge[[x]]$lambda.min)[-1])

      weights_ridge <- lapply(1:length(1:h), function(x) as.matrix(coeffients[[x]]/sum(coeffients[[x]])))

      weights1_ridge <-  lapply(1:length(1:h), function(x) data.frame(h = x, weights = weights_ridge[[x]]))

      weightsDF_ridge <- dplyr::bind_rows(lapply(weights1_ridge, function(x) x %>%dplyr::mutate(model = unlist(Specified_Models))))

      result <- list(weights =  weightsDF_ridge, metalearner = "Ridge", cvmse =  CVerror)

      class(result) <- "stack"

      return(result)

    }

    # learning weights using lasso regression

    else if (metalearner=="Lasso"){

      Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

      lasso.model <- lapply(1:length(1:h), function(x) glmnet::glmnet(xtrain[[x]][-naRows[[x]],], ytrain[[x]][[1]][-naRows[[x]],], alpha = 1))

      cv_lasso <- lapply(1:length(1:h), function(x) glmnet::cv.glmnet(xtrain[[x]][-naRows[[x]],], ytrain[[x]][[1]][-naRows[[x]],], alpha = 1))

      coeffients <- lapply(1:length(1:h), function(x) coef(lasso.model[[x]], s = cv_lasso[[x]]$lambda.min)[-1])

      weights_lasso <- lapply(1:length(1:h), function(x) as.matrix(coeffients[[x]]/sum(coeffients[[x]])))

      weights1_lasso <-  lapply(1:length(1:h), function(x) data.frame(h = x, weights = weights_lasso[[x]]))

      weightsDF_lasso <- dplyr::bind_rows(lapply(weights1_lasso, function(x) x %>% dplyr::mutate(model = unlist(Specified_Models))))

      result <- list(weights =  weightsDF_lasso, metalearner = "Lasso", cvmse =  CVerror)

      class(result) <- "stack"

      return(result)

    }

    # learning weights using elastic net regression

    else if (metalearner=="Elastic"){

      Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

      elastic.model <- lapply(1:length(1:h), function(x) glmnet::glmnet(xtrain[[x]][-naRows[[x]],], ytrain[[x]][[1]][-naRows[[x]],], alpha = 0.5))

      cv_elnet <- lapply(1:length(1:h), function(x) glmnet::cv.glmnet(xtrain[[x]][-naRows[[x]],], ytrain[[x]][[1]][-naRows[[x]],], alpha = 0.5))

      coeffients <- lapply(1:length(1:h), function(x) coef(elastic.model[[x]], s = cv_elnet[[x]]$lambda.min)[-1])

      weights_elastic <- lapply(1:length(1:h), function(x) as.matrix(coeffients[[x]]/sum(coeffients[[x]])))

      weights1_elastic <-  lapply(1:length(1:h), function(x) data.frame(h = x, weights = weights_elastic[[x]]))

      weightsDF_elastic <- dplyr::bind_rows(lapply(weights1_elastic, function(x) x %>% dplyr::mutate(model = unlist(Specified_Models))))

      result <- list(weights =  weightsDF_elastic, metalearner = "Elastic",  cvmse =  CVerror)

      class(result) <- "stack"

      return(result)

    }

    # learning weights using non-negative least square regression

    else if(metalearner=="nnls") {

      Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

      coeffients <-lapply(1:length(1:h), function(x) nnls::nnls(xtrain[[x]][-naRows[[x]],], ytrain[[x]][[1]][-naRows[[x]],])$x)

      weights_nnls <- lapply(1:length(1:h), function(x) as.matrix(coeffients[[x]]/sum(coeffients[[x]])))

      weights1_nnls <-  lapply(1:length(1:h), function(x) data.frame(h = x, weights = weights_nnls[[x]]))

      weightsDF_nnls <- dplyr::bind_rows(lapply(weights1_nnls, function(x) x %>% dplyr::mutate(model = unlist(Specified_Models))))

      result <- list(weights = weightsDF_nnls, metalearner = "nnls", cvmse =  CVerror)

      class(result) <- "stack"

      return(result)

    }

    # learning weights using linear squared regression

    else if (metalearner=="Linear"){

      Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

      y <- lapply(1:length(1:h), function(x) lapply(1:length(models), function(m) ytrain[[x]][[m]][-naRows[[x]],]))


      data_train <- lapply(1:length(1:h), function(x) data.frame(y = y[[x]][[1]], xtrain[[x]][-naRows[[x]],]))

      linear.model <- lapply(1:length(1:h), function(x) lm(y~., data = data_train[[x]]))

      coeffients <- lapply(1:length(1:h), function(x) coef(linear.model[[x]])[-1])

      weights_linear <- lapply(1:length(1:h), function(x) as.matrix(coeffients[[x]]/sum(coeffients[[x]])))

      weights1_linear <- lapply(1:length(1:h), function(x) data.frame(h = x, weights = weights_linear[[x]]))

      weightsDF_linear <- dplyr::bind_rows(lapply(weights1_linear, function(x) x %>% dplyr::mutate(model = unlist(Specified_Models))))

      result <- list(weights = weightsDF_linear,  metalearner = "Linear",  cvmse =  CVerror)

      class(result) <- "stack"

      return(result)

    }

  }

  else if (!normalize)

  {

    if (metalearner=="Ridge"){

      # names of specified models

      Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

      ridge.model <- lapply(1:length(1:h), function(x) glmnet::glmnet(xtrain[[x]][-naRows[[x]],], ytrain[[x]][[1]][-naRows[[x]],], alpha = 0))

      cv_ridge <- lapply(1:length(1:h), function(x) glmnet::cv.glmnet(xtrain[[x]][-naRows[[x]],], ytrain[[x]][[1]][-naRows[[x]],], alpha = 0))

      coeffients <- lapply(1:length(1:h), function(x) coef(ridge.model[[x]], s = cv_ridge[[x]]$lambda.min)[-1])

      weights_ridge <- lapply(1:length(1:h), function(x) as.matrix(coeffients[[x]]))

      weights1_ridge <-  lapply(1:length(1:h), function(x) data.frame(h = x, weights = weights_ridge[[x]]))

      weightsDF_ridge <- dplyr::bind_rows(lapply(weights1_ridge, function(x) x %>%dplyr::mutate(model = unlist(Specified_Models))))

      result <- list(weights =  weightsDF_ridge, metalearner = "Ridge",  cvmse =  CVerror)

      class(result) <- "stack"

      return(result)

    }

    else if (metalearner=="Lasso"){

      # names of specified models

      Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

      lasso.model <- lapply(1:length(1:h), function(x) glmnet::glmnet(xtrain[[x]][-naRows[[x]],], ytrain[[x]][[1]][-naRows[[x]],], alpha = 1))

      cv_lasso <- lapply(1:length(1:h), function(x) glmnet::cv.glmnet(xtrain[[x]][-naRows[[x]],], ytrain[[x]][[1]][-naRows[[x]],], alpha = 1))

      coeffients <- lapply(1:length(1:h), function(x) coef(lasso.model[[x]], s = cv_lasso[[x]]$lambda.min)[-1])

      weights_lasso <- lapply(1:length(1:h), function(x) as.matrix(coeffients[[x]]))

      weights1_lasso <-  lapply(1:length(1:h), function(x) data.frame(h = x, weights = weights_lasso[[x]]))

      weightsDF_lasso <- dplyr::bind_rows(lapply(weights1_lasso, function(x) x %>% dplyr::mutate(model = unlist(Specified_Models))))

      result <- list(weights =  weightsDF_lasso, metalearner = "Lasso", cvmse =  CVerror)

      class(result) <- "stack"

      return(result)

    }

    else if (metalearner=="Elastic"){

      # names of specified models

      Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

      elastic.model <- lapply(1:length(1:h), function(x) glmnet::glmnet(xtrain[[x]][-naRows[[x]],], ytrain[[x]][[1]][-naRows[[x]],], alpha = 0.5))

      cv_elnet <- lapply(1:length(1:h), function(x) glmnet::cv.glmnet(xtrain[[x]][-naRows[[x]],], ytrain[[x]][[1]][-naRows[[x]],], alpha = 0.5))

      coeffients <- lapply(1:length(1:h), function(x) coef(elastic.model[[x]], s = cv_elnet[[x]]$lambda.min)[-1])

      weights_elastic <-lapply(1:length(1:h), function(x) as.matrix(coeffients[[x]]))

      weights1_elastic <-  lapply(1:length(1:h), function(x) data.frame(h = x, weights = weights_elastic[[x]]))

      weightsDF_elastic <- dplyr::bind_rows(lapply(weights1_elastic, function(x) x %>% dplyr::mutate(model = unlist(Specified_Models))))

      result <- list(weights =  weightsDF_elastic, metalearner = "Elastic", cvmse =  CVerror)

      class(result) <- "stack"

      return(result)
    }

    else if(metalearner=="nnls") {

      # names of specified models

      Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

      coeffients <-lapply(1:length(1:h), function(x) nnls::nnls(xtrain[[x]][-naRows[[x]],], ytrain[[x]][[1]][-naRows[[x]],])$x)

      weights_nnls <- lapply(1:length(1:h), function(x) as.matrix(coeffients[[x]]))

      weights1_nnls <-  lapply(1:length(1:h), function(x) data.frame(h = x, weights = weights_nnls[[x]]))

      weightsDF_nnls <- dplyr::bind_rows(lapply(weights1_nnls, function(x) x %>% dplyr::mutate(model = unlist(Specified_Models))))

      result <- list(weights = weightsDF_nnls, metalearner = "nnls", cvmse =  CVerror)

      class(result) <- "stack"

      return(result)

    }

    else if (metalearner=="Linear"){

      # names of specified models

      Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

      linear.model <- lapply(1:length(1:h), function(x) lm(y~., data = data_train[[x]]))

      coeffients <- lapply(1:length(1:h), function(x) coef(linear.model[[x]])[-1])

      weights_linear <- lapply(1:length(1:h), function(x) as.matrix(coeffients[[x]]))

      weights1_linear <- lapply(1:length(1:h), function(x) data.frame(h = x, weights = weights_linear[[x]]))

      weightsDF_linear <- dplyr::bind_rows(lapply(weights1_linear, function(x) x %>% dplyr::mutate(model = unlist(Specified_Models))))

      result <- list(weights = weightsDF_linear, metalearner = "Linear",  cvmse =  CVerror)

      class(result) <- "stack"

      return(result)

    }

  }

}

#' @export
#'
#' @param alpha is the level of the test
#'
#' @param B is the bootstrap samples.
#'
#' @param l is the block length.
#'
#' @return Returns an object of class \code{mcs} with the following components:
#'
#' \item{Weights}{Returns the combination weights for different horizon.}
#'
#' \item{Method}{Returns the trainining method either cv or sv.}
#'
#' \item{Selected models}{Returns the superior models selected.}
#'
#'  @examples weights <- mcs(modlist, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 5, B = 5000, l=3, alpha = 0.1,  method = "cv")
#'

mc <- function(Loss, B, l){

  # Matrix of boostrap indixes

  if (ncol(Loss) == 1) stop("MCS: Need more data. Only one model entered.")
  if(any(is.na(Loss))) stop("NAs in Loss are not allowed")
  if(any(abs(Loss) == Inf)) stop("Inf in Loss are not allowed")

  LbarB <- boot::tsboot(Loss, colMeans, R = B, sim = "fixed",l = l)$t
  Lbar <- colMeans(Loss)
  zeta.Bi <- t(t(LbarB)-Lbar)
  save.res <- c()

  for(j in 1:(ncol(Loss) - 1)){
    Lbardot <- mean(Lbar)
    zetadot <- rowMeans(zeta.Bi)
    vard <- colMeans((zeta.Bi-zetadot)^2)
    t.stat <- (Lbar-Lbardot)/sqrt(vard)
    t.max <- max(t.stat)
    model.t.max <- which(t.stat==t.max)
    t.stat.b <- t(t(zeta.Bi-zetadot)/sqrt(vard))
    t.max.b <- apply(t.stat.b,1,max)
    p <- length(which(t.max < t.max.b))/B
    save.res <- c(save.res,p)
    names(save.res)[j] <- names(model.t.max)
    Lbar <- Lbar[-model.t.max]
    zeta.Bi <- zeta.Bi[,-model.t.max]
  }

  save.res <- c(save.res,1)
  names(save.res)[j+1] <- names(Lbar)
  save.p <- save.res

  for(i in 2:(length(save.res)-1)){
    save.p[i] <- max(save.res[i-1],save.res[i])
  }

  aux <- match(colnames(Loss),names(save.p))
  save.p <- save.p[aux]
  save.res <- save.res[aux]

  return(list(test = save.p, individual = save.res))
}

#' @export

mcs <- function(models,  method = "cv", data = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL,
                ages = NULL, years = NULL, h = NULL, B = 5000, l = 3, alpha = 0.1, holdout = round(length(years.fit)/3,0))

{

  # Check if more than one model is supplied

  if(length(models)<2) stop("Argument models needs to contain more than one model.")

  # check if the supplied models are different

  if( length(unique(unname(models)))==1) stop("Models must be different.")

  # Check the forecast horizon

  if (h>0 && h!= as.integer(h)) stop("The forecast horizon h must be a positive integer.")

  if ( holdout < round(length(years.fit)/3,0))
  {
    cat(paste("Warning: holdout is small \n"))
  }

  if (B < 1000)
  {
    cat(paste("Warning: B is small \n"))
  }

  modelNames <- names(models)

  if (method!="cv" && method!="sv") stop("This is undefined method")

  if (method == "sv")

  {
    # loss

    forecasts <-valloss(models = models, data = data, ages.fit = ages.fit, years.fit = years.fit, holdout = holdout, h = h, Dxt = Dxt, Ext = Ext, ages = ages, years = years)$rates

    modelforecasts <- dplyr::bind_cols(lapply(1:length(forecasts), function(x) forecasts[[x]]["rate"]))

    colnames(modelforecasts) <- names(models)

    obsRatesDF <-valloss(models = models, data = data, ages.fit = ages.fit, years.fit = years.fit, holdout = holdout, h = h, Dxt = Dxt, Ext = Ext, ages = ages, years = years)$obsrate


    observedrates <- dplyr::filter(obsRatesDF, ages%in%ages.fit, year%in%((max(years.fit) - holdout + 1):max(years.fit)))%>% dplyr::select(obsrate)

    obsRatesDF0 <- dplyr::bind_cols(rep(observedrates, ncol(modelforecasts)))

    colnames(obsRatesDF0) <- names(models)

    # weights

    MCS_vald <- mc((log(modelforecasts)-log(obsRatesDF0))^2, B = B, l = l)$test

    selected <- which(MCS_vald>alpha)

    weights <- c()

    for (model in modelNames) {
      if (model %in% modelNames[selected]) {weights = c(weights, 1/length(modelNames[selected])) } else {weights = c(weights, 0)}

    }

    out <- (dplyr::bind_rows(lapply(rep(list(as.data.frame(weights)), h), function(x) x%>%dplyr::mutate(model = modelNames)))%>%dplyr::mutate(h = rep(1:h, each = length(models))))[,c("h","weights","model")]

    res <- list(weights = out,  selected =  selected)

    class(res) <- "mcs"

    return(res)

  }

  else if ( method=="cv")

  {
    output0 <- cvloss(models = models, data = data, ages.fit = ages.fit, years.fit = years.fit, h = h, Dxt = Dxt, Ext = Ext, ages = ages, years = years)

    MCS_cv <- lapply(1:length(1:h), function(x) mc((output0$cvloss[[x]])^2,  B = B, l = l)$test)

    selected <- lapply(MCS_cv, function(x)  which(x > alpha))

    weights <- c()

    for (i in 1:length(1:h))

    {
      for (model in modelNames) {
        if (model %in% modelNames[selected[[i]]]) {weights = c(weights, 1/length(modelNames[selected[[i]]])) } else {weights = c(weights,0)

        }
      }
    }

    output <- data.frame(h = rep(1:length(1:h), each = length(modelNames)), weights, model = rep(unlist(modelNames)), stringsAsFactors = F)

    result <- list(weights = output, method = "cv", selected =  selected,  cvmse =  output0$CVE)

    class(result) <- "mcs"

    return(result)

  }

}
