#' @title Weights estimation using different model combinations methods.
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
#' @details Model Confidence Set chooses the subset of superior mortality models to combine
#' where each model is assigned equal weight
#'
#' @details Stacked Regression Ensemble combines point forecasts from
#' multiple base learners using the weights that optimize a
#' cross-validation criterion.
#'
#' @param models are the specified list of models to be combined.
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
#' @export
#'
valloss <- function(models, data = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL, ages = NULL, years = NULL, holdout = NULL, h = NULL)

{
  # Check the forecast horizon

  if (h<0 || h!= as.integer(h)) stop("The forecast horizon h must be a positive integer.")

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

  output <-structure(list(rates = rates, obsrate = rates_observed))

}

#' @export

cvloss <- function(models, data = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL,
                   ages = NULL, years = NULL, h = NULL)

{
  # Check the forecast horizon

  if (h<0 || h!= as.integer(h)) stop("The forecast horizon h must be a positive integer.")

  # Check if more than one model is supplied

  if(length(models)<2) stop("Argument models needs to contain more than one model.")

  # models must be different

  if(length(unique(unname(models))) < length(models)) stop("Models must be different.")

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
    if (is.null(ages)) ages <- 1:nrow(Dxt)
    if (is.null(years)) years <- 1:ncol(Dxt)
  }

  if (is.null(ages.fit)) ages.fit <- ages
  if (is.null(years.fit)) years.fit <- years

  # names of specified models

 Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

 cvModels1 <- lapply(1:h, function(h) lapply(models, function(x) StMoMo::cv.StMoMo(x, h = h, data = data, ages.train = ages.fit, years.train = years.fit, returnY = TRUE, Dxt = Dxt, Ext = Ext, ages = ages, years = years)))

 cvRates1 <- lapply(cvModels1, function(m) lapply(m, function (x) {tidyr::pivot_longer(dplyr::bind_cols(ages = ages.fit, as.data.frame(x$cv.rates)), cols = 2:(length(years.fit) + 1),
                                                                                            names_to = "year", values_to = "rate")%>% dplyr::select(rate)}))

 xtrain <- lapply(1:h, function(x) log(as.matrix(dplyr::bind_cols(cvRates1[[x]]))))

 for (i in 1:h)

 { colnames(xtrain[[i]]) <- unlist(Specified_Models)

 }

 naRows <- lapply(1:h, function(x) which(is.na(rowSums(xtrain[[x]]))))

 cvRates_mcs <- lapply(1:h, function(x) dplyr::bind_cols(cvRates1[[x]]))

 for (i in 1:h)

  { colnames(cvRates_mcs[[i]]) <- unlist(Specified_Models)

  }

  Forecasts_mcs <- lapply(1:h, function(x)  cvRates_mcs[[x]][-naRows[[x]],])

  ytrain <- lapply(1:h, function(x) lapply(1:length(models), function(m) tidyr::pivot_longer(dplyr::bind_cols(ages = ages.fit, as.data.frame(cvModels1[[x]][[m]]$rates)), cols = 2:(length(years.fit) + 1),
                                                                                            names_to = "year", values_to = "rate")%>% dplyr::select(rate)))

  observedrates_mcs <- lapply(1:h, function(x) lapply(1:length(models), function(m) ytrain[[x]][[m]][-naRows[[x]],]))

  obsRatesDF4 <- lapply(1:h, function(x)  dplyr::bind_cols(rep(observedrates_mcs[[x]][[1]], ncol(Forecasts_mcs[[1]]))))

  for (i in 1:h)

  { colnames(obsRatesDF4[[i]]) <- unlist(Specified_Models)

   }

  # Cvmse

  cvmse <- lapply(1:h, function(h) lapply(1:length(models), function(x) cvModels1[[h]][[x]]["cv.mse"]))

  cverr<-  lapply(1:h, function(x) data.frame(h = x, cv = dplyr::bind_rows(cvmse[[x]])))

  CVerror0 <- dplyr::bind_rows(lapply(cverr, function(x) x %>%dplyr::mutate(model = unlist(Specified_Models))))

  CVerror <-  CVerror0 %>% dplyr::rename(model.cvmse ="cv.mse")

  cvloss <- lapply(1:h, function(x)  log(Forecasts_mcs[[x]]) - log(obsRatesDF4[[x]]))

  out <- structure(list(cvloss = cvloss, CVE =  CVerror, cvModels = cvModels1, cvmse = cvmse, cvRates = cvRates1,  xtrain =  xtrain, naRows = naRows, ytrain = ytrain))

}

#' @param holdout optional scalar of number of years withheld to calculate the projection bias.
#' Must be a subset of \code{years}.
#'
#' @param method optional paramater specifying how the models are trained. Cross-validation (cv) is default option
#' Single-validation (sv) option can also be specified.
#'
#' @return Returns an object of class \code{weight} with the following components:
#'
#' \item{weights}{Returns the combination weights for different horizons.}
#' \item{comb.method}{Returns the combination method}
#' \item{cvmse}{Returns the cross-validation mean squared error for each all the indiviual models for different horizons.}
#' \item{method}{Returns the trainining method either cv or sv.}
#'                 
#' @references
#'
#' Kessy, Salvatory, Michael Sherris, Andrés Villegas, and Jonathan Ziveyi. 2021. 
#' “Mortality Forecasting Using Stacked Regression Ensembles.” SSRN Electronic Journal. https://doi.org/10.2139/ssrn.3823511.
#'
#' @examples
#'
#' define models
#' Dont run 
#' LC <- lc()
#' APC <- apc()
#' modlist <- list("LC"= LC, "APC" = APC)
#' bma_weight_val <- bma(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit,h = 15, method = "sv")
#' bma_weight_cv <- bma(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 15, method = "cv")                  
#'
#' @export

bma <- function(models, method = "cv", data = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL, ages = NULL, years = NULL, holdout = round(length(years.fit)/3,0), h = NULL)

{ # the function calculate the Bayesian weights. 

  if ( holdout < round(length(years.fit)/3,0))
  {
    cat(paste("Warning: holdout is small \n"))
  }


  if (method!="cv" && method!="sv") stop("This is undefined method")

  Specified_Models <- lapply(1:length(models), function(x) names(models[x]))

  if (method == "cv")

  {

    output0 <- cvloss(models = models, data = data, ages.fit = ages.fit, years.fit = years.fit, h = h, Dxt = Dxt, Ext = Ext, ages = ages, years = years)

    weights_bma<-  lapply(1:h, function(x) data.frame(h = x, weights = exp(-0.5*dplyr::bind_rows(output0$cvmse[[x]]))/sum(exp(-0.5*dplyr::bind_rows(output0$cvmse[[x]])))))

    out <- dplyr::bind_rows(lapply(weights_bma, function(x) x %>%dplyr::mutate(model = unlist(Specified_Models))))

    output <- out%>%dplyr::rename(model.weights = cv.mse)

    result <- structure(list(weights = output, method = "cv", cvmse =  output0$CVE, comb.method = "BMA"))

    class(result)<- "weight"

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

    colnames(output) <- c("h","model.weights", "model")

    output$model <- factor(output$model, levels = Specified_Models, labels = Specified_Models)

    res <- structure(list(weights = output, method = "sv", comb.method = "BMA"))

    class(res) <- "weight"

    return(res)

  }
}

#'
#' @param metalearner a supervised machine learning algorithm for learning the weights in the stacked regression ensemble.
#' Default option is non-negative least square regression (nnls)
#' Other options are standard linear regression (Linear), lasso regression (Lasso)
#' Ridge regression (Ridge) and Elastic net (Elastic)
#'
#' @param normalize allows the user to specify if the weights are to sum to one or not. The default option normalize = TRUE
#' makes all the weights sum to a unit, otherwise when normalize = FALSE the weights do not sum to one.
#'
#' @return Returns an object of class \code{weight} with the following components:
#'                                    
#' \item{Weights}{Returns the combination weights for different horizons.}
#' \item{metalearner}{Returns the meta-learner used to learn the weights.}                                     
#' \item{comb.method}{Returns the combination method}
#'
#' @examples
#'
#' # define models
#' Dont run 
#' LC <- lc()
#' APC <- apc()
#' modlist <- list("LC"= LC, "APC" = APC)
#' metaData <- stackMetadata(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 15)                     
#' weight_stack <- stack(metaData, metalearner = "nnls", normalize = TRUE)                                      
#' @export
#'
#'
stack  <- function(stackmeta,...) {
  UseMethod("stack")
}

#'
#' @export
#'
stack.stackmeta <- function(stackmeta, metalearner = "nnls", normalize = TRUE)

{ # the function calculate the weights using stacked regression. 

  if (metalearner!="Lasso" && metalearner!="Ridge" && metalearner!="Elastic" && metalearner!="nnls" && metalearner!="Linear") stop("unknown metalearner.")
  
  # Train meta model

  if (normalize)
    
    # make all combining weights sum to a unit 

  {  # learning weights using ridge regression

    if (metalearner=="Ridge"){

      ridge.model <- lapply(1:length(stackmeta$metadata), function(x) glmnet::glmnet(stackmeta$metadata[[x]][,1:length(stackmeta$models)], stackmeta$metadata[[x]][,ncol(stackmeta$metadata[[x]])], alpha = 0))

      cv_ridge <- lapply(1:length(stackmeta$metadata), function(x) glmnet::cv.glmnet(stackmeta$metadata[[x]][,1:length(stackmeta$models)], stackmeta$metadata[[x]][,ncol(stackmeta$metadata[[x]])], alpha = 0))

      coeffients <- lapply(1:length(stackmeta$metadata), function(x) coef(ridge.model[[x]], s = cv_ridge[[x]]$lambda.min)[-1])

      weights_ridge <- lapply(1:length(stackmeta$metadata), function(x) as.matrix(coeffients[[x]]/sum(coeffients[[x]])))

      weights1_ridge <-  lapply(1:length(stackmeta$metadata), function(x) data.frame(h = x, model.weights = weights_ridge[[x]]))

      weightsDF_ridge <- dplyr::bind_rows(lapply(weights1_ridge, function(x) x %>%dplyr::mutate(model = stackmeta$models)))

      result <- structure(list(weights =  weightsDF_ridge, metalearner = "Ridge", comb.method = "stack"))

      class(result) <- "weight"

      return(result)

    }

    # learning weights using lasso regression

    else if (metalearner=="Lasso"){

      lasso.model <- lapply(1:length(stackmeta$metadata), function(x) glmnet::glmnet(stackmeta$metadata[[x]][,1:length(stackmeta$models)], stackmeta$metadata[[x]][,ncol(stackmeta$metadata[[x]])], alpha = 1))

      cv_lasso <- lapply(1:length(stackmeta$metadata), function(x) glmnet::cv.glmnet( stackmeta$metadata[[x]][,1:length(stackmeta$models)], stackmeta$metadata[[x]][,ncol(stackmeta$metadata[[x]])], alpha = 1))

      coeffients <- lapply(1:length(stackmeta$metadata), function(x) coef(lasso.model[[x]], s = cv_lasso[[x]]$lambda.min)[-1])

      weights_lasso <- lapply(1:length(stackmeta$metadata), function(x) as.matrix(coeffients[[x]]/sum(coeffients[[x]])))

      weights1_lasso <-  lapply(1:length(stackmeta$metadata), function(x) data.frame(h = x, model.weights = weights_lasso[[x]]))

      weightsDF_lasso <- dplyr::bind_rows(lapply(weights1_lasso, function(x) x %>% dplyr::mutate(model = stackmeta$models)))

      result <- structure(list(weights =  weightsDF_lasso, metalearner = "Lasso", comb.method = "stack"))

      class(result) <- "weight"

      return(result)

    }

    # learning weights using elastic net regression

    else if (metalearner=="Elastic"){

      elastic.model <- lapply(1:length(stackmeta$metadata), function(x) glmnet::glmnet(stackmeta$metadata[[x]][,1:length(stackmeta$models)], stackmeta$metadata[[x]][,ncol(stackmeta$metadata[[x]])], alpha = 0.5))

      cv_elnet <- lapply(1:length(stackmeta$metadata), function(x) glmnet::cv.glmnet(stackmeta$metadata[[x]][,1:length(stackmeta$models)], stackmeta$metadata[[x]][,ncol(stackmeta$metadata[[x]])], alpha = 0.5))

      coeffients <- lapply(1:length(stackmeta$metadata), function(x) coef(elastic.model[[x]], s = cv_elnet[[x]]$lambda.min)[-1])

      weights_elastic <- lapply(1:length(stackmeta$metadata), function(x) as.matrix(coeffients[[x]]/sum(coeffients[[x]])))

      weights1_elastic <-  lapply(1:length(stackmeta$metadata), function(x) data.frame(h = x, model.weights = weights_elastic[[x]]))

      weightsDF_elastic <- dplyr::bind_rows(lapply(weights1_elastic, function(x) x %>% dplyr::mutate(model = stackmeta$models)))

      result <- structure(list(weights =  weightsDF_elastic, metalearner = "Elastic", comb.method = "stack"))

      class(result) <- "weight"

      return(result)

    }

    # learning weights using non-negative least square regression

    else if(metalearner=="nnls") {

      coeffients <-lapply(1:length(stackmeta$metadata), function(x) nnls::nnls(stackmeta$metadata[[x]][,1:length(stackmeta$models)], stackmeta$metadata[[x]][,ncol(stackmeta$metadata[[x]])])$x)

      weights_nnls <- lapply(1:length(1:length(stackmeta$metadata)), function(x) as.matrix(coeffients[[x]]/sum(coeffients[[x]])))

      weights1_nnls <-  lapply(1:length(1:length(stackmeta$metadata)), function(x) data.frame(h = x,  model.weights = weights_nnls[[x]]))

      weightsDF_nnls <- dplyr::bind_rows(lapply(weights1_nnls, function(x) x %>% dplyr::mutate(model = stackmeta$models)))

      result <- structure(list(weights = weightsDF_nnls, metalearner = "nnls", comb.method = "stack"))

      class(result) <- "weight"

      return(result)

    }

    # learning weights using standard linear squared regression

    else if (metalearner=="Linear"){

      linear.model <- lapply(1:length(stackmeta$metadata), function(x) lm(rate~., data =  as.data.frame(stackmeta$metadata[[x]])))

      coeffients <- lapply(1:length(stackmeta$metadata), function(x) unname(coef(linear.model[[x]])[-1]))

      weights_linear <- lapply(1:length(stackmeta$metadata), function(x) as.matrix(coeffients[[x]]/sum(coeffients[[x]])))

      weights1_linear <- lapply(1:length(stackmeta$metadata), function(x) data.frame(h = x, model.weights = weights_linear[[x]]))

      weightsDF_linear <- dplyr::bind_rows(lapply(weights1_linear, function(x) x %>% dplyr::mutate(model = stackmeta$models)))

      result <- structure(list(weights = weightsDF_linear,  metalearner = "Linear", comb.method = "stack"))

      class(result) <- "weight"

      return(result)

    }

  }

  else if (!normalize)

  { # not necessary the combining weights sum to a unit 
    
    # learning weights using ridge regression

    if (metalearner=="Ridge"){

      ridge.model <- lapply(1:length(stackmeta$metadata), function(x) glmnet::glmnet(stackmeta$metadata[[x]][,1:length(stackmeta$models)], stackmeta$metadata[[x]][,ncol(stackmeta$metadata[[x]])], alpha = 0))

      cv_ridge <- lapply(1:length(stackmeta$metadata), function(x) glmnet::cv.glmnet(stackmeta$metadata[[x]][,1:length(stackmeta$models)], stackmeta$metadata[[x]][,ncol(stackmeta$metadata[[x]])], alpha = 0))

      coeffients <- lapply(1:length(stackmeta$metadata), function(x) coef(ridge.model[[x]], s = cv_ridge[[x]]$lambda.min)[-1])

      weights_ridge <- lapply(1:length(stackmeta$metadata), function(x) as.matrix(coeffients[[x]]))

      weights1_ridge <-  lapply(1:length(stackmeta$metadata), function(x) data.frame(h = x, model.weights = weights_ridge[[x]]))

      weightsDF_ridge <- dplyr::bind_rows(lapply(weights1_ridge, function(x) x %>%dplyr::mutate(model = stackmeta$models)))

      result <- structure(list(weights =  weightsDF_ridge, metalearner = "Ridge", comb.method = "stack"))

      class(result) <- "weight"

      return(result)

    }
                                                 
   # learning weights using lasso regression
                                                 
    else if (metalearner=="Lasso"){

      lasso.model <- lapply(1:length(stackmeta$metadata), function(x) glmnet::glmnet(stackmeta$metadata[[x]][,1:length(stackmeta$models)], stackmeta$metadata[[x]][,ncol(stackmeta$metadata[[x]])], alpha = 1))

      cv_lasso <- lapply(1:length(stackmeta$metadata), function(x) glmnet::cv.glmnet(stackmeta$metadata[[x]][,1:length(stackmeta$models)], stackmeta$metadata[[x]][,ncol(stackmeta$metadata[[x]])], alpha = 1))

      coeffients <- lapply(1:length(stackmeta$metadata), function(x) coef(lasso.model[[x]], s = cv_lasso[[x]]$lambda.min)[-1])

      weights_lasso <- lapply(1:length(stackmeta$metadata), function(x) as.matrix(coeffients[[x]]))

      weights1_lasso <-  lapply(1:length(stackmeta$metadata), function(x) data.frame(h = x, model.weights = weights_lasso[[x]]))

      weightsDF_lasso <- dplyr::bind_rows(lapply(weights1_lasso, function(x) x %>% dplyr::mutate(model = stackmeta$models)))

      result <- structure(list(weights =  weightsDF_lasso, metalearner = "Lasso", comb.method = "stack"))

      class(result) <- "weight"

      return(result)

    }
                                                 
# learning weights using elastic net regression
                                                 
    else if (metalearner=="Elastic"){

      elastic.model <- lapply(1:length(stackmeta$metadata), function(x) glmnet::glmnet(stackmeta$metadata[[x]][,1:length(stackmeta$models)], stackmeta$metadata[[x]][,ncol(stackmeta$metadata[[x]])], alpha = 0.5))

      cv_elnet <- lapply(1:length(stackmeta$metadata), function(x) glmnet::cv.glmnet(stackmeta$metadata[[x]][,1:length(stackmeta$models)], stackmeta$metadata[[x]][,ncol(stackmeta$metadata[[x]])], alpha = 0.5))

      coeffients <- lapply(1:length(stackmeta$metadata), function(x) coef(elastic.model[[x]], s = cv_elnet[[x]]$lambda.min)[-1])

      weights_elastic <- lapply(1:length(stackmeta$metadata), function(x) as.matrix(coeffients[[x]]))

      weights1_elastic <-  lapply(1:length(stackmeta$metadata), function(x) data.frame(h = x, model.weights = weights_elastic[[x]]))

      weightsDF_elastic <- dplyr::bind_rows(lapply(weights1_elastic, function(x) x %>% dplyr::mutate(model =stackmeta$models)))

      result <- structure(list(weights =  weightsDF_elastic, metalearner = "Elastic",  comb.method = "stack"))

      class(result) <- "weight"

      return(result)

    }
                                                   
# learning weights using non-negative least squared regression
                                                   
    else if(metalearner=="nnls") {

      coeffients <-lapply(1:length(stackmeta$metadata), function(x) nnls::nnls(stackmeta$metadata[[x]][,1:length(stackmeta$models)], stackmeta$metadata[[x]][,ncol(stackmeta$metadata[[x]])])$x)

      weights_nnls <- lapply(1:length(stackmeta$metadata), function(x) as.matrix(coeffients[[x]]))

      weights1_nnls <-  lapply(1:length(stackmeta$metadata), function(x) data.frame(h = x,  model.weights = weights_nnls[[x]]))

      weightsDF_nnls <- dplyr::bind_rows(lapply(weights1_nnls, function(x) x %>% dplyr::mutate(model = stackmeta$models)))

      result <- structure(list(weights = weightsDF_nnls, metalearner = "nnls", comb.method = "stack"))

      class(result) <- "weight"

      return(result)
    }
                                                
# learning weights using standard linear squared regression
                                                
    else if (metalearner=="Linear"){

      linear.model <- lapply(1:length(stackmeta$metadata), function(x) lm(rate~., data =  as.data.frame(stackmeta$metadata[[x]])))

      coeffients <- lapply(1:length(stackmeta$metadata), function(x) unname(coef(linear.model[[x]])[-1]))

      weights_linear <- lapply(1:length(stackmeta$metadata), function(x) as.matrix(coeffients[[x]]))

      weights1_linear <- lapply(1:length(stackmeta$metadata), function(x) data.frame(h = x, model.weights = weights_linear[[x]]))

      weightsDF_linear <- dplyr::bind_rows(lapply(weights1_linear, function(x) x %>% dplyr::mutate(model =stackmeta$models)))

      result <- structure(list(weights = weightsDF_linear,  metalearner = "Linear", comb.method = "stack"))

      class(result) <- "weight"

      return(result)

    }

  }

}

#'
#' @param alpha is the level of the test with the default value of alpha = 0.1.
#'
#' @param B is the bootstrap samples with default sample of B = 5000.
#'
#' @param l is the block length with the default value of l=3.
#'
#' @return Returns an object of class \code{CoMoMo} with the following components:
#'
#' \item{weights}{Returns the combination weights for different horizons.}
#' \item{comb.method}{Returns the combination method}
#' \item{method}{Returns the trainining method either cv or sv.}
#' \item{cvmse}{Returns the cross-validation mean squared error for each all the indiviual models for different horizons.}                                            
#' \item{Selected models}{Returns the superior models selected.}
#'
#' @examples
#'
#' # define models
#' Dont run 
#' LC <- lc()
#' APC <- apc()
#' modlist <- list("LC"= LC, "APC" = APC)
#' mcs_weight_val <- mcs(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 15, method = "sv")
#' mcs_weight_cv <- mcs(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 15,  method = "cv")
                                              
                    
#' @export

mc <- function(Loss, B, l){

  # Matrix of boostrap indixes

  if (ncol(Loss) == 1) stop("MCS: Need more data. Only one model entered.")
  if(any(is.na(Loss))) stop("NAs in Loss are not allowed")
  if(any(abs(Loss) == Inf)) stop("Inf in Loss are not allowed")

  LbarB <- boot::tsboot(Loss, colMeans, R = B, sim = "fixed", l = l)$t
  Lbar <- colMeans(Loss)
  zeta.Bi <- t(t(LbarB)-Lbar)
  save.res <- c()

  for(j in 1:(ncol(Loss) - 1))

    {

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

{ # function to compute the combination weights using model confidence set

  if ( holdout < round(length(years.fit)/3,0))
  {
    cat(paste("Warning: holdout is small \n"))
  }

  if (B < 1000)
  {
    cat(paste("Warning: B is small \n"))
  }
  
  if ( alpha < 0 ||  alpha > 1) {
    stop(" alpha must be in (0,1)")
  }


  if (method!="cv" && method!="sv") stop("This is undefined method")

  # single validation set 
  
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

    model.weights <- c()

    for (model in names(models)) {
      if (model %in% names(models)[selected]) {model.weights = c(model.weights, 1/length( names(models)[selected])) } else {model.weights = c(model.weights, 0)}

    }

    out <- (dplyr::bind_rows(lapply(rep(list(as.data.frame(model.weights)), h), function(x) x%>%dplyr::mutate(model = names(models))))%>%dplyr::mutate(h = rep(1:h, each = length(models))))[,c("h","model.weights","model")]

    res <- structure(list(weights = out, method = "sv", selected =  selected, comb.method = "MCS"))

    class(res) <- "weight"

    return(res)

  }
                                    
# cross-validation
                                    
  else if ( method=="cv")

  {
    output0 <- cvloss(models = models, data = data, ages.fit = ages.fit, years.fit = years.fit, h = h, Dxt = Dxt, Ext = Ext, ages = ages, years = years)

    MCS_cv <- lapply(1:h, function(x) mc((output0$cvloss[[x]])^2,  B = B, l = l)$test)

    selected <- lapply(MCS_cv, function(x)  which(x > alpha))

    model.weights <- c()

    for (i in 1:h)

    {
      for (model in  names(models)) {
        if (model %in%  names(models)[selected[[i]]]) {model.weights = c( model.weights, 1/length( names(models)[selected[[i]]])) } else { model.weights = c(model.weights, 0)

        }
      }
    }

    output <- data.frame(h = rep(1:h, each = length(models)),  model.weights, model = rep(unlist(names(models))), stringsAsFactors = F)

    result <- structure(list(weights = output, method = "cv", selected =  selected,  cvmse =  output0$CVE, comb.method = "MCS"))

    class(result) <- "weight"

    return(result)

  }

}
 
#' @export                       
                       
frequentist <- function(models, data = NULL, Dxt = NULL, Ext = NULL, ages.fit = NULL, years.fit = NULL,
                        ages = NULL, years = NULL, h = NULL)
  
{  # estimates the weights using the cvmse for each horizon 
   # weight <- cvme/sum(cvmse)
  
  Specified_Models <- lapply(1:length(models), function(x) names(models[x]))
  
  output0 <- cvloss(models = models, data = data, ages.fit = ages.fit, years.fit = years.fit, h = h, Dxt = Dxt, Ext = Ext, ages = ages, years = years)
  
  fweights<-  lapply(1:h, function(x) data.frame(h = x, weights = dplyr::bind_rows(output0$cvmse[[x]]/sum(dplyr::bind_rows(output0$cvmse[[x]])))))
  
  out <- dplyr::bind_rows(lapply(fweights, function(x) x %>%dplyr::mutate(model = unlist(Specified_Models))))
  
  output <- out%>%dplyr::rename(model.weights = cv.mse)
  
  result <- structure(list(weights = output, cvmse =  output0$CVE, comb.method = "frequentist"))
  
  class(result)<- "weight"

}
