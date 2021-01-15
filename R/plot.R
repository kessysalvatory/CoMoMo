#'
#' @title Plot the weights from different model combination methods
#'
#' @param object an object containing the weights for different horizons.
#'
#'  @examples
#'
#'  weights <- bma(modlist, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 5, method = "cv")
#'  plot(weights)

#' @export

plot.CoMoMo <- function(object, ...) {

  wplot <- ggplot2::ggplot(object$weights) + ggplot2::geom_line(aes(x = h, y = weights, group = model, colour = model, linetype = model), size = 2) +

    ggplot2::theme(plot.title = element_text(hjust = 0.5)) + ggplot2::labs(x = "Forecasting Horizon", y = "Weights", title = "Model Weights by Forecasting Horizon")

  return(wplot)

}


#' @export

plot.CoMoMo.CoMoMo <- function(object, ...) {


  wplot  <- ggplot2::ggplot( object$mse) + ggplot2::geom_line(aes(x = h , y = mse, group = model, colour = model, linetype = model), size = 2) + ggplot2::theme(plot.title = element_text(hjust = 0.5)) + ggplot2::labs(x = "Forecasting Horizon", y = "MSE", title = "Mean Square Error by Forecasting Horizon")

  return(wplot)

}








