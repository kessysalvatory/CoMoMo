#'
#' @title Plot the weights from different model combination methods
#'
#' @param object an object containing the weights for different horizons.
#'
#'  @examples
#'
#' weights <- bma(modlist, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit, h = 5, method = "cv")
#'  plot(weights)

#' @export

plot.weight <- function(object, size = 2, ...) {

  wplot <- ggplot2::ggplot(object$weights) + ggplot2::geom_line(aes(x = h, y = model.weights, group = model, colour = model, linetype = model), size = size) +

    ggplot2::theme(plot.title = element_text(hjust = 0.5)) + ggplot2::labs(x = "Forecasting Horizon", y = "Weights", title = "Model Weights by Forecasting Horizon")

  return(wplot)

}








