
#' @usage combines the fitted mortality models and combination weights.
#'
#' @return Returns an object of class \code{CoMoMo} with the following components:
#'
#' \item{Comb}{Returns the combinated forecasts for different horizon.}
#'
#'  @examples
#'
#'  modelFits <- fitCoMoMo(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit)
#'
#'  modcom <- CoMoMo(modelFits, weight = stack_nnls_weight)
#'
#'
#' @export

CoMoMo  <- function(object, weight = NULL)

{

  if(is.null(weight))

  {
    output <- structure(list(model.fits = object))

    class(output) <- "CoMoMo.simple"

    return(output)

  }

  else {

  output <- structure(list(model.fits = object, weight = weight))

  class(output) <- "CoMoMo"

  return(output)

  }

}



















