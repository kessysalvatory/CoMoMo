
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
#'  @examples
#'
#'  comb <- CoMoMo(modlist, data = DataStMoMo, ages.fit = ages.fit, years.fit = years.fit, h = 5)
#'
#'  weights <- bma(modlist, data = DataStMoMo, ages.fit = ages.fit, years.fit = years.fit, h = 5, method = "cv")
#'
#'  comb <- CoMoMo(modlist, weight =  weights, data = DataStMoMo,  ages.fit = ages.fit, years.fit = years.fit, h = 5)
#'
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



















