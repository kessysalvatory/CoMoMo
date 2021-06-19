
#' @usage combines the fitted mortality models and combination weights.
#'
#' @return Returns an object of class \code{CoMoMo} with the following components:
#'
#' \item{Comb}{Returns the fitted GAPC mortality models.}
#' \item{Comb}{Returns the combination weights.}
#' \item{Comb}{Returns the combination method.}
#' \item{Comb}{Returns the meta-learner if the combination method is the stacked regression ensemble.}
#'
#'  @examples
#'  Dont run 
#'  modelFits <- fitCoMoMo(models, data = DataStMoMo, ages.fit = agesFit, years.fit = yearsFit)
#'  modcom <- CoMoMo(modelFits, weight = stack_nnls_weight) 
#'
#' @export

CoMoMo  <- function(object, weight = NULL)

{ # when model combination is simple model averaging. 

  if(is.null(weight))

  {
    output <- structure(list(model.fits = object))

    class(output) <- "CoMoMo.simple"

    return(output)

  }
   # Other model combinations 

  else {

  output <- structure(list(model.fits = object, weight = weight))

  class(output) <- "CoMoMo"

  return(output)

  }

}



















