#' @title Assigning names to different specified models.
#'
#' @models is the list of the mortality models
#'
#' @modelNames is the list of the names of the models
#'
#' @return A list of named mortality models.

#' @export

model <- function (models, modelNames) {

  for (i in 1:length(modelNames))

  {
    if (class(modelNames[[i]])!="character") stop("Names should be character")
  }

  if(length(models)!=length(modelNames)) stop("Not of equal length")
  if(length(models) < 2) stop("Specify more than one model.")

  output <- setNames(models, modelNames)

  return(output)

}




