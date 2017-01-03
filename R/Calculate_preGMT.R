#' Calculate pre-GMT
#'
#' \code{Calculate_preGMT} calculates the log-transformed pre-vaccination geometric mean titer (pre-GMT)
#'
#' Non-logged HAI titers for each strain are used to calculate the geometric mean and
#' the geometric mean for each subject is subsequently log2-transformed.
#' 
#' @param dat_list a named list like the one returned by \code{\link{FormatTiters}}.
#' @param subjectCol the name of the column specifying a subject ID. Default is "SubjectID".
#' @return A named vector containing the pre-GMT for each subject
#' @author Stefan Avey
#' @import dplyr
#' @export
#' @examples
#' ## Prepare the data
#' titer_list <- FormatTiters(Year2_Titers)
#'
#' Calculate_preGMT(titer_list)
Calculate_preGMT <- function(dat_list, subjectCol = "SubjectID") {
  .GeometricMean <- function(x, na.rm = TRUE){
    if(any(x <= 0)) {
      warning("Some values of x are not positive, only positive values will be used")
    }
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  result <- do.call(rbind.data.frame, dat_list) %>%
    group_by_(subjectCol) %>%
    mutate(preGMT = log2(.GeometricMean(2^Pre))) %>%
    ungroup() %>%
    select_(subjectCol, "preGMT") %>%
    unique()
  response <- result$preGMT
  names(response) <- result[[subjectCol]]
  return(response)
}
