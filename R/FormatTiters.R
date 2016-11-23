#' Format antibody titers.
#'
#' \code{FormatTiters} formats titers into a list with one tidy data frame per viral strain
#'
#'
#'
#' @param titers a data frame containing one row per subject per strain. The following columns are required:
#' \describe{
#'   \item{SubjectID}{Subject IDs (column name can vary)}
#'   \item{Strain}{The name of the viral strain for the observation}
#'   \item{Pre}{The pre-vaccination (or pre-infection) titer}
#'   \item{Post}{The post-vaccination (or post-infection) titer}
#'   \item{...}{Other columns which will be preserved}
#' }
#' @param log2Transform logical specifying whether titer values should be log2 transformed
#' @param fcMinZero should negative fold changes be set to 0? Default is \code{TRUE}
#'
#' @return a list of data frames with one data frame per viral strain containing the "Pre" and "Post" titer measurements (row names are removed).
#' 
#' @author Stefan Avey
#' @import dplyr
#' @export
#' @examples
#' titer_list <- FormatTiters(Year1_Titers, log2Transform = TRUE, fcMinZero = TRUE)
FormatTiters <- function(titers, log2Transform = TRUE, fcMinZero = TRUE)
{
  if(log2Transform) {
    message("- Log transforming Pre and Post columns")
    trans <- log2
  } else {
      trans <- identity
    }
  if(fcMinZero) {
    message("- Setting any negative log fold changes to 0")
  }
  titer_list <- list()
  strains <- sort(unique(titers$Strain))
  for(i in seq_along(strains)) {
    result <- titers %>%
      dplyr::filter(Strain == strains[i]) %>%
      dplyr::mutate(Post = trans(Post), Pre = trans(Pre)) %>%
      dplyr::mutate(FC = Post - Pre) %>%
      dplyr::arrange(Pre, FC) %>%
      distinct()
    if(fcMinZero) {
      result <- result %>% dplyr::mutate(FC = pmax(FC, 0))
    }
    rownames(result) <- NULL    # remove rownames which may cause problems later
    titer_list[[strains[i]]] <- result
  }
  return(titer_list)
}
