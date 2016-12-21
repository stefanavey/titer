#' Calculate whoResp
#'
#' \code{CalculatewhoResp} calculates a response definition similar to the WHO defintion using a 4-fold cutoff.
#' 
#' Subjects are responders ("R") if they acheive a 4-fold or greater fold change in
#' titer to at least 2 strains, nonresponders ("NR") if they do not acheive a 4-fold
#' or greater fold change in titer to any strain, and intermediate ("X") otherwise.
#' Missing (\code{NA}) values are handled by being returned as missing in the
#' endpoints in the output
#'
#' @param dat_list a named list like the one returned by \code{\link{FormatTiters}}.
#' @param subjectCol the name of the column specifying a subject ID. Default is "SubjectID".
#' @return A named list with 1 element named "whoResp" containing the response ("NR", "X", or "R").
#' @author Stefan Avey
#' @import dplyr
#' @export
#' @examples
#' ## Prepare the data
#' titer_list <- FormatTiters(Year2_Titers)
#'
#' CalculatewhoResp(titer_list)
CalculatewhoResp <- function(dat_list, subjectCol = "SubjectID") {
  fourFC <- data.frame(matrix(nrow = nrow(dat_list[[1]]),
                              ncol = length(dat_list),
                              dimnames = list(NULL, names(dat_list))),
                       check.names = FALSE)
  if(length(unique(lapply(dat_list, dim))) != 1) {
    stop("Each data frame in `dat_list` must have the same dimensions")
  }
  result <- do.call(rbind.data.frame, dat_list) %>%
    mutate(fourFC = FC >= 2) %>%
    group_by_(subjectCol) %>%
    summarize(numFourFC = sum(fourFC)) %>%
    ungroup() %>%
    mutate(whoResp = NA) %>%
    mutate(whoResp = ifelse(numFourFC >= 2, "R", whoResp)) %>%
    mutate(whoResp = ifelse(numFourFC == 1, "X", whoResp)) %>%    
    mutate(whoResp = ifelse(numFourFC == 0, "NR", whoResp))
  whoResp <- result %>%
    select(whoResp) %>%
    unlist()
  names(whoResp) <- result[[subjectCol]]
  return(list(whoResp = whoResp))
}
