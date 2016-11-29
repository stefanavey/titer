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
#' @return A named vector containing the response ("NR", "X", or "R") for each subject
#' @author Stefan Avey
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
  for(i in seq_along(dat_list)) {
    dat <- dat_list[[i]]
    fourFC[, i] <- dat$FC >= 2
  }
  counts <- apply(fourFC, 1, sum)
  R <- counts >= 2
  NR <- counts == 0
  response <- rep("X", length(counts))
  missing <- apply(fourFC, 1, function(row) any(is.na(row)))
  response[missing] <- NA
  response[which(R)] <- "R"
  response[which(NR)] <- "NR"
  names(response) <- dat[[subjectCol]]
  return(response)
}
