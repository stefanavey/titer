#' Calculate Thakar2015
#'
#' \code{Calculate_Thakar2015} calculates a response definition based on a Four Fold change to At Least Two Strains (Thakar2015).
#' 
#' Subjects are responders (default "R") if they acheive a 4-fold or greater fold change in
#' titer to at least 2 strains, nonresponders (default "NR") if they do not acheive a 4-fold
#' or greater fold change in titer to any strain, and intermediate (default "X") otherwise.
#' Missing (\code{NA}) values are handled by being returned as missing in the
#' endpoints in the output
#'
#' @param dat_list a named list like the one returned by \code{\link{FormatTiters}}.
#' @param subjectCol the name of the column specifying a subject ID. Default is "SubjectID".
#' @param responseLabels names for low, middle, and high responses 
#' @return A named list with 1 element named "Thakar2015" containing the response ("NR", "X", or "R").
#' @references Thakar J, et al. (2015) Aging-dependent alterations in gene expression and a mitochondrial signature of responsiveness to human influenza vaccination. Aging (Albany NY) 7(1):38-52.
#' @author Stefan Avey
#' @import dplyr
#' @export
#' @examples
#' ## Prepare the data
#' titer_list <- FormatTiters(Year2_Titers)
#'
#' Calculate_Thakar2015(titer_list)
Calculate_Thakar2015 <- function(dat_list, subjectCol = "SubjectID",
                            responseLabels = c("NR", "X", "R")) {
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
    mutate(ffalts = NA) %>%
    mutate(ffalts = ifelse(numFourFC >= 2, "R", ffalts)) %>%
    mutate(ffalts = ifelse(numFourFC == 1, "X", ffalts)) %>%    
    mutate(ffalts = ifelse(numFourFC == 0, "NR", ffalts))
  ffalts <- result %>%
    select(ffalts) %>%
    unlist()
  names(ffalts) <- result[[subjectCol]]
  return(list(Thakar2015 = ffalts))
}
