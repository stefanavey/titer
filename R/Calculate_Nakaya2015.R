#' Calculate Nakaya2015
#'
#' \code{Calculate_Nakaya2015} calculates the endpoint used in Nakaya et al. 2015
#' 
#' First calculate the maximum fold change (MFC) derived titer metric described
#' in Nakaya et al. 2015. Then check whether both of these conditions are satisfied:
#'     i) MFC is at least a 4-fold increase
#'    ii) The "Post" antibody titer is 1:40 or more for at least 1 strain
#' Subjects are classified as high responders if they satisfy both conditions and
#' low responders otherwise.
#'
#' Missing (\code{NA}) values are handled by being returned as missing in the
#' endpoints in the output
#'
#' @param dat_list a named list like the one returned by \code{\link{FormatTiters}}.
#' @param subjectCol the name of the column specifying a subject ID. Default is "SubjectID".
#' @param responseLabels names for low and high responses 
#' @param na_action how should missing \code{NA} values be treated. Default is "na.fail"
#' @param ... Additional arguments passed to \code{lm}
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{data}{a data frame containing the MFC and indicator variables that determine whether subject is a low or high responder (see details)}
#'   \item{Nakaya2015}{a named vector containing the discretized endpoint}
#' }
#' @seealso \code{CalculateMFC}
#' @references Nakaya HI, et al. (2015) Systems Analysis of Immunity to Influenza Vaccination across Multiple Years and in Diverse Populations Reveals Shared Molecular Signatures. Immunity 43(6):1186-1198.
#' @author Stefan Avey
#' @export
#' @examples
#' ## Prepare the data
#' titer_list <- FormatTiters(Year2_Titers)
#'
#' ## Calculate the endpoint
#' endpoints <- Calculate_Nakaya2015(titer_list)
#' summary(endpoints)
Calculate_Nakaya2015 <- function(dat_list, subjectCol = "SubjectID",
                                responseLabels = paste0(c("low", "high"),
                                    "Responder"), na_action = "na.fail", ...) {
  if(length(unique(lapply(dat_list, dim))) != 1) {
    stop("Each data frame in `dat_list` must have the same dimensions")
  }
  mfc <- CalculateMFC(dat_list, subjectCol = subjectCol, discretize = NULL)
  mfc_dat <- data.frame(mfc, subject = names(mfc$MFC), stringsAsFactors = FALSE)
  cond1 <- do.call(rbind.data.frame, dat_list) %>%
    mutate_at(subjectCol, as.character) %>%
    group_by_(subjectCol) %>%
    summarize(PostOver40 = any(Post >= log2(40)))
  cond2 <- data.frame(subject = names(mfc$MFC), FC4 = mfc$MFC >= log2(4),
                      stringsAsFactors = FALSE)
  result <- full_join(cond2, cond1, by = c(subject = subjectCol)) %>%
    left_join(mfc_dat, by = "subject") %>%
    mutate(Response = ifelse(FC4 & PostOver40, responseLabels[2], responseLabels[1]))
  ## Rename columns to use user-supplied subjectID
  colnames(result) <- gsub("^subject$", subjectCol, colnames(result))
  nakaya2015 <- result$Response
  names(nakaya2015) <- result[[subjectCol]]
  return(list(data = result, Nakaya2015 = nakaya2015))
}
