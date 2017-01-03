#' Calculate MFC
#'
#' \code{Calculate_MFC} calculates the (log-transformed) maximum fold change over all strains.
#'
#'
#' 
#' @param dat_list a named list like the one returned by \code{\link{FormatTiters}}.
#' @param subjectCol the name of the column specifying a subject ID. Default is "SubjectID".
#' @param discretize a vector of quantiles in (0, 0.5] specifying where to make the cutoff for low, moderate and high responses. Default is 20\% and 30\%.
#' @param responseLabels names for low, moderate and high responses 
#' @return A list with the following elements:
#' \describe{
#'   \item{MFC}{a named vector containing the continuous MFC endpoints}
#'   \item{MFC_d<X>}{a named vector containing the discrete MFC endpoint with a cutoff at <X>}
#'   \item{...}{Other named vectors containing discrete MFC endpoints}
#' }
#' @return A named list containing the MFC for each subject and any discretized metrics
#' @author Stefan Avey
#' @import dplyr
#' @importFrom stats quantile setNames
#' @export
#' @examples
#' ## Prepare the data
#' titer_list <- FormatTiters(Year2_Titers)
#'
#' Calculate_MFC(titer_list)
Calculate_MFC <- function(dat_list, subjectCol = "SubjectID", discretize = c(0.2, 0.3),
                         responseLabels = paste0(c("low", "moderate", "high"),
                             "Responder")) {
  if(length(unique(lapply(dat_list, dim))) != 1) {
    stop("Each data frame in `dat_list` must have the same dimensions")
  }
  result <- do.call(rbind.data.frame, dat_list) %>%
    group_by_(subjectCol) %>%
    mutate(MFC = max(FC)) %>%
    ungroup() %>%
    select_(subjectCol, "MFC") %>%
    unique()
  MFC <- result$MFC
  names(MFC) <- result[[subjectCol]]
  ## Calculated discretized metrics
  if(!is.null(discretize)) {
    disList <- vector(mode = "list", length = length(discretize))
    names(disList) <- discretize
    for(dis in discretize) {
      tmp <- rep(NA, nrow(dat_list[[1]]))
      names(tmp) <- names(MFC)
      lowR <- MFC <= quantile(MFC, dis, na.rm = TRUE)
      highR <- MFC >= quantile(MFC, 1 - dis, na.rm = TRUE)
      modR <- !(lowR | highR | is.na(MFC))
      tmp[lowR] <- 0
      tmp[modR] <- 1
      tmp[highR] <- 2
      disList[[as.character(dis)]] <- factor(tmp, labels = responseLabels, levels = 0:2)
    }
    disList <- setNames(disList, paste0("MFC_d",
                                        as.character(as.numeric(names(disList))*100)))
    out <- c(MFC = list(MFC), disList)
  } else {
      out <- list(MFC = MFC)
    }
  return(out)
}
