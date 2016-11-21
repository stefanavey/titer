##' @title FormatTiters
##'
##' @description
##' \code{FormatTiters} formats titers into a list with one data frame per viral strain
##'
##' @param titers a data frame containing the titer information
##' @param strains the names of the virus strains
##' @param subjectCol the name of the column specifying a subject ID. Default is "SubjectID". 
##' @param otherCols a character vector specifying which additional columns of titers to retain. (Defaults to an empty character vector).
##' @param d0Cols the column names of day 0 (baseline) columns 
##' @param fcCols the column names of fold change columns 
##' @param fcMinZero should negative fold changes be set to 0? Default is \code{TRUE}
##' @param log2Transform logical specifying whether titer values should be log2 transformed
##' @details
##'
##' @return a list of data frames with one data frame per viral strain containing the baseline ("d0"), fold change ("fc") and any other columns specified by the \code{otherColumns} argument.
##' 
##' 
##' @author Stefan Avey
##' @keywords HIPC
##' @export
##' @examples
##' \dontrun{
##' ## Example using the master phenotype file
##' library(dplyr)
##' titers <- master %>%
##'   filter(Year == 1, AgeGroup %in% "Young", !is.na(whoResp))
##' titer_list <- FormatTiters(titers,
##'                            strains = c("A_California_7_2009",
##'                                "A_Perth_16_2009",
##'                                "B_Brisbane_60_2008"))
##' }
FormatTiters <- function(titers, strains,
                         subjectCol = "SubjectID",
                         otherCols = vector(mode = "character"),
                         d0Cols = paste0("d0_", strains),
                         fcCols = paste0("fc_", strains),
                         fcMinZero = TRUE, log2Transform = TRUE)
{
  if(log2Transform) {
    trans <- log2
  } else {
      trans <- identity
    }
  titer_list <- list()
  for(i in seq_along(strains)) {
    titers$d0 <- trans(titers[[d0Cols[i]]])
    titers$fc <- trans(titers[[fcCols[i]]])
    if(fcMinZero) {
      titers$fc <- pmax(0, titers$fc)
    }
    selCols <- c(subjectCol, otherCols, "d0", "fc")
    result <- titers[,selCols]
    result <- result[!duplicated(result),]
    ord <- order(result$d0, result$fc)
    result <- result[ord,]
    rownames(result) <- NULL
    titer_list[[strains[i]]] <- result
  }
  return(titer_list)
}

