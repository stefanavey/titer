#################################################################################
## Calculate_StdNorm                                                            ##
#################################################################################
#' Calculate Normalized Titers
#'
#' \code{Calculate_StdNorm} calculates the standardized d0 or fc titers
#'
#' This must be run on only 1 cohort at a time because titers will be normalized
#' across all subjects. The median is used but unlike the original reference,
#' the standard deviation is calculated rather than the maximum absolute deviation.
#' 
#' @param dat       Data frame containing \code{fcStdCols}
#' @param type      What should be standarized. Either "d0", or "fc".
#' @param fcToOne   Logical. Are titer fold changes allowed to be less than 1
#'                  or should these be changed to 1 before standardization?
#'                  Default is FALSE and no changes will be made. Only relevant
#'                  when \code{type == "fc"}
#' @param idCol     Name of column containing subject IDs
#' @param cols      column names containing the titer measurements
#'                  for each strain
#' @return          A data frame like \code{dat} but with standarized columns added
#'
#' @import dplyr
#' @importFrom stats median sd setNames
#' 
#' @author Stefan Avey
#' @references Tsang JS, et al. (2014) Global analyses of human immune variation reveal baseline predictors of postvaccination responses. Cell 157(2):499-513.
#' @export
#' @examples
#' ## First Example
#' 
Calculate_StdNorm <- function(dat, type, fcToOne = FALSE, idCol = "SubjectID",
                             cols = grep(paste0(type, "_[AB]"),
                                 colnames(dat), value = TRUE)) {
  ## Functions for calculating the std values
  .std <- function(x) { (x - median(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) }
  ## .mad <- function(x) {median(abs(x - median(x, na.rm = TRUE)), na.rm =TRUE)}
  ## .std <- function(x) { (x - median(x, na.rm = TRUE)) / .mad(x) }
  newCols <- gsub(pattern = type, replacement = paste0(type, "_std_norm"), cols)
  cols <- setNames(cols, newCols)
  ## Set fold change values less than 1 to 1
  if(fcToOne && type == "fc") {
    dat[,cols] <- apply(dat[,cols], 2, function(col) ifelse(col < 1, 1, col))
  }
  ## Remove any duplicated subjects and calculate std values
  res <- dat %>%
    select(which(colSums(!is.na(.)) > 0)) %>% # remove columns with all NA
    select(one_of(idCol, cols)) %>%        
    distinct() %>%                      # remove duplicates
    mutate_each_(funs = funs(.std), vars = cols) %>%
    select(one_of(idCol, newCols))
  ret <- full_join(dat, res, by = idCol)
  return(ret)
}

#################################################################################
## Calculate_D0NormPaired                                                       ##
#################################################################################
#' Calculate_D0NormPaired
#'
#' \code{Calculate_D0NormPaired} calculates the normalized day 0 titer paired with the titer with maximum normalized fold change 
#'
#' If there are multiple strains that have the maximal fold change,
#' choose the day 0 titer that is higher since this will allow
#' for a greater adjustment and better chance of being a high responder.
#'
#' Column names containing the day 0 titers for each strain
#' standardized across subjects are assumed to follow the same
#' pattern as \code{fcStdCols} with "d0" replacing "fc" in the name.
#'
#' @param dat       data frame containing \code{fcStdCols}
#' @param fcStdCols column names containing the titer fold changes for each strain
#'                  standardized across subjects
#' @return          a numeric vector containing the values from d0StdCols that
#'                  correspond to the maximum over the strains of fcStdCols
#' @author Stefan Avey
#' @export
#' @examples
#' ## First Example
#' 
Calculate_D0NormPaired <- function(dat,
                                  fcStdCols = grep("fc_std_norm",
                                      colnames(dat), value = TRUE)) {
  maxStrains <- apply(dat[,fcStdCols], 1, max)
  d0Paired <- rep(NA, nrow(dat))
  for(i in 1:nrow(dat)) {
    if(!is.na(maxStrains[i])) {
      maxCols <- fcStdCols[which(dat[i, fcStdCols] == maxStrains[i])]
      d0Cols <- gsub("fc", "d0", maxCols)
      d0Paired[i] <- max(dat[i,d0Cols])
    } else {
        d0Paired[i] <- NA
      }
  }
  return(data.frame(dat, d0_norm_paired = d0Paired))
}


#################################################################################
## CalculatePadjMFC                                                            ##
#################################################################################
#' Calculate_padjMFC
#'
#' \code{Calculate_padjMFC} calculates the paired, adjusted maximum fold change (padjMFC)
#' 
#' Calculate the paired, adjusted maximum fold change (padjMFC) from
#' fc_norm_max_ivt and d0_norm_paired using linear regression to
#' remove the effect of baseline titers. Missing (\code{NA}) values are handled
#' and any missing values in fcCol and d0Col will also be missing in the output.
#'
#' @param dat the data containing the columns \code{fcCol} and \code{d0Col}
#' @param fcCol character string specifying the name of the fold change column from \code{dat}
#' @param d0Col character string specifying the name of the day 0 column from \code{dat}
#' @param discretize a vector of quantiles in (0, 0.5] specifying where to make the cutoff for low, moderate and high responses. Default is 20\% and 30\%.
#' @param scaleResiduals Logical. Should residuals be scaled inversely by the
#'                       square of the confidence intervals from the linear model.
#' @param responseLabels names for low, moderate and high responses 
#' @param ... Additional arguments passed to \code{lm}
#'
#' @return A list with the first element named "linearModel" for the linear model and then "padjMFC" containing the continuous padjMFC metric and one additional element for each value of discretize giving the discrete labels.
#'
#' @seealso \code{lm}
#' 
#' @author Stefan Avey
#' @importFrom stats as.formula lm quantile setNames predict
#' @export
#' @examples
#' ## First Example
#'
Calculate_padjMFC <- function(dat, fcCol = "fc_norm_max_ivt", d0Col = "d0_norm_paired",
                             discretize = c(0.2, 0.3), scaleResiduals = FALSE,
                             responseLabels = paste0(c("low", "moderate", "high"),
                                 "Responder"), ...) {
  if(! fcCol %in% colnames(dat)) {
    stop("fcCol ('", fcCol, "') not found in column names of dat")
  }
  if(! d0Col %in% colnames(dat)) {
    stop("d0Col ('", d0Col, "') not found in column names of dat")
  }
  ## Calculate padjMFC
  form <- as.formula(paste(fcCol, "~", d0Col))
  model <- lm(data = dat, formula = form, na.action = "na.omit", ...)
  padjMFC <- rep(NA, nrow(dat))
  if(scaleResiduals) {
    cis <- predict(model, na.action = "na.omit",
                   newdata = dat, se.fit = FALSE,
                   level = 0.95, interval = "confidence")
    intervals <- apply(cis, 1, function(row) { return(row["upr"] - row["lwr"]) })
    padjMFC[as.numeric(names(model$residuals))] <- model$residuals / intervals^2    
  } else {
      padjMFC[as.numeric(names(model$residuals))] <- model$residuals
    }
  ## Calculated discretized metrics
  disList <- vector(mode = "list", length = length(discretize))
  names(disList) <- discretize
  for(dis in discretize) {
    tmp <- rep(NA, nrow(dat))
    lowR <- padjMFC <= quantile(padjMFC, dis, na.rm = TRUE)
    highR <- padjMFC >= quantile(padjMFC, 1 - dis, na.rm = TRUE)
    modR <- !(lowR | highR | is.na(padjMFC))
    tmp[lowR] <- 0
    tmp[modR] <- 1
    tmp[highR] <- 2
    disList[[as.character(dis)]] <- factor(tmp, labels = responseLabels, levels = 0:2)
  }
  disList <- setNames(disList, paste0("padjMFC_d",
                                      as.character(as.numeric(names(disList))*100)))
  return(c(linearModel = list(model), padjMFC = list(padjMFC), disList))
}
