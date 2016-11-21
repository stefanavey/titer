#' Calculate SAdjMFC
#'
#' \code{CalculateSAdjMFC} calculates the baseline-adjusted maximum fold change (MFC)
#' for each viral strain
#' 
#' Calculates the baseline-adjusted fold change for each strain of virus
#' using (unnormalized) fold change and baseline titers. Linear regression or
#' an exponential curve is used to remove the effect of baseline titers on fold changes.
#' The score function (\code{scoreFun}) is used to combine the adjusted fold change across
#' multiple strains.
#' Missing (\code{NA}) values are handled by being returned as missing in the
#' endpoints in the output
#'
#' @param datList a list with one data frame for each strain and each data frame containing the columns \code{fcCol} and \code{d0Col}. The order of each data frame must be the same and they must be the same dimensions. In addition, each data frame must be sorted by \code{d0Col} from low to high.
#' @param subjectCol the name of the column specifying a subject ID. Default is "SubjectID".
#' @param method a character string specifying the method used to model the relationship between day 0 and fold change values. One of either "lm" for a linear model or "exp" for an exponential model.
#' @param scoreFun a function applied to all (potentially scaled) residuals for each subject to determine the endpoint. Default is \code{max} but \code{sum} may also be useful to quantify the total response.
#' @param fcCol character string specifying the name of the fold change column in each element of \code{datList}
#' @param d0Col character string specifying the name of the day 0 column in each element of \code{datList}
#' @param normalize Logical specifying whether residuals should be normalized with the inverse normal transform. Default is \code{TRUE}.
#' @param discretize a vector of quantiles in (0, 0.5] specifying where to make the cutoff for low, moderate and high responses. Default is 20\% and 30\%.
#' @param scaleResiduals Logical. Should residuals be scaled inversely by the
#'                       square of the confidence intervals from the linear model.
#' @param responseLabels names for low, moderate and high responses 
#' @param na_action how should missing \code{NA} values be treated. Default is "na.fail"
#' @param ... Additional arguments passed to \code{lm} if method == "lm" or \code{nls} if method == "exp"
#'
#' @return A list with the following elements:
#'         "models":         the models calculated on each strain separately (with names the same as on \code{datList})
#'         "residualMatrix": the matrix of residuals
#'         "SAdjMFC":        a list containing the continuous and discrete SAdjMFC metrics
#' @seealso \code{lm, nls}
#' @author Stefan Avey
#' @keywords HIPC
#' @export
#' @examples
#' ## First Example
#'
CalculateSAdjMFC <- function(datList, subjectCol = "SubjectID",
                             method = c("lm", "exp"),
                             scoreFun = max,
                             fcCol = "fc", d0Col = "d0", normalize = TRUE,
                             discretize = c(0.2, 0.3), scaleResiduals = FALSE,
                             responseLabels = paste0(c("low", "moderate", "high"),
                                 "Responder"), na_action = "na.fail",
                             ...) {
  method <- match.arg(method)
  if(length(unique(lapply(datList, dim))) != 1) {
    stop("Each data frame in `datList` must have the same dimensions")
  }
  if(! fcCol %in% Reduce(intersect, lapply(datList, colnames))) {
    stop("fcCol ('", fcCol, "') not found in column names of dat")
  }
  if(! d0Col %in% Reduce(intersect, lapply(datList, colnames))) {
    stop("d0Col ('", d0Col, "') not found in column names of dat")
  }
  if(method == "exp" && scaleResiduals) {
    warning("Scaling of residuals is not implemented for method == 'exp'.")
    scaleResiduals <- FALSE
  }
  ## Inverse Normal Transform
  .INT <- function(x, na.last = "keep",
                   ties.method = c("average", "first", "last",
                       "random", "max", "min"), ...) {
    ties.method <- match.arg(ties.method)
    xranks <- rank(x, na.last = na.last, ties.method = ties.method)
    tempp <- (xranks - 0.5)/length(xranks)
    return(qnorm(tempp, ...))
  }
  residuals_list <- model_list <- list()
  ## Calculate residual for each strain
  for(i in seq_along(datList)) {
    dat <- datList[[i]]
    rownames(dat) <- NULL               # most rownames will break this code
    ## Check if arranged in order of d0 column
    ord <- order(dat[[d0Col]])
    if(!all(ord == 1:nrow(dat))) {
      stop("`datList[[", i, "]]` is not ordered by ", d0Col, " from low to high!")
    }
    if(method == "lm") {
      form <- as.formula(paste(fcCol, "~", d0Col))
      model <- lm(data = dat, formula = form, na.action = na_action, ...)
    } else if (method == "exp") {
        form <- as.formula(paste(fcCol, "~ exp(a + b *", d0Col, ")"))
        model <- nls(data = dat, formula = form, start = list(a = 0, b = 0),
                     na.action = na_action, ...)
        ## dups <- duplicated(dat[[d0Col]])
        ## expFit <- data.frame(d0 = dat[[d0Col]][!dups],
        ##                      fc = predict(model)[!dups])
        ## residuals <- apply(dat, 1, function(row) {
        ##                           as.numeric(row[fcCol]) -
        ##                             as.numeric(expFit[expFit[["d0"]] == as.numeric(row[d0Col]),"fc"])
                                ## })
      }
    residuals <- residuals(model)
    names(residuals) <- dat[[subjectCol]]
    if(scaleResiduals) {
      cis <- stats::predict(model, na.action = na_action,
                            newdata = dat, se.fit = FALSE,
                            level = 0.95, interval = "confidence")
      intervals <- apply(cis, 1, function(row) { return(row["upr"] - row["lwr"]) })
      residuals <- residuals / intervals^2    
    }
    residuals_list[[i]] <- residuals[order(names(residuals))]
    model_list[[i]] <- model
  }
  if(length(residuals_list) > 1) {
    residual_mat <- Reduce(cbind, residuals_list)
  } else {
      residual_mat <- matrix(residuals_list[[1]], ncol = 1)
    }
  colnames(residual_mat) <- names(datList)
  if(normalize) {
    residual_mat <- apply(residual_mat, 2, .INT)
  }
  SAdjMFC <- apply(residual_mat, 1, scoreFun, na.rm = TRUE)
  ## Calculated discretized metrics
  disList <- vector(mode = "list", length = length(discretize))
  names(disList) <- discretize
  for(dis in discretize) {
    tmp <- rep(NA, nrow(dat))
    names(tmp) <- names(SAdjMFC)
    lowR <- SAdjMFC <= quantile(SAdjMFC, dis, na.rm = TRUE)
    highR <- SAdjMFC >= quantile(SAdjMFC, 1 - dis, na.rm = TRUE)
    modR <- !(lowR | highR | is.na(SAdjMFC))
    tmp[lowR] <- 0
    tmp[modR] <- 1
    tmp[highR] <- 2
    disList[[as.character(dis)]] <- factor(tmp, labels = responseLabels, levels = 0:2)
  }
  disList <- setNames(disList, paste0("SAdjMFC_d",
                                      as.character(as.numeric(names(disList))*100)))
  names(model_list) <- names(datList)
  return(c(models = list(model_list), residualMatrix = list(residual_mat),
           SAdjMFC = list(SAdjMFC), disList))
}
