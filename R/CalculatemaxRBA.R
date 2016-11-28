#' Calculate maxRBA
#'
#' \code{CalculatemaxRBA} calculates the maximum residual after baseline-adjustment for each viral strain
#' 
#' Calculates the baseline-adjusted fold change for each strain of virus
#' using (unnormalized) fold change and baseline titers. Linear regression or
#' an exponential curve is used to remove the effect of baseline titers on fold changes.
#' The score function (\code{scoreFun}) is used to combine the adjusted fold change across
#' multiple strains.
#' Missing (\code{NA}) values are handled by being returned as missing in the
#' endpoints in the output
#'
#' @param dat_list a named list like the one returned by \code{\link{FormatTiters}}.
#' @param subjectCol the name of the column specifying a subject ID. Default is "SubjectID".
#' @param method a character string specifying the method used to model the relationship between day 0 and fold change values. One of either "lm" for a linear model or "exp" for an exponential model.
#' @param yMinZero a logical specifying whether fitted y values below 0 should be set to 0.
#' @param scoreFun a function applied to all (potentially scaled) residuals for each subject to determine the endpoint. Default is \code{max} but \code{sum} may also be useful to quantify the total response.
#' @param normalize Logical specifying whether residuals should be normalized with the inverse normal transform. Default is \code{TRUE}.
#' @param discretize a vector of quantiles in (0, 0.5] specifying where to make the cutoff for low, moderate and high responses. Default is 20\% and 30\%.
#' @param scaleResiduals Logical. Should residuals be scaled inversely by the
#'                       square of the confidence intervals from the linear model.
#' @param responseLabels names for low, moderate and high responses 
#' @param na_action how should missing \code{NA} values be treated. Default is "na.fail"
#' @param ... Additional arguments passed to \code{lm} if method == "lm" or \code{nls} if method == "exp"
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{models}{the models calculated on each strain separately (with names the same as on \code{dat_list})}
#'   \item{residualMatrix}{the matrix of residuals}
#'   \item{maxRBA}{a list containing the continuous and discrete maxRBA metrics}
#' }
#' @seealso \code{lm, nls}
#' @author Stefan Avey
#' @export
#' @examples
#' ## Prepare the data
#' titer_list <- FormatTiters(Year2_Titers)
#'
#' ## Using a linear fit
#' endpoints <- CalculatemaxRBA(titer_list, method = "lm")
#' summary(endpoints)
#' ## Get discrete endpoints using upper/lower 30%
#' endpoints$maxRBA_d30
#'
#' ## Get endpoints with a 50% split into high and low
#' endpoints <- CalculatemaxRBA(titer_list, method = "exp", discretize = 0.5)
#' endpoints$maxRBA_d50
CalculatemaxRBA <- function(dat_list, subjectCol = "SubjectID",
                             method = c("exp", "lm"), yMinZero = FALSE,
                             scoreFun = max, discretize = c(0.2, 0.3),
                             normalize = TRUE, scaleResiduals = FALSE,
                             responseLabels = paste0(c("low", "moderate", "high"),
                                 "Responder"), na_action = "na.fail",
                             ...) {
  method <- match.arg(method)
  if(length(unique(lapply(dat_list, dim))) != 1) {
    stop("Each data frame in `dat_list` must have the same dimensions")
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
  for(i in seq_along(dat_list)) {
    dat <- dat_list[[i]]
    ## Check if arranged in order of Pre column
    ord <- order(dat$Pre)
    if(!all(ord == 1:nrow(dat))) {
      stop("`dat_list[[", i, "]]` is not ordered by 'Pre' column. Use only output from `FormatTiters`!")
    }
    if(method == "lm") {
      model <- lm(data = dat, formula = "FC ~ Pre", na.action = na_action, ...)
    } else if (method == "exp") {
        model <- nls(data = dat, formula = "FC ~ exp(a + b * Pre)",
                     start = list(a = 0, b = 0),
                     na.action = na_action, ...)
      }
    if(yMinZero && method == "lm") {
      residuals <- residuals(model)
      setToZero <- fitted(model) < 0
      residuals[setToZero] <- 0
    } else {
        residuals <- residuals(model)
      }
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
  colnames(residual_mat) <- names(dat_list)
  if(normalize) {
    residual_mat <- apply(residual_mat, 2, .INT)
  }
  maxRBA <- apply(residual_mat, 1, scoreFun, na.rm = TRUE)
  ## Calculated discretized metrics
  disList <- vector(mode = "list", length = length(discretize))
  names(disList) <- discretize
  for(dis in discretize) {
    tmp <- rep(NA, nrow(dat))
    names(tmp) <- names(maxRBA)
    lowR <- maxRBA <= quantile(maxRBA, dis, na.rm = TRUE)
    highR <- maxRBA >= quantile(maxRBA, 1 - dis, na.rm = TRUE)
    modR <- !(lowR | highR | is.na(maxRBA))
    tmp[lowR] <- 0
    tmp[modR] <- 1
    tmp[highR] <- 2
    disList[[as.character(dis)]] <- factor(tmp, labels = responseLabels, levels = 0:2)
  }
  disList <- setNames(disList, paste0("maxRBA_d",
                                      as.character(as.numeric(names(disList))*100)))
  names(model_list) <- names(dat_list)
  return(c(models = list(model_list), residualMatrix = list(residual_mat),
           maxRBA = list(maxRBA), disList))
}
