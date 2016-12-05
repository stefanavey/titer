#' Calculate TRI
#'
#' \code{CalculateTRI} calculates the Titer Response Index (TRI)
#' 
#' Calculates the Titer Response Index (TRI) defined in Bucasas et al. 2011
#' Missing (\code{NA}) values are handled by being returned as missing in the
#' endpoints in the output
#'
#' @param dat_list a named list like the one returned by \code{\link{FormatTiters}}.
#' @param subjectCol the name of the column specifying a subject ID. Default is "SubjectID".
#' @param discretize a vector of quantiles in (0, 0.5] specifying where to make the cutoff for low, moderate and high responses. Default is 20\% and 30\%.
#' @param responseLabels names for low, moderate and high responses 
#' @param na_action how should missing \code{NA} values be treated. Default is "na.fail"
#' @param ... Additional arguments passed to \code{lm}
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{models}{the models calculated on each strain separately (with names the same as on \code{dat_list})}
#'   \item{residualMatrix}{the matrix of residuals}
#'   \item{scores}{a data frame containing the four scores (before scaling)}
#'   \item{TRI}{a named vector containing the continuous TRI endpoint}
#'   \item{TRI_d<X>}{a named vector containing the discrete TRI endpoint with cutoffs defined by the <X>% quantile (may be more than 1, see \code{discretize})}
#' }
#' @seealso \code{lm}
#' @references Bucasas KL, et al. (2011) Early patterns of gene expression correlate with the humoral immune response to influenza vaccination in humans. J Infect Dis 203(7):921-9.
#' @author Stefan Avey
#' @importFrom stats lm fitted quantile median
#' @export
#' @examples
#' ## Prepare the data
#' titer_list <- FormatTiters(Year2_Titers)
#'
#' ## Calculate the titer response index (TRI)
#' endpoints <- CalculateTRI(titer_list)
#' summary(endpoints)
#'
#' ## Get discrete endpoints using upper/lower 30%
#' endpoints$TRI_d30
#'
#' ## Recreate Supp. Fig. S1
#' pairs(endpoints$scores, col = endpoints$TRI_d30)
CalculateTRI <- function(dat_list, subjectCol = "SubjectID", discretize = c(0.2, 0.3),
                         responseLabels = paste0(c("low", "moderate", "high"),
                             "Responder"), na_action = "na.fail", ...) {
  if(length(unique(lapply(dat_list, dim))) != 1) {
    stop("Each data frame in `dat_list` must have the same dimensions")
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
    ## Fit linear model
    model <- lm(data = dat, formula = "FC ~ Pre", na.action = na_action, ...)
    residuals <- residuals(model) / sd(residuals(model))
    names(residuals) <- dat[[subjectCol]]
    residuals_list[[i]] <- residuals[order(names(residuals))]
    model_list[[i]] <- model
  }
  if(length(residuals_list) > 1) {
    residual_mat <- Reduce(cbind, residuals_list)
  } else {
      residual_mat <- matrix(residuals_list[[1]], ncol = 1)
    }
  colnames(residual_mat) <- names(dat_list)
  residual_rank <- apply(residual_mat, 2, rank, na.last = "keep")
  scores <- data.frame(mean = apply(residual_mat, 1, mean, na.rm = TRUE),
                       median = apply(residual_mat, 1, median, na.rm = TRUE),
                       mean_rank = apply(residual_rank, 1, mean, na.rm = TRUE),
                       median_rank = apply(residual_rank, 1, median, na.rm = TRUE))
  scores_scaled <- scale(scores, center = TRUE, scale = TRUE)
  TRI <- apply(scores, 1, mean, na.rm = TRUE)
  names(model_list) <- names(dat_list)
  ## Calculated discretized metrics
  disList <- vector(mode = "list", length = length(discretize))
  names(disList) <- discretize
  for(dis in discretize) {
    tmp <- rep(NA, nrow(dat))
    names(tmp) <- names(TRI)
    lowR <- TRI <= quantile(TRI, dis, na.rm = TRUE)
    highR <- TRI >= quantile(TRI, 1 - dis, na.rm = TRUE)
    modR <- !(lowR | highR | is.na(TRI))
    tmp[lowR] <- 0
    tmp[modR] <- 1
    tmp[highR] <- 2
    disList[[as.character(dis)]] <- factor(tmp, labels = responseLabels, levels = 0:2)
  }
  disList <- setNames(disList, paste0("TRI_d",
                                      as.character(as.numeric(names(disList))*100)))
  return(c(models = list(model_list),
           residualMatrix = list(residual_mat),
           scores = list(scores),
           TRI = list(TRI),
           disList))
}
