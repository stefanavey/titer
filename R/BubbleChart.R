#' Bubble Chart
#'
#' \code{BubbleChart} visualizes baseline vs fold change in titers
#'
#' This plot was designed for HAI titer data with baseline columns and fold change columns for multiple strains.
#'
#' @param dat_list a named list like the one returned by \code{\link{FormatTiters}}. Values are assumed to be log2-transformed.
#' @param fit what type of fit to add. Current options are "lm" for linear model, "exp" for exponential, or \code{NULL} for no smoothing.
#' @param yMinZero a logical specifying whether fitted y values below 0 should be set to 0.
#' @param eqSize Text size of the equation. Only relevant if \code{fit} is not \code{NULL}
#' @param subjectCol the name of the column specifying a subject ID. Default is "SubjectID". 
#' @param colorBy a character string specifying an endpoint to colorBy or \code{NULL} (default) for no coloring.
#' @param xlimits the x-axis limits (passed to \code{scale_x_continuous})
#' @param xbreaks the x-axis breaks (passed to \code{scale_x_continuous})
#' @param ylimits the y-axis limits (passed to \code{scale_y_continuous})
#' @param ybreaks the y-axis breaks (passed to \code{scale_y_continuous})
#' @param plot logical indicating whether to plot or not. Default is TRUE
#' @param cols numeric specifying how many columns to layout plot
#' @param ... other arguments besides \code{method} and \code{subjectCol} passed to \code{\link{CalculatemaxRBA}}.
#' @return (invisibly) a list of ggplot2 objects.
#' 
#' @import grid ggplot2
#' @importFrom stats na.omit coef quantile
#' 
#' @author Stefan Avey
#' @keywords HIPC
#' @seealso \code{FormatTiters}
#' @export
#' @examples
#' ## Prepare the data
#' titer_list <- FormatTiters(Year2_Titers)
#'
#' ## Basic plot without any fitted model
#' BubbleChart(titer_list)
#'
#' ## Change layout to plot all in a single column
#' BubbleChart(titer_list, cols = 1)
#'
#' ## Add a linear fit
#' BubbleChart(titer_list, fit = "lm")
#' 
#' ## Add an exponential fit
#' BubbleChart(titer_list, fit = "exp")
#'
#' ## Add coloring by age
#' BubbleChart(titer_list, fit = "exp", colorBy = "AgeGroup")
BubbleChart <- function(dat_list, subjectCol = "SubjectID",
                        fit = NULL, yMinZero = FALSE,
                        eqSize = 6 / log2(length(dat_list)+1),
                        colorBy = NULL,
                        xlimits = c(1.5, 10.5), xbreaks = 2:10,
                        ylimits = c(-0.5, 10), ybreaks = seq(0, 10, 2),
                        plot = TRUE, cols = 2, ...) {
  if (sum(subjectCol == unlist(lapply(dat_list, colnames))) != length(dat_list)) {
    stop("Must specify a valid subject column name using the `subjecCol` argument")
  }
  plotList <- list()
  ## Determine upper limit for size of counts
  upLim <- max(sapply(dat_list, function(plotDat) {
                        coords <- paste(na.omit(plotDat[,"Pre"]),
                                        na.omit(plotDat[,"FC"]), sep = ',')
                        max(table(coords))
                      }))
  ## Plots for each individual strain
  for (strain in names(dat_list)) {
    plotDat <- dat_list[[strain]]
    gg <- ggplot(plotDat, aes(x = Pre, y = FC)) +
      geom_hline(aes(yintercept = log2(1)), color = "black") + 
      geom_hline(aes(yintercept = log2(4)), color = "grey20", alpha = 0.5) +
      geom_vline(aes(xintercept = log2(40)), color = "grey20", alpha = 0.5) +
      scale_size(range = c(2,7), limits = c(1, upLim)) +
      scale_y_continuous(limits = ylimits, labels = 2^ybreaks,
                         breaks = ybreaks) +      
      scale_x_continuous(limits = xlimits, labels = 2^xbreaks,
                         breaks = xbreaks, minor_breaks = NULL) +
      ## xlab(expression("log"[2]("day 0 titer"))) +
      xlab("Baseline Titer") +      
      ## ylab(expression("log"[2]("day28 / day0 titer"))) +
      ylab("Titer Fold Change") +       
      ggtitle(strain) + 
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    if (!is.null(fit)) {
      ## Fit models and save endpoints
      endpoints <- CalculatemaxRBA(dat_list, subjectCol = subjectCol,
                                    method = fit, yMinZero = yMinZero,  ...)
      mod <- endpoints$models[[strain]]
      ## Plot the exponential curve or linear fit
      if (fit == "exp") {
        ggSmooth <-  geom_smooth(method = "nls",
                                 formula = y ~ exp(a + b * x),
                                 method.args = list(start = c(a = 0, b = 0)),
                                 se = FALSE,
                                 color = "blue")
      } else if (fit == "lm") {
          xintercept <- -1 * coef(mod)["(Intercept)"] / coef(mod)["Pre"]
          if (yMinZero && (xintercept < max(plotDat$Pre))) {
            plotVals <- data.frame(Pre = c(min(plotDat$Pre), xintercept, max(plotDat$Pre)),
                                   FC = c(mod$fitted.values[1], 0, 0))
            ggSmooth <- geom_path(data = plotVals, color = "blue",
                                  size = 1.1)
          } else {
              ggSmooth <- geom_smooth(method = "lm", aes(x = Pre, y = FC),
                                      color = "blue", se = FALSE)
            }
        } else {
            stop("fit must be 'exp' or 'lm'")
          }
      ## Add text to plot with formula
      gg <- gg + geom_text(aes(x = quantile(xlimits[1]:xlimits[2], 0.50),
                               y = quantile(ylimits[1]:ylimits[2], 0.90)),
                           label = GetEqn(mod),
                           parse = TRUE, size = eqSize, color = "black")
      if (!is.null(colorBy)) { 
        if (colorBy %in% names(endpoints)) {
          plotDat[[colorBy]] <- endpoints[[colorBy]][as.character(plotDat[[subjectCol]])]
        }
      }
    }
    if (!is.null(colorBy)) {
      if (colorBy %in% colnames(plotDat)) {
        gg <- gg + geom_count(data = plotDat, mapping = aes_string(color = colorBy),
                              position = position_jitter(width = 0.2, height = 0.2))
      } else {
          stop("`colorBy` must be the name of a column in dat_list or a valid endpoint name from CalculatemaxRBA()")
        }
    } else {
        gg <- gg + geom_count()
      }
    ## Add smoothing after geom_count
    if (!is.null(fit)) {
      gg <- gg + ggSmooth
    }
    plotList[[strain]] <- gg
  }
  if (plot) {
    Multiplot(plotlist = plotList, cols = cols)
  }
  return(invisible(plotList))
}
