#' Bubble Chart
#'
#' \code{BubbleChart} visualizes baseline vs fold change in titers
#'
#' This plot was designed for HAI titer data with baseline columns and fold change columns for multiple strains.
#'
#' @param dat_list a named list like the one returned by \code{\link{FormatTiters}}
#' @param subjectCol the name of the column specifying a subject ID. Default is "SubjectID". 
#' @param colorBy a character string specifying an endpoint to colorBy or \code{NULL} (default) for no coloring.
#' @param xlimits the x-axis limits (passed to \code{scale_x_continuous})
#' @param xbreaks the x-axis breaks (passed to \code{scale_x_continuous})
#' @param plot logical indicating whether to plot or not. Default is TRUE
#' @param cols numeric specifying how many columns to layout plot
#' @param scale_y a character string specifying whether the y axis should be "fixed" for all strains or "free".
#' @param ... other arguments passed to \code{\link{CalculateSAdjMFC}}. Specifying \code{fit} will add to the plot 
#' @return (invisibly) a list of ggplot2 objects.
#' 
#' @import grid ggplot2
#' @author Stefan Avey
#' @keywords HIPC
#' @seealso \code{FormatTiters}
#' @export
#' @examples
#' ## Prepare the data
#' library(dplyr)
#' library(ggplot2)
#' titers <- filter(Year2_Titers, AgeGroup == "Young")
#' strains <- c("A_California_7_2009", "A_Perth_16_2009", "B_Brisbane_60_2008")
#' titer_list <- FormatTiters(titers, strains, subjectCol = "YaleID")
#'
#' ## Basic plot without any fitted model
#' BubbleChart(titer_list)
#'
#' ## Change layout to plot 3 strains in a single column
#' BubbleChart(titer_list, cols = 1)
#'
#' ## Add a linear fit
#' BubbleChart(titer_list, method = "lm", subjectCol = "YaleID")
#' 
#' ## Add an exponential fit
#' BubbleChart(titer_list, method = "exp", subjectCol = "YaleID")
BubbleChart <- function(dat_list, subjectCol = "SubjectID", colorBy = NULL,
                        xlimits = c(1.5, 10.5), xbreaks = 2:10,
                        plot = TRUE, cols = 2, ...) {
  plotList <- list()
  ## Determine upper limit for size of counts
  upLim <- max(sapply(dat_list, function(plotDat) {
                        coords <- paste(na.omit(plotDat[,"d0"]),
                                        na.omit(plotDat[,"fc"]), sep = ',')
                        max(table(coords))
                      }))
  ## Plots for each individual strain
  for(strain in names(dat_list)) {
    d0Col <- paste0("d0_", strain)
    fcCol <- paste0("fc_", strain)
    plotDat <- dat_list[[strain]]
    gg <- ggplot(plotDat, aes(x = d0, y = fc)) +
      geom_hline(aes(yintercept = log2(1)), color = "black") + 
      geom_hline(aes(yintercept = log2(4)), color = "grey20", alpha = 0.5) +
      geom_vline(aes(xintercept = log2(40)), color = "grey20", alpha = 0.5) +
      scale_size(range = c(2,7), limits = c(1, upLim)) +
      scale_x_continuous(limits = xlimits, breaks = xbreaks) +
      xlab(expression("log"[2]("day 0 titer"))) +
      ylab(expression("log"[2]("day28 / day0 titer"))) + 
      ggtitle(strain) + 
      theme_bw()
    if(!is.null(fit)) {
      ## Plot the exponential curve or linear fit
      if(fit == "exp") {
        gg <- gg + geom_smooth(method = "nls",
                               formula = y ~ exp(a + b * x),
                               method.args = list(start = c(a = 0, b = 0)),
                               se = FALSE,
                               color = "blue")
      } else if (fit == "lm") {
          gg <- gg + geom_smooth(method = "lm", aes(x = d0, y = fc), color = "blue")
        } else {
            stop("fit must be 'exp' or 'lm'")
          }
      ## Add text to plot with formula
      ## Calling this function is an easy way to get the formulas
      endpoints <- CalculateSAdjMFC(dat_list, subjectCol = subjectCol, ...)
      mod <- endpoints$models[[strain]]
      gg <- gg + geom_text(aes(x = 6,#mean(unique(d0), na.rm = TRUE),
                               y = quantile(unique(fc), 0.95)),
                           label = GetEqn(mod),
                           parse = TRUE, size = 2.5, color = "black")
      if (!is.null(colorBy)) {
        if (colorBy %in% names(endpoints)) {
          plotDat[[colorBy]] <- endpoints[[colorBy]][as.character(plotDat[[subjectCol]])]
        }
        if (colorBy %in% colnames(plotDat)) {
          gg <- gg + geom_count(data = plotDat, mapping = aes_string(color = colorBy),
                                position = position_jitter(width = 0.2, height = 0.2))
        } else {
            stop("`colorBy` must be the name of a column in dat_list or a valid endpoint name from CalculateSAdjMFC()")
          }
      } else {
          gg <- gg + geom_count()
        }
    } else {
        gg <- gg + geom_count()
      }
    plotList[[strain]] <- gg
  }
  if(plot) {
    Multiplot(plotlist = plotList, cols = cols)
  }
  return(invisible(plotList))
}
