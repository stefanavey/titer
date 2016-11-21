##' @title BubbleChart
##'
##' @description
##' \code{BubbleChart} wraps \code{ggplot2} to plot baseline vs fold changes
##'
##' @param dat_list a list like the one returned by \code{FormatTiters}
##' @param fit what type of fit to add. Current options are "lm" for linear model, "exp" for exponential, or \code{NULL} for no smoothing.
##' @param xlimits the x-axis limits (passed to scale_x_continuous)
##' @param xbreaks the x-axis breaks (passed to scale_x_continuous)
##' @param plot logical indicating whether to plot or not. Default is TRUE
##'
##' @param cols numeric specifying how many columns to layout plot
##' @param scale_y a character string specifying whether the y axis should be "fixed" for all strains or "free".
##' @details
##' This plot was designed for HAI titer data with baseline columns and fold change columns for multiple strains.
##'
##' @return a list of ggplot2 objects.
##'
##' 
##' @import grid ggplot2
##' @importFrom aveytoolkit Multiplot GetEqn
##' @author Stefan Avey
##' @keywords HIPC
##' @seealso \code{FormatTiters}
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
##' BubbleChart(titer_list)
BubbleChart <- function(dat_list, fit = NULL,
                        xlimits = c(1.5, 10.5), xbreaks = 2:10,
                        plot = TRUE, cols = 2) {
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
      ## geom_count(position = position_jitter(width = 0.2, height = 0.2)) +
      geom_count() +
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
      endpoints <- CalculateSAdjMFC(plotDatList, method = fit)
      mod <- endpoints$models[[strain]]
      gg <- gg + geom_text(aes(x = 6,#mean(unique(d0), na.rm = TRUE),
                               y = quantile(unique(fc), 0.95)),
                           label = GetEqn(mod),
                           parse = TRUE, size = 2.5, color = "black")
    }
    plotList[[strain]] <- gg
  }
  if(plot) {
    Multiplot(plotlist = plotList, cols = cols)
  }
  return(invisible(plotList))
}
