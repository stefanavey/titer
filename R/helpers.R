#################################################################################
## Multiplot                                                                   ##
#################################################################################
#' Multiple ggplot2 plots on the same page
#'
#' Multiple Plot Function for ggplot
#'
#' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right,
#' and 3 will go all the way across the bottom.
#'
#' @param ... ggplot objects
#' @param plotlist a list of ggplot objects
#' @param cols Number of columns in layout
#' @param layout A matrix specifying the layout. If present, 'cols' is ignored
#' @importFrom grid grid.newpage grid.layout pushViewport viewport
#' @author R Cookbook
#' @references \url{http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_\%28ggplot2\%29/}
#' @examples
#' library(ggplot2)
#' 
#' ## This example uses the ChickWeight dataset, which comes with ggplot2
#' ## First plot
#' p1 <- ggplot(ChickWeight, aes(x=Time, y=weight, colour=Diet, group=Chick)) +
#'   geom_line() +
#'   ggtitle("Growth curve for individual chicks")
#' 
#'                                         # Second plot
#' p2 <- ggplot(ChickWeight, aes(x=Time, y=weight, colour=Diet)) +
#'   geom_point(alpha=.3) +
#'   geom_smooth(alpha=.2, size=1) +
#'   ggtitle("Fitted growth curve per diet")
#' 
#'                                         # Third plot
#' p3 <- ggplot(subset(ChickWeight, Time==21), aes(x=weight, colour=Diet)) +
#'   geom_density() +
#'   ggtitle("Final weight, by diet")
#' 
#'                                         # Fourth plot
#' p4 <- ggplot(subset(ChickWeight, Time==21), aes(x=weight, fill=Diet)) +
#'   geom_histogram(colour="black", binwidth=50) +
#'   facet_grid(Diet ~ .) +
#'   ggtitle("Final weight, by diet") +
#'   theme(legend.position="none")        # No legend (redundant in this graph)
#' 
#' Multiplot(p1, p2, p3, p4, cols=2)
Multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {

  ## Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  ## If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    ## Make the panel
    ## ncol: Number of columns of plots
    ## nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
      ## Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

      ## Make each plot, in the correct location
      for (i in 1:numPlots) {
        ## Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                              layout.pos.col = matchidx$col))
      }
    }
}

#################################################################################
## GetEqn                                                                      ##
#################################################################################
#' Get Formatted Model Equation
#'
#' \code{GetEqn} gets the equation for various models in a human readable format
#'
#'
#' 
#' @references original lm_eqn and inspiration from this SO post \url{http://stackoverflow.com/questions/7549694/ggplot2-adding-regression-line-equation-and-r2-on-graph}.
#' @author Stefan Avey
#' @keywords aveytoolkit
#' @examples
#' ## First Example
#' 
#' @param m a model object
GetEqn <- function(m)
{
  .lm_eqn <- function(m) {
    l <- list(a = format(coef(m)[1], digits = 2),
              b = format(abs(coef(m)[2]), digits = 2),
              r2 = format(summary(m)$r.squared, digits = 2));
    
    if (coef(m)[2] >= 0)  {
      eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
    } else {
        eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
      }
    as.character(as.expression(eq));
  }

  .nls_eqn <- function(m) {
    l <- list(a = format(coef(m)[1], digits = 2),
              b = format(abs(coef(m)[2]), digits = 2),
              rss = format(deviance(m), digits = 2))

    if (coef(m)[2] >= 0)  {
      eq <- substitute(italic(y) ==~expr(a + b %.% italic(x))*", RSS ="~rss,l)
    } else {
        eq <- substitute(italic(y) ==~exp(a - b %.% italic(x))*", RSS ="~rss,l)
      }
    
    as.character(as.expression(eq));
  }


  if(class(m) == "nls")
    .nls_eqn(m)
  else if (class(m) == "lm")
    .lm_eqn(m)
}
