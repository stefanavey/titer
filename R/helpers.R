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
#' @import grid
#' @author R Cookbook
#' @references \url{http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_\%28ggplot2\%29/}
Multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {

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
#' @param m a model object
#' 
#' @references original lm_eqn and inspiration from this SO post \url{http://stackoverflow.com/questions/7549694/ggplot2-adding-regression-line-equation-and-r2-on-graph}.
#' @author Stefan Avey
#' @importFrom stats coef deviance
#' @keywords aveytoolkit
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

#################################################################################
## `+.uneval`                                                                  ##
#################################################################################
#' Addition for \code{aes()} and \code{aes_string()}
#'
#' \code{+.uneval} is a helper function to allow adding aes and aes_string in ggplot2
#'
#'
#' @param a first argument
#' @param b second argument
#' 
#' @references \url{http://stackoverflow.com/questions/28777626/how-do-i-combine-aes-and-aes-string-options}
`+.uneval` <- function(a,b) {
  `class<-`(utils::modifyList(a,b), "uneval")
}
