#' Titer bar plots.
#'
#' \code{Barplot} plots the baseline and day 28 titers
#'
#'
#' 
#' @param dat_list a named list like the one returned by \code{\link{FormatTiters}}.
#' @param subjectCol 
#' @param cols numeric specifying how many columns to layout plot
#' @param colors a vector of colors specifying bar colors. If \code{dat_list} contains more than 4 elements, you must specify your own colors.
#' @return (invisibly) a list of ggplot2 object(s).
#'
#' @import grid ggplot2 dplyr tidyr
#' @author Stefan Avey
#' @export
#' @examples
#' ## Prepare the data
#' strains <- c("A_California_7_2009", "A_Perth_16_2009", "B_Brisbane_60_2008")
#' titer_list <- FormatTiters(Year1_Titers, strains, subjectCol = "YaleID", otherCols = "AgeGroup")
#'
#' ## Bar plot of a single strain
#' Barplot(titer_list[strains[1]], subjectCol = "YaleID")
#'
#' ## Bar plot of all 3 strains
#' Barplot(titer_list, subjectCol = "YaleID")
#'
#' ## Can improve readability of previous plot by separating into groups
#' ## For example, group by AgeGroup
#' Barplot(titer_list, subjectCol = "YaleID", groupVar = "AgeGroup")
Barplot <- function(dat_list, subjectCol = "SubjectID", cols = 1,
                    groupVar = NULL,
                    colors =
                      c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
                        "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00"))
{
  plotList <- list()
  ## Convert list to data frame
  dat_df <- do.call(rbind.data.frame, dat_list) %>%
    mutate(strain = rep(names(dat_list),
               each = nrow(dat_list[[1]])),
           d28 = d0 + fc) %>%
    gather("condition", "titer", matches("(d0)|(28)|(fc)"))

  fcDat <- dat_df %>%
    filter(condition == "fc") %>%
    mutate(fourFC = ifelse(titer >= log2(4), TRUE, FALSE)) %>%
    select_(subjectCol, "strain", "fourFC")
  plotDat <- full_join(dat_df, fcDat, by = c(subjectCol, "strain")) %>%
    filter(condition != "fc")
  plotDat[[subjectCol]] <- factor(plotDat[[subjectCol]])
  lims <- c(min(plotDat$titer, na.rm = TRUE), max(plotDat$titer, na.rm = TRUE))
  if(!is.null(groupVar)) {
    if (!is.factor(plotDat[[groupVar]])) {
      factor(plotDat[[groupVar]])
    }
    groupLevels <- levels(plotDat[[groupVar]])
    for(group in groupLevels) {
      ## Create a plot for each group
      toKeep <- plotDat[[groupVar]] == group
      pd <- plotDat[toKeep, ]
      gg <- ggplot(pd, aes_string(x = subjectCol) +
                     aes(y = titer,
                         group = interaction(condition, strain),
                         fill = interaction(condition, strain), color = fourFC)) +
        geom_hline(aes(yintercept = log2(40)), color = "grey20", alpha = 0.5) +
        geom_bar(stat = "identity", position = "dodge") +
        coord_cartesian(ylim = lims) +
        scale_y_continuous(breaks = lims[1]:lims[2]) +
        scale_color_manual(values = c("white", "black"),
                           name = "4 Fold Change", guide = FALSE) +
        scale_fill_manual(values = colors[1:(length(dat_list)*2)], name = "Day.Strain") +
        ylab(expression("log"[2]("HAI Titer"))) +
        theme_bw() +
        theme(strip.text =element_text(size = 16),
              axis.text = element_text(size = 14),
              axis.text.x = element_text(angle = 60, hjust = 1),          
              title=element_text(size=20, face="bold"))      
      if(any(plotDat[[groupVar]] == group)) {
        gg <- gg + facet_grid(as.formula(paste("~", groupVar)),
                              scales = "free_x", drop = TRUE)
      } else {
          gg <- gg + geom_blank()
        }
      if(group != groupLevels[length(groupLevels)]) {
        gg <- gg + xlab("")
      }
      plotList[[group]] <- gg
    }
  } else {
      gg <- ggplot(plotDat,
                   aes_string(x = subjectCol) +
                     aes(y = titer,
                         group = interaction(condition, strain),
                         fill = interaction(condition, strain), color = fourFC)) +
        geom_hline(aes(yintercept = log2(40)), color = "grey20", alpha = 0.5) +
        geom_bar(stat = "identity", position = "dodge") +
        coord_cartesian(ylim = lims) +
        scale_y_continuous(breaks = lims[1]:lims[2]) +
        scale_color_manual(values = c("white", "black"),
                           name = "4 Fold Change", guide = FALSE) +
        scale_fill_manual(values = colors[1:(length(dat_list)*2)], name = "Day.Strain") +
        ylab(expression("log"[2]("HAI Titer"))) +
        theme_bw() +
        theme(strip.text =element_text(size = 16),
              axis.text = element_text(size = 14),
              axis.text.x = element_text(angle = 60, hjust = 1),          
              title=element_text(size=20, face="bold"))
      plotList[[1]] <- gg
    }
  Multiplot(plotlist = plotList, cols = cols)
  return(invisible(plotList))
}
