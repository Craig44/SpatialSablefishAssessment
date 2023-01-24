#' Vertically align a list of plots.
#'
#' @description This function aligns the given list of plots so that the x axis are aligned.
#' It assumes that the graphs share the same range of x data.
#'
#' @param ... The list of plots to align.
#' @param globalTitle The title to assign to the newly created graph.
#' @param keepTitles TRUE if you want to keep the titles of each individual
#' plot.
#' @param keepXAxisLegends TRUE if you want to keep the x axis labels of each
#' individual plot. Otherwise, they are all removed except the one of the graph
#' at the bottom.
#' @param nb.columns The number of columns of the generated graph.
#' @details from https://stackoverflow.com/questions/41569817/align-multiple-plots-in-ggplot2-when-some-have-legends-and-others-dont
#' @return The gtable containing the aligned plots.
#' @export
#' @examples
#' x = seq(0, 10, length.out = 200)
#' y1 = sin(x)
#' y2 = cos(x)
#' y3 = sin(x) * cos(x)
#' df1 <- data.frame(x, y1, y2)
#' df1 <- reshape2::melt(df1, id.vars = "x")
#' g1 <- ggplot2::ggplot(df1, ggplot2::aes(x, value, color = variable)) + ggplot2::geom_line()
#' df2 <- data.frame(x, y3)
#' g2 <- ggplot2::ggplot(df2, ggplot2::aes(x, y3)) + ggplot2::geom_line()
#' g <- VAlignPlots(g1, g2, globalTitle = "Alignment test")
#' grid::grid.newpage()
#' grid::grid.draw(g)
VAlignPlots <- function(...,
                        globalTitle = "",
                        keepTitles = FALSE,
                        keepXAxisLegends = FALSE,
                        nb.columns = 1) {
  # Retrieve the list of plots to align
  plots.list <- list(...)

  # Remove the individual graph titles if requested
  if (!keepTitles) {
    plots.list <- lapply(plots.list, function(x) x <- x + ggtitle(""))
    plots.list[[1]] <- plots.list[[1]] + ggtitle(globalTitle)
  }

  # Remove the x axis labels on all graphs, except the last one, if requested
  if (!keepXAxisLegends) {
    plots.list[1:(length(plots.list)-1)] <-
      lapply(plots.list[1:(length(plots.list)-1)],
             function(x) x <- x + theme(axis.title.x = element_blank()))
  }

  # Builds the grobs list
  grobs.list <- lapply(plots.list, ggplotGrob)

  #max_width(grobs.list, value_only = FALSE)
  #unit.pmax(lapply(grobs.list, "[[", 'widths')[[1]],lapply(grobs.list, "[[", 'widths')[[2]])

  # Get the max width
  all_widths = lapply(grobs.list, "[[", 'widths')
  max_widths = all_widths[[1]]
  for(i in 1:length(all_widths[[1]])) {
    temp = do.call(unit.pmax,lapply(all_widths, FUN = function(x) {x[[i]]}))
    max_widths[1] = temp
  }
  #widths.list <- do.call(unit.pmax, lapply(grobs.list, "[[", 'widths'))

  #widths.list <- do.call(max_width, lapply(grobs.list, FUN = max_width)

  #grobWidth(grobs.list)

  # Assign the max width to all grobs
  grobs.list <- lapply(grobs.list, function(x) {
    x[['widths']] = max_widths
    x})

  # Create the gtable and display it
  g <- grid.arrange(grobs = grobs.list, ncol = nb.columns)
  # An alternative is to use arrangeGrob that will create the table without
  # displaying it
  #g <- do.call(arrangeGrob, c(grobs.list, ncol = nb.columns))

  return(g)
}
