#' Custom `ggplot2` theme
#'
#' A custom `ggplot2` theme (based on `ggplot2::theme_bw()`).
#'
#' @param axis_col The colour of axis lines, ticks and labels.
#' @examples
#' theme_custom()
#' @export
theme_custom <- function (axis_col = "black") { 
  ggplot2::theme_bw() %+replace% 
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(1, 1, 1, 1),
      axis.line = ggplot2::element_line(colour = axis_col, size = 0.25),
      axis.ticks = ggplot2::element_line(colour = axis_col, size = 0.25),
      axis.text = ggplot2::element_text(colour = axis_col, size = 8),
      axis.title = ggplot2::element_text(colour = axis_col, size = 8),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = 11),
      strip.placement = "outside",
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = ggplot2::unit(0, "cm"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 8)
    )
}