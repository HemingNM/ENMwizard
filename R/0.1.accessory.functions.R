### Accessory functions

#' Find repeated charaters in  a character vector
#'
#' This function will take the shorter element and the characters repeated
#' in the other elements of the character vector
#'
#' @param x Character vector
# #' @examples
# #' \dontrun{
# #' find_repeated_characters(c("aaa", "aab", "aa"))
# #' }
#' @author Neander M. Heming
#' @export
find_repeated_characters <- function(x){
  min.l <- min(sapply(x, nchar))
  for(i in 2:min.l){
    le <- length(unique(substr(x, start = 1, stop = i)))
    if(le>1){
      rptd <- unique(substr(x, start = 1, stop = i-1))
      return(rptd)
    } else if(le == min.l){
      rptd <- unique(substr(x, start = 1, stop = i-1))
      return(rptd)
    } else if(le == 1){
      rptd <- unique(substr(x, start = 1, stop = min.l))
      return(rptd)
    }
  }
}

#' Draw north arrow in plotted map
#'
#' This function will draw a simple north arrow in a plotted map
#'
#' @param xb Position at x axis
#' @param yb Position at y axis
#' @param len Arrow length
#' @param lab Label
#' @param cex.lab A numerical vector giving the amount by which plotting
#' label characters should be scaled relative to the default.
#' @param fill The color for filling polygons of both sides of the north arrow.
#' The default, c("white", border), fills one side with white and the other side
#' with the color of the border. Two colors should be provided.
#'
#' @inheritParams graphics::polygon
# #' @examples
# #' \dontrun{
# #' find_repeated_characters(c("aaa", "aab", "aa"))
# #' }
#' @author Neander M. Heming
#' @export
north_arrow <- function (xb, yb, len, lab = "NORTH", cex.lab = 1, border ="black", fill=c("white", border), ...){
  arrow.x1 = c(1, 0, 0, 1)
  arrow.x2 = arrow.x1*(-1)
  arrow.y = c(-0.5, 3.5, 0.5, -0.5)
  graphics::polygon(xb + arrow.x1 * len, yb + arrow.y * len, col=fill[1], border=border, ...)
  graphics::polygon(xb + arrow.x2 * len, yb + arrow.y * len, col=fill[2], border=border, ...)
  graphics::text(xb, yb - 1.1*graphics::strheight(lab, cex = cex.lab), lab, cex = cex.lab, col = border)
}
