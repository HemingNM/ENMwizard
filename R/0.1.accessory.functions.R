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
#' @param x Position at x axis
#' @param y Position at y axis
#' @param size Arrow size
#' @param lab Label
#' @param cex.lab A numerical vector giving the amount by which plotting
#' label characters should be scaled relative to the default.
#' @param fill The color for filling polygons of both sides of the north arrow or compass rose.
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
north_arrow <- function (x, y, size, lab = "North", cex.lab = 1, border ="black", fill=c("white", border), ...){
  arrow.x1 = c(1, 0, 0, 1)
  arrow.x2 = arrow.x1*(-1)
  arrow.y = c(-0.5, 3.5, 0.5, -0.5)
  graphics::polygon(x + arrow.x1 * size, y + arrow.y * size, col=fill[1], border=border, ...)
  graphics::polygon(x + arrow.x2 * size, y + arrow.y * size, col=fill[2], border=border, ...)
  graphics::text(x, y - 1.1*graphics::strheight(lab, cex = cex.lab), lab, cex = cex.lab, col = border)
}


#' Draws circle over compass rose in plotted map
#'
#' @param radius circle radius
#' @param nv number of 'vertices', points that maken the circle
#' @inheritParams compass_rose
#' @inheritParams graphics::polygon
#'
#' @keywords internal
circle4compass_rose <- function (x, y, radius=c(1, 1.3), nv = 200, bearing=0, border = NULL, fill = NA, lty = 1,
                               density = NULL, angle = 45, lwd = 1) {
  nv <- round(nv/16, 0)*16
  xylim <- par("usr")
  plotdim <- par("pin")
  ymult <- plotrix::getYmult()
  angle.inc <- 2 * pi/nv
  angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
  if (length(fill) < length(radius))
    fill <- rep(fill, length.out = length(radius))
  xv <- matrix(nrow = nv, ncol = length(radius))# numeric(length(radius)*nv)
  yv <- xv
  for (circle in 1:length(radius)) {
    xv[,circle] <- cos(angles + bearing) * radius[circle] + x
    yv[,circle] <- sin(angles + bearing) * radius[circle] * ymult + y
  }
  gap <- nv/16 - 1
  nsq <- (1:nv)[c(T, rep(F, gap))]
  nsq <- nsq[c(F,T)]
  polygon(c(xv[1:nv,1], xv[nv:1,2], xv[1,1]) , c(yv[1:nv,1], yv[nv:1,2], yv[1,1]),
          border = border, col = fill[1], lty = lty,
          density = density, angle = angle, lwd = lwd)
  for(i in nsq){
    polygon(c(xv[c(i:(i+gap)),1], xv[c((i+gap):i),2], xv[i,1]) , c(yv[c(i:(i+gap)),1], yv[c((i+gap):i),2], yv[i,1]) ,
            border = border, col = fill[2], lty = lty,
            density = density, angle = angle, lwd = lwd)
  }
  invisible(list(x = xv, y = yv))
}

#' Draw north arrow in plotted map
#'
#' This function will draw a compass rose in the plotted map
#'
#' @param size Compass size
#' @param bearing bearing direction in degrees
#' @param d
#' @param draw.circle Logical. Draw circle around compass?
#' @inheritParams north_arrow
#' @inheritParams graphics::polygon
#' @references Tanimura et al. 2007. Auxiliary Cartographic Functions in R: North
#' Arrow, Scale Bar, and Label with a Leader Arrow. Journal of Statistical
#' Software. Volume 19, Code Snippet 1. https://www.jstatsoft.org/article/view/v019c01/v19c01.pdf
#'
# https://gis.stackexchange.com/questions/243743/alternative-north-arrow-using-r-gistools-package
compass_rose <- function(x, y, size, bearing=0, d=.8, draw.circle=F, fill = c("white","black"), cex=1, lwd=1,...) {
  bearing <- bearing*pi/180
  # checking arguments
  if(missing(x) | missing(y)) stop("x or y is missing")
  if(missing(size)) stop("size is missing")
  # default colors are white and black
  if(length(fill)!=2) stop("two colors must be provided")
  fill <- rep(fill, 8)
  loc <- c(x, y)
  ymult <- plotrix::getYmult()

  if(draw.circle){
    # plot(mods.thrshld.lst$Mole[[2]]$continuous[[1]], col="gray" )
    circle4compass_rose(loc[1], loc[2], radius=c(size*1.45, size*1.65), nv = 192, bearing=bearing, border = NULL, fill = fill[2:1], lty = 1,
                      density = NULL, angle = 45, lwd = lwd)

    # calculating coordinates of polygons
    # small triangles
    radiit <- rep(size/c(1.25, 1.25, 1.25, 1.15)*2.15, 8)
    pos <- (c(0:2,1)/2.2-2.42)+rep(seq(0, 28, 4), each=4)
    xt <- radiit*cos((pos)*pi/16+bearing) + loc[1]
    yt <- radiit*sin((pos)*pi/16+bearing) * ymult + loc[2]

    for (i in seq(1,32, 4)) {
      x1 <- c(xt[i], xt[i+1], xt[i+3])
      y1 <- c(yt[i], yt[i+1], yt[i+3])
      polygon(x1,y1,col=fill[1], lwd=lwd*.8)

      x1 <- c(xt[i+2], xt[i+1], xt[i+3])
      y1 <- c(yt[i+2], yt[i+1], yt[i+3])
      polygon(x1,y1,col=fill[2],lwd=lwd*.8)
    }
  }

  # calculating coordinates of main polygons
  radii <- rep(size/c(1.2,3,1.75,3)*2.2,4) # c(1,4,2,4) # c(1,3,1.8,3)
  x <- radii[(0:15)+1]*cos((0:15)*pi/8+bearing) + loc[1]
  y <- radii[(0:15)+1]*sin((0:15)*pi/8+bearing) * ymult + loc[2]
  # drawing polygons
  for (i in 1:16) {
    i2 <- ifelse(i+1<=16, i+1, 1)
    x1 <- c(x[i],x[i2],loc[1])
    y1 <- c(y[i],y[i2],loc[2])
    polygon(x1,y1,col=fill[i],lwd=lwd)
  }
  # drawing letters
  b <- c("E","N","W","S")
  for (i in 0:3) text((size/1.2*2.2*d+par("cxy")[1])*cos(bearing+i*pi/2)+loc[1],
                      (size/1.2*2.2*d+par("cxy")[1])*sin(bearing+i*pi/2) * ymult + loc[2],b[i+1],
                      cex=cex)
}
