% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/0.1.accessory.functions.R
\name{circle4compass_rose}
\alias{circle4compass_rose}
\title{Draws circle over compass rose in plotted map}
\usage{
circle4compass_rose(
  x,
  y,
  radius = c(1, 1.3),
  nv = 200,
  bearing = 0,
  border = NULL,
  fill = NA,
  lty = 1,
  density = NULL,
  angle = 45,
  lwd = 1
)
}
\arguments{
\item{x}{Position at x axis}

\item{y}{Position at y axis}

\item{radius}{circle radius}

\item{nv}{number of 'vertices', points that maken the circle}

\item{bearing}{bearing direction in degrees}

\item{border}{the color to draw the border.  The default, \code{NULL},
    means to use \code{\link[graphics]{par}("fg")}.  Use \code{border = NA} to
    omit borders.

    For compatibility with S, \code{border} can also be logical, in
    which case \code{FALSE} is equivalent to \code{NA} (borders omitted)
    and \code{TRUE} is equivalent to \code{NULL} (use the foreground colour),
  }

\item{fill}{The color for filling polygons of both sides of the north arrow or compass rose.
The default, c("white", border), fills one side with white and the other side
with the color of the border. Two colors should be provided.}

\item{lty}{the line type to be used, as in \code{\link[graphics]{par}}.}

\item{density}{the density of shading lines, in lines per inch.  The
    default value of \code{NULL} means that no shading lines are drawn.
    A zero value of \code{density} means no shading nor filling whereas
    negative values and \code{NA} suppress shading (and so allow
    color filling).}

\item{angle}{the slope of shading lines, given as an angle
     in degrees (counter-clockwise).}

\item{lwd}{A numerical value giving the amount by which plotting text
and symbols should be magnified relative to the default.
See \code{\link[graphics]{par}} for details}
}
\description{
Draws circle over compass rose in plotted map
}
\keyword{internal}
