% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/7.f.plot.mdl.sel.crit.diff.R
\name{plot_mdl}
\alias{plot_mdl}
\title{internal function for \code{\link{plot_mdl_diff}}, \code{\link{plot_scn_diff_b}}}
\usage{
plot_mdl(
  sc,
  mcmp,
  mtp,
  basemap,
  thrshld.path,
  thrshld.nms.mod,
  save,
  mfrow = c(0, 0),
  msnm,
  msr,
  scnm,
  scr
)
}
\arguments{
\item{sc}{Index of climatic scenario to be plotted}

\item{mtp}{Stack or brick of thresholded predictions}

\item{basemap}{Shapefile to be plotted with. Usually a continent or country shapefile}

\item{thrshld.path}{path to threshold projections}

\item{thrshld.nms.mod}{names of threshold models}

\item{save}{Logical. If TRUE will save plots in pdf.}

\item{mfrow}{Number of rows and cols to draw plot (see ?par for details). If}

\item{msnm}{Character vector. Short names of model selection criteria to be replaced. Same as given in
Model names (ex. AUC10, AUCmtp, OR10, ORmtp)}

\item{msr}{Character vector. Long names of model selection criteria to replace the short names. Must
be in same order of argument 'msnm'}

\item{scnm}{Character vector. Short names of climatic scenarios to be replaced. Ex. "cc26bi70", "cc45bi70"}

\item{scr}{Character vector. Long names of climatic scenarios to replace the short names. Must
be in same order of argument 'scnm'. Ex. "2070-CCSM4-rcp2.6", "2070-CCSM4-rcp4.5"}
}
\description{
internal function for \code{\link{plot_mdl_diff}}, \code{\link{plot_scn_diff_b}}
}
\keyword{internal}
