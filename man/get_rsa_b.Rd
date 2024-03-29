% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/8.1.f.metrics-area_gain_loss.R
\name{get_rsa_b}
\alias{get_rsa_b}
\title{Compute changes in total suitable area for multiple species}
\usage{
get_rsa_b(mrs.l, area.raster = NULL, digits = 2, numCores = 1)
}
\arguments{
\item{mrs.l}{List containing stack or brick of range shifts. See \code{\link{range_shift_b}}}

\item{area.raster}{A raster containing the cell areas to be summed across
the suitable pixels. This allows summing areas of habitat when the pixel is
partially occupied with the habitat of interest.}

\item{digits}{integer indicating the number of decimal places. see ?round for details.}

\item{numCores}{Number of cores to use for parallelization. If set to 1, no paralellization is performed}
}
\value{
An object of class "data.table" "data.frame" containing calculated areas of
 range shift for each species between a selected climatic scenario and all other climatic scenarios
 for threshold, model and masked locations (if projections were masked,
 see \code{\link{mask_thr_projs_mscn_b}}).
}
\description{
Compute changes in total suitable area at multiple climatic scenarios,
threshold and model criteria for each species.
}
\examples{
\dontrun{
spp_diff_a <- get_rsa_b(mods.thrshld.lst)
}
}
\seealso{
\code{\link{range_shift}}, \code{\link{get_rsa}}
}
