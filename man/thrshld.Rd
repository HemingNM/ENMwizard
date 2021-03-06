% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/6.f.threshold.R
\name{thrshld}
\alias{thrshld}
\title{Apply threshold for MaxEnt projections of a species}
\usage{
thrshld(mcmp, thrshld.i = 4:6, t.all = FALSE, sp.nm = "species", numCores = 1)
}
\arguments{
\item{mcmp}{An object returned by \code{\link{proj_mdl}}, containing calibrated models and model
projections for each species; or species "i" of a object returned by "proj_mdl_b", containing a list of
calibrated models and model projections for each species}

\item{thrshld.i}{List of threshold criteria to be applied. Use numbers to choose the desired one(s). Current options:
1. Fixed.cumulative.value.1 (fcv1);
2. Fixed.cumulative.value.5 (fcv5);
3. Fixed.cumulative.value.10 (fcv10);
4. Minimum.training.presence (mtp);
5. 10.percentile.training.presence (x10ptp);
6. Equal.training.sensitivity.and.specificity (etss);
7. Maximum.training.sensitivity.plus.specificity (mtss);
8. Balance.training.omission.predicted.area.and.threshold.value (bto);
9. Equate.entropy.of.thresholded.and.original.distributions (eetd).}

\item{t.all}{logical. Should threshold be applied on individual and consensus projections?
Default is FALSE. Ignored if consensus projections are not found.}

\item{sp.nm}{Species name. Used to name the output folder}

\item{numCores}{Number of cores to use for parallelization. If set to 1, no paralellization is performed}
}
\value{
Stack or brick of thresholded predictions
}
\description{
This function will apply the selected threshold criterias to MaxEnt model projection(s) of a 'mcmp' object
and save on the folder "3_out.MaxEnt/Mdls.[species name]/Mdls.thrshld". For each projection (species and climatic
scenario), two layers will be generated, one with suitability above the threshold value and another with presence/absence only.
}
\examples{
\dontrun{
mods.thrshld <- thrshld(mcmp=mxnt.mdls.preds, thrshld.i = 4:6)
plot(mods.thrshld[[1]][[2]]) # continuous
plot(mods.thrshld[[2]][[2]]) # binary
}
}
\seealso{
\code{\link{thrshld_b}}
}
