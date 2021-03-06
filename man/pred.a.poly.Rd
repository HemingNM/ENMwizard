% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/xx.Deprecated.functions.R
\name{pred.a.poly}
\alias{pred.a.poly}
\title{Select area for projection based on the extent of occ points}
\usage{
pred.a.poly(...)
}
\arguments{
\item{...}{additional arguments}
}
\description{
This function will create SpatialPolygons that will be used to crop/mask raster/brick objects to be used on model projections. It has several options
The user can crop a squared area with an extent larger than the extent of occ.poly. By default, the "extent increase"
is the maximum of latitudinal and longitudinal extent "max(abs(ext.proj[1]-ext.proj[2]), abs(ext.proj[3]-ext.proj[4]))".
The result is added to each side of the occ.poly extent. This may be changed by setting "mult", which will be multiplied
by the "extent increase". Latitudinal and longitudinal increase may also vary independently by setting "same=FALSE".
}
\details{
The user can also set a buffer around occ.poly to cut the raster/brick object. Buffer value is, by default,
calculated in the same way as "extent increase", using "max(abs(ext.proj[1]-ext.proj[2]), abs(ext.proj[3]-ext.proj[4]))".
But the exact value can be defined through "deg.incr" and "mult". This method takes longer to run.
}
