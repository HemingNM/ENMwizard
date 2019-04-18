#' Tuning and evaluation of ENMs with Maxent for several species using ENMeval
#'
#' This function is a wrapper for ENMeval::ENMevaluate. See ?ENMeval::ENMevaluate for details
#' It works with a named list of species occurrence data (occ.l) and a list of
#' cropped environmental variables (a.calib.l) for model tuning.
#'
#' @param occ.l list of species occurrence data
#' @param a.calib.l list of predictors (cropped environmental variables) for model tuning. Used in model calibration. Argument 'x' of dismo::maxent. Raster* object or SpatialGridDataFrame, containing grids with
#' predictor variables. These will be used to extract values from for the point locations. Can
#' also be a data.frame, in which case each column should be a predictor variable and each row
#' a presence or background record.
#' @param bg.coords.l list of background localities. Two-column matrix or data.frame of longitude and latitude (in that order) of background localities (required for 'user' method).
#' @param bg.grp.l list containing a vector of bins of occurrence localities (required for 'user' method) for each species.
#' @param occ.grp.l list containing a vector of bins of occurrence localities (required for 'user' method) for each species.
#' @param resultsOnly logical; If TRUE, only results, occ.pts, and bg.pts are returned.
#' 'predictions', 'models', 'occ.grp', and 'bg.grp' are set to NULL.
#' It will not be possible to check MaxEnte models, predictions and grouping of occ and bg points.
#' Can be used to optimize allocated RAM memory when 'ENMevaluate' objects are too large.
#' @inheritParams  ENMeval::ENMevaluate
#'
#' @seealso \code{\link[ENMeval]{ENMevaluate}}
#' @examples
#' ENMeval.res.lst <- ENMevaluateB(occ.locs, occ.b.env, parallel = T , numCores = 7)
#' @export
ENMevaluateB <- function(occ.l, a.calib.l, bg.coords.l = NULL,
                         occ.grp.l = NULL, bg.grp.l = NULL,
                         RMvalues = seq(0.5, 4.5, 0.5),
                              fc = c("L", "P", "Q", "H",
                                     "LP", "LQ", "LH",
                                     "PQ", "PH", "QH",
                                     "LPQ", "LPH", "LQH", "PQH",
                                     "LPQH"),
                              categoricals = NULL, n.bg = 10000, method = "block",
                              algorithm = 'maxnet', overlap = FALSE, aggregation.factor = c(2, 2),
                              kfolds = NA, bin.output = FALSE, clamp = TRUE,
                              rasterPreds = TRUE, parallel = FALSE, numCores = NULL,
                              progbar = TRUE, updateProgress = FALSE,
                              resultsOnly = F, ...){

  if(is.null(bg.coords.l)) bg.coords.l <- vector("list", length(occ.l))
  if(length(bg.coords.l)!=length(occ.l)){ stop("bg.coords.l must contain one dataset of 'bg.coords' for each species")}
  ENMeval.res <- vector("list", length(occ.l))
  names(ENMeval.res) <- names(occ.l)

  if(resultsOnly){
    for(i in 1:length(occ.l)){
      cat(c( "sp", i, "\n", names(occ.l)[i]), "\n")
      eval <- ENMeval::ENMevaluate(occ.l[[i]], a.calib.l[[i]], bg.coords = bg.coords.l[[i]],
                                   occ.grp = occ.grp.l[[i]], bg.grp = bg.grp.l[[i]], RMvalues=RMvalues,
                                   fc = fc, categoricals = categoricals, n.bg = n.bg, method = method,
                                   algorithm = algorithm, overlap = overlap, aggregation.factor = aggregation.factor,
                                   kfolds = kfolds, bin.output = bin.output, clamp = clamp,
                                   rasterPreds = rasterPreds, parallel = parallel, numCores = numCores,
                                   progbar = progbar, updateProgress = updateProgress)
      ENMeval.res[[i]] <- list( # methods::new("ENMevaluate.opt",
                              results = eval@results,
                              occ.pts = eval@occ.pts,
                              bg.pts = eval@bg.pts)
    }
  } else {
      for(i in 1:length(occ.l)){
        ENMeval.res[[i]] <- ENMeval::ENMevaluate(occ.l[[i]], a.calib.l[[i]], bg.coords = bg.coords.l[[i]],
                                                 occ.grp = occ.grp.l[[i]], bg.grp = bg.grp.l[[i]], RMvalues=RMvalues,
                                                 fc = fc, categoricals = categoricals, n.bg = n.bg, method = method,
                                                 algorithm = algorithm, overlap = overlap, aggregation.factor = aggregation.factor,
                                                 kfolds = kfolds, bin.output = bin.output, clamp = clamp,
                                                 rasterPreds = rasterPreds, parallel = parallel, numCores = numCores,
                                                 progbar = progbar, updateProgress = updateProgress)
      }
    }
  return(ENMeval.res)
}



#' Optimize the size of ENMevaluateB objects
#'
#' This function will set to NULL (erase) the largest slots ('predictions', 'models', 'occ.grp', and 'bg.grp')
#' of ENMeval::ENMevaluate objects. Only results, occ.pts, and bg.pts are returned.
#' Use with care. It will not be possible to check MaxEnte models, predictions and grouping of occ and bg points.
#' Can be used to optimize allocated RAM memory when 'ENMevaluate' objects are too large.
#'
#' @inheritParams mxntCalibB
#'
#' @seealso \code{\link{ENMevaluateB}}, \code{\link[ENMeval]{ENMevaluate}},
#' @keywords internal
#' @export
optENMevalObjL <- function (ENMeval.o.l) {
  ENMeval.res <- vector("list", length(ENMeval.o.l))
  names(ENMeval.res) <- names(ENMeval.o.l)
  for (i in seq_along(ENMeval.o.l)) {
    ENMeval.res[[i]] <- list( # methods::new("ENMevaluate.opt",
                             results = ENMeval.o.l[[i]]@results,
                             occ.pts = ENMeval.o.l[[i]]@occ.pts,
                             bg.pts = ENMeval.o.l[[i]]@bg.pts)
  }
  return(ENMeval.res)
}


