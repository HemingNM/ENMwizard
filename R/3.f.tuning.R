#' Tuning and evaluation of ENMs with Maxent for several species using ENMeval
#'
#' This function is a wrapper for ENMeval::ENMevaluate. See ?ENMeval::ENMevaluate for details
#' It works with a named list of species occurrence data (occ.l) and a list of
#' cropped environmental variables (a.calib.l) for model tuning.
#' @param occ.l list of species occurrence data
#' @param a.calib.l list of predictors (cropped environmental variables) for model tuning. Used in model calibration. Argument 'x' of dismo::maxent. Raster* object or SpatialGridDataFrame, containing grids with
#' predictor variables. These will be used to extract values from for the point locations. Can
#' also be a data.frame, in which case each column should be a predictor variable and each row
#' a presence or background record.
#' @param bg.coords.l list of background localities. Two-column matrix or data.frame of longitude and latitude (in that order) of background localities (required for 'user' method).
#' @inheritParams  ENMeval::ENMevaluate
#' @examples
#' ENMeval.res.lst <- ENMevaluate.batch(occ.locs, occ.b.env, parallel = T , numCores = 7)
#' @export
ENMevaluate.batch <- function(occ.l, a.calib.l, bg.coords.l = NULL, occ.grp = NULL,
                              bg.grp = NULL, RMvalues = seq(0.5, 4, 0.5),
                              fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                              categoricals = NULL, n.bg = 10000, method = NULL,
                              overlap = FALSE, aggregation.factor = c(2, 2),
                              kfolds = NA, bin.output = FALSE, clamp = TRUE,
                              rasterPreds = TRUE, parallel = FALSE, numCores = NULL,
                              progbar = TRUE, updateProgress = FALSE, ...){

  if(is.null(bg.coords.l)) bg.coords.l <- vector("list", length(occ.l))
  if(length(bg.coords.l)!=length(occ.l)){ stop("bg.coords.l must contain one dataset of 'bg.coords' for each species")}
  ENMeval.res <- vector("list", length(occ.l))
  names(ENMeval.res) <- names(occ.l)
  # ENMeval.occ.results <- ENMeval.res
  for(i in 1:length(occ.l)){
    cat(c( "sp", i, "\n", names(occ.l)[i]), "\n")
    ENMeval.res[[i]] <- ENMeval::ENMevaluate(occ.l[[i]], a.calib.l[[i]], bg.coords = bg.coords.l[[i]],
                                             occ.grp = occ.grp, bg.grp = bg.grp, RMvalues=RMvalues,
                                    fc = fc, categoricals = categoricals, n.bg = n.bg, method = method,
                                    overlap = overlap, aggregation.factor = aggregation.factor,
                                    kfolds = kfolds, bin.output = bin.output, clamp = clamp,
                                    rasterPreds = rasterPreds, parallel = parallel, numCores = numCores,
                                    progbar = progbar, updateProgress = updateProgress)
    # ENMeval.occ.results[[i]] <- ENMeval.res[[i]]@results
  }
  return(ENMeval.res)
}