#' Tuning and evaluation of ENMs with Maxent for several species using ENMeval
#'
#' This function is a wrapper for ENMeval::ENMevaluate. See ?ENMeval::ENMevaluate for details
#' It works with a named list of species occurrence data (occ_locs) and a list of
#' cropped environmental variables (occ_b_env) for model tuning.
#' @param occ_locs list of species occurrence data
#' @param occ_b_env list of cropped environmental variables for model tuning
#' @inheritParams  ENMeval::ENMevaluate
#' @examples
#' ENMeval_res.lst <- ENMevaluate.batch(occ_locs, occ_b_env, parallel = T , numCores = 7)
#' @export
ENMevaluate.batch <- function(occ_locs, occ_b_env, bg.coords.l = NULL, occ.grp = NULL,
                              bg.grp = NULL, RMvalues = seq(0.5, 4, 0.5),
                              fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                              categoricals = NULL, n.bg = 10000, method = NULL,
                              overlap = FALSE, aggregation.factor = c(2, 2),
                              kfolds = NA, bin.output = FALSE, clamp = TRUE,
                              rasterPreds = TRUE, parallel = FALSE, numCores = NULL,
                              progbar = TRUE, updateProgress = FALSE, ...){

  if(is.null(bg.coords.l)) bg.coords.l <- vector("list", length(occ_locs))
  if(length(bg.coords.l)!=length(occ_locs)){ stop("bg.coords.l must contain one dataset of 'bg.coords' for each species")}
  ENMeval_res <- vector("list", length(spp.occ.list))
  names(ENMeval_res) <- names(spp.occ.list)
  # ENMeval_occ_results <- ENMeval_res
  for(i in 1:length(spp.occ.list)){
    cat(c( "sp", i, "\n", names(spp.occ.list)[i]), "\n")
    ENMeval_res[[i]] <- ENMeval::ENMevaluate(occ_locs[[i]], occ_b_env[[i]], bg.coords = bg.coords.l[[i]],
                                             occ.grp = occ.grp, bg.grp = bg.grp, RMvalues=RMvalues,
                                    fc = fc, categoricals = categoricals, n.bg = n.bg, method = method,
                                    overlap = overlap, aggregation.factor = aggregation.factor,
                                    kfolds = kfolds, bin.output = bin.output, clamp = clamp,
                                    rasterPreds = rasterPreds, parallel = parallel, numCores = numCores,
                                    progbar = progbar, updateProgress = updateProgress)
    # ENMeval_occ_results[[i]] <- ENMeval_res[[i]]@results
  }
  return(ENMeval_res)
}
