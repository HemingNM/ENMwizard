#### 4.3.3 aplicar threshold
# TODO examples
# name of arg "mxnt.mdls.preds.sp[...]" shortened to "mcmp"

#' Apply threshold for MaxEnt projections of a species
#'
#' This function will apply the selected threshold criterias to MaxEnt model projection(s) of a 'mcmp' object
#' and save on the folder "3_out.MaxEnt/Mdls.[species name]/Mdls.thrshld". For each projection (species and climatic
#' scenario), two layers will be generated, one with suitability above the threshold value and another with presence/absence only.
#'
#' @inheritParams calib_mdl
#' @param mcmp Species "i" of a object returned by "proj_mdl_b", containing a list of
#' calibrated models and model projections for each species
#' @param scn.nm Name of climatic scenario to be looked for
# #' @param path.mdls Path where thresholded rasters will be saved
# #' @param pred.nm name of prediction to be appended to the final name. Usually "pres", "past" or "fut".
#' @param thrshld.i List of threshold criteria to be applied. Use numbers to choose the desired one(s). Current options:
#' 1. Fixed.cumulative.value.1 (fcv1);
#' 2. Fixed.cumulative.value.5 (fcv5);
#' 3. Fixed.cumulative.value.10 (fcv10);
#' 4. Minimum.training.presence (mtp);
#' 5. 10.percentile.training.presence (x10ptp);
#' 6. Equal.training.sensitivity.and.specificity (etss);
#' 7. Maximum.training.sensitivity.plus.specificity (mtss);
#' 8. Balance.training.omission.predicted.area.and.threshold.value (bto);
#' 9. Equate.entropy.of.thresholded.and.original.distributions (eetd).
#' @seealso \code{\link{thrshld_b}}
#' @return Stack or brick of thresholded predictions
#' @examples
#' \dontrun{
#' mods.thrshld <- thrshld(mcmp=mxnt.mdls.preds, thrshld.i = 4:6, pred.args, path.mdls)
#' plot(mods.thrshld[[1]][[2]]) # continuous
#' plot(mods.thrshld[[2]][[2]]) # binary
#' }
#' @export
thrshld <- function(mcmp, thrshld.i = 4:6, sp.nm="species", numCores=1) { # path.mdls = NULL,
  ###- get necessary objects
  mxnt.mdls <- mcmp[["mxnt.mdls"]]
  pred.args <- mcmp$pred.args
  mod.nms <- mcmp[["selected.mdls"]]$sel.cri # gsub("Mod.", "", names(mcmp[["mxnt.preds"]][[1]]))
  names(mod.nms) <- apply(mcmp[["selected.mdls"]][,c("rm", "features")], 1,
                          function(x){
                            paste0("Mod_", paste(x, collapse = "_"))
                          })
  # mod.nms  <- names(mcmp$mxnt.mdls)
  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

  ###- define and create output path
  path.mdls <- paste("3_out.MaxEnt", paste0("Mdls.",sp.nm), sep = "/")
  thrshld.path <- paste(path.mdls, outpt, "Mdls.thrshld", sep='/')
  if(dir.exists(thrshld.path)==FALSE) {dir.create(thrshld.path, recursive = T)}

  ###- get threshold name from maxent output
  # thrshld.nms <- c("fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd")[thrshld.i]
  thrshld.nms <- tnm[thrshld.i]
  thrshld.crit <- rownames(mxnt.mdls[[1]]@results)[grepl(outpt, rownames(mxnt.mdls[[1]]@results), ignore.case = T) & # TODO use "outpt" variable
                                                     grepl("threshold", rownames(mxnt.mdls[[1]]@results))][thrshld.i]

  ###- extract threshold values from selected models and criteria
  thrshld.crit.v <- as.data.frame(matrix(data=sapply(thrshld.crit,
                                                     function(y){
                                                       sapply(mxnt.mdls[names(mod.nms)],
                                                              function(x){
                                                                x@results[rownames(x@results) == y]
                                                              }) }),
                                         ncol = length(thrshld.i)))
  colnames(thrshld.crit.v) <- thrshld.nms
  rownames(thrshld.crit.v) <- mod.nms
  # thrshld.crit.v

  ####- Ensemble models
  ###- get threshold values of individual models
  mod.sel.crit <- mcmp$mSel #names(pred.r)
  if(sum(grepl("AvgAIC", mod.sel.crit))>0) {
    wv.aic <- mcmp[["selected.mdls"]][grep("AIC_", mcmp[["selected.mdls"]]$sel.cri),"w.AIC"]
  }
  if(sum(grepl("WAAUC", mod.sel.crit))>0) {
    wv.wa <- mcmp[["selected.mdls"]][grep("WAAUC_", mcmp[["selected.mdls"]]$sel.cri),"avg.test.AUC"]
  }
  if(sum(grepl("EBPM", mod.sel.crit))>0) {
    wv.bp <- rep(1, length(grep("EBPM_", mcmp[["selected.mdls"]]$sel.cri)))
  }
  if(sum(grepl("ESORIC", mod.sel.crit))>0) {
    wv.es <- rep(1, length(grep("ESORIC_", mcmp[["selected.mdls"]]$sel.cri)))
  }

  ###- compute threshold values of ensemble models using individual models
  thrshld.mod.crt <- rbind(
    if(sum(grepl("AIC_", mod.nms))>1){
      matrix(apply(data.frame(thrshld.crit.v[grep("AIC_", mcmp[["selected.mdls"]]$sel.cri),]), 2, function(x, wv) {
        stats::weighted.mean(x, wv)
      }, wv.aic), nrow = 1, dimnames = list("AvgAIC", thrshld.nms) )
    } , # else {thrshld.crit <- thrshld.crit.v}
    if(sum(grepl("WAAUC_", mod.nms))>0){
      matrix(apply(data.frame(thrshld.crit.v[grep("WAAUC_", mcmp[["selected.mdls"]]$sel.cri),]), 2, function(x, wv) {
        stats::weighted.mean(x, wv)
      }, wv.wa), nrow = 1, dimnames = list("WAAUC", thrshld.nms) )
    } ,
    if(sum(grepl("EBPM_", mod.nms))>0){
      matrix(apply(data.frame(thrshld.crit.v[grep("EBPM_", mcmp[["selected.mdls"]]$sel.cri),]), 2, function(x, wv) {
        stats::weighted.mean(x, wv)
      }, wv.bp), nrow = 1, dimnames = list("EBPM", thrshld.nms) )
    } ,
    if(sum(grepl("ESORIC_", mod.nms))>0){
      matrix(apply(data.frame(thrshld.crit.v[grep("ESORIC_", mcmp[["selected.mdls"]]$sel.cri),]), 2, function(x, wv) {
        stats::weighted.mean(x, wv)
      }, wv.es), nrow = 1, dimnames = list("ESORIC", thrshld.nms) )
    } ,
    thrshld.crit.v)

  ###- remove individual models used to build ensemble
  s.nms <- c("LowAIC", "ORmtp", "OR10", "AUCmtp", "AUC10", "^AvgAIC$", "^EBPM$", "^WAAUC$", "^ESORIC$")
  thrshld.mod.crt <- subset(thrshld.mod.crt, grepl(paste0(s.nms, collapse = "|"), rownames(thrshld.mod.crt)))

  ###- get projections to apply thresholds
  # TODO - choose between consensus, mxnt.preds, or both
  scn.nms <- c(names(mcmp$mxnt.preds), names(mcmp$scn.consensus))

  ###- workhorse function
  fthr <- function(j, scn.nms, mcmp, thrshld.i, sp.nm,
                   # thrshld.crit,
                   mod.sel.crit, thrshld.path, thrshld.mod.crt){
    scn.nm <- scn.nms[j]
    if(j <= length(mcmp$mxnt.preds)){
      pred.r <- mcmp$mxnt.preds[[scn.nm]] # mcmp[[match(pred.nm, names(mcmp))]] # , fixed=TRUE # [pred.i]
    } else {
      pred.r <- mcmp$scn.consensus[[scn.nm]]
    }

    thrshld.nms <- colnames(thrshld.mod.crt)
    brick.nms.t <- paste0("mxnt.pred.", scn.nm, ".", thrshld.nms)
    brick.nms.t.b <- paste0("mxnt.pred.", scn.nm, ".", thrshld.nms, ".b")

    mt.lst <- vector("list", length = length(thrshld.nms))
    names(mt.lst) <- thrshld.nms
    mt <- list(continuous=mt.lst, binary=mt.lst)

    for(t in base::seq_along(thrshld.nms)){

      mod.sel.crit.t <- paste(paste0(mod.sel.crit, ".", scn.nm), thrshld.nms[t], sep=".")
      mod.sel.crit.b <- paste(paste0(mod.sel.crit, ".", scn.nm, ".b"), thrshld.nms[t], sep=".")

      pred.t <- pred.r
      pred.t <- raster::stack(lapply(mod.sel.crit, # seq_along(mod.sel.crit)
                                     function(m, pred.t, thrshld.mod.crt, t) {
                                       mp <- grep(m, names(pred.t))
                                       mt <- grep(m, rownames(thrshld.mod.crt))
                                       pred.t[[mp]][pred.t[[mp]] < thrshld.mod.crt[mt,t]] <- 0
                                       return(pred.t[[mp]])
                                     }, pred.t, thrshld.mod.crt, t))

      names(pred.t) <- mod.sel.crit.t

      mt$continuous[[thrshld.nms[t]]] <- raster::writeRaster(x = pred.t,
                                                             filename = paste(thrshld.path, paste0("mxnt.pred", gsub(".mxnt.pred","", paste0(".",scn.nm)), ".", thrshld.nms[t], ".grd"), sep='/'),
                                                             format = "raster", overwrite = T) #)

      # create presence only raster
      pred.t <- raster::stack(lapply(mod.sel.crit, # seq_along(mod.sel.crit)
                                     function(m) {
                                       mp <- grep(m, names(pred.t))
                                       mt <- grep(m, rownames(thrshld.mod.crt))
                                       pred.t[[mp]][pred.t[[mp]] >= thrshld.mod.crt[mt,t]] <- 1
                                       return(pred.t[[mp]])
      }))
      mt$binary[[thrshld.nms[t]]] <- raster::writeRaster(x = pred.t,
                                                         filename = paste(thrshld.path, paste0("mxnt.pred", gsub(".mxnt.pred","", paste0(".",scn.nm)), ".", thrshld.nms[t], ".b", ".grd"), sep='/'),
                                                         format = "raster", overwrite = T)
    }
    return(mt)
  }

  ###- raster computation
  {
    if(numCores>1){
      cl <- parallel::makeCluster(numCores)
      parallel::clusterExport(cl, list("thrshld")) # , "scn.nms"

      mods.thrshld.spi <- parallel::clusterApply(cl, base::seq_along(scn.nms), # scn.nms
                                                 fthr,
                                                 scn.nms, mcmp, thrshld.i, sp.nm,
                                                 # thrshld.crit,
                                                 mod.sel.crit, thrshld.path, thrshld.mod.crt)
      parallel::stopCluster(cl)

    } else {
      mods.thrshld.spi <- lapply(base::seq_along(scn.nms), # scn.nms
                                 fthr,
                                 scn.nms, mcmp, thrshld.i, sp.nm,
                                 # thrshld.crit,
                                 mod.sel.crit, thrshld.path, thrshld.mod.crt)
    }
    names(mods.thrshld.spi) <- scn.nms
  }
  return(mods.thrshld.spi)
  # return(mt)
}


#### 4.8.5 threshold for past and future pred
#' Apply threshold for MaxEnt projections for multiple species
#'
#' This function will apply the selected threshold criterias to MaxEnt model projection(s) of a 'mcmp.l' object
#' and save on the folder "3_out.MaxEnt/Mdls.[species name]/Mdls.thrshld". For each projection (species and climatic
#' scenario), two layers will be generated, one with suitability above the threshold value and another with presence/absence only.
#'
#' @param mcmp.l Object returned by \code{\link{proj_mdl_b}}, containing a list of calibrated models
#' and model projections for each species.
#' @inheritParams thrshld
#' @inheritParams calib_mdl
#' @seealso \code{\link{thrshld}}
#' @return List of stack or brick of thresholded predictions
#' @examples
#' \dontrun{
#' mods.thrshld.lst <- thrshld_b(mcmp.l=mxnt.mdls.preds.cf)
#' }
#' @export
thrshld_b <- function(mcmp.l, thrshld.i = 4:6, numCores = 1) {
  # thrshld for each species
  mods.thrshld <- vector("list", length(mcmp.l))
  names(mods.thrshld) <- names(mcmp.l)
  for(i in names(mcmp.l)){ # species i
    scn.nms <- c(names(mcmp.l[[i]]$mxnt.preds), names(mcmp.l[[i]]$scn.consensus)) # # scn.nms <- names(mcmp.l[[i]]$mxnt.preds)
    mods.thrshld.spi <- thrshld(mcmp.l[[i]], thrshld.i, sp.nm = i, numCores) # sp.nm = names(mcmp.l)[i]
    names(mods.thrshld.spi) <- scn.nms
    mods.thrshld[[i]] <- mods.thrshld.spi
  }
  return(mods.thrshld)
}


#' Threshold names to sub
#' @keywords internal
tnm <- c("fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd")
#' Threshold names to sub
#' @keywords internal
tr <- c("FCV1", "FCV5", "FCV10", "LPT (mtp)", "10P (x10ptp)", "ETSS (etss)", "MTSS", "BTO", "EETD")



