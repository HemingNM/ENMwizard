#### 4.3.3 aplicar threshold
# TODO examples
# name of arg "mxnt.mdls.preds.sp[...]" shortened to "mcmp"

#' Apply threshold for MaxEnt projections of a species
#'
#' This function will apply the selected threshold criterias to MaxEnt model projection(s) of a 'mcmp' object
#' and save on the folder "3_out.MaxEnt/Mdls.[species name]/Mdls.thrshld". For each projection (species and climatic
#' scenario), two layers will be generated, one with suitability above the threshold value and another with presence/absence only.
#'
#' @inheritParams mxntCalib
#' @param mcmp Species "i" of a object returned by "mxntProjB", containing a list of
#' calibrated models and model projections for each species
#' @param scn.nm Name of climatic scenario to be looked for
#' @param path.mdls Path where thresholded rasters will be saved
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
#' @seealso \code{\link{thrB}}
#' @return Stack or brick of thresholded predictions
#' @examples
#' mods.thrshld <- thr(mcmp=mxnt.mdls.preds, thrshld.i = 4:6, pred.args, path.mdls)
#' plot(mods.thrshld[[1]][[2]]) # continuous
#' plot(mods.thrshld[[2]][[2]]) # binary
#' @export
thr <- function(mcmp, scn.nm = "", thrshld.i = 4:6, path.mdls = NULL, sp.nm="species") {
  if(is.null(path.mdls)){
    path.mdls <- paste("3_out.MaxEnt", paste0("Mdls.",sp.nm), sep = "/")
  }
  if(dir.exists(path.mdls)==FALSE) {dir.create(path.mdls)}
  mxnt.mdls <- mcmp[["mxnt.mdls"]]

  #### TODO use the "slot" to find and loop through predictions
  ##  check here
  pred.r <- mcmp$mxnt.preds[[scn.nm]] # mcmp[[match(pred.nm, names(mcmp))]] # , fixed=TRUE # [pred.i]
  mod.sel.crit <- names(pred.r)
  pred.args <- mcmp$pred.args
  mod.nms  <- mcmp[["selected.mdls"]]$sel.cri

  if(sum(grepl("AvgAIC", mod.sel.crit))>0) {
    wv.aic <- mcmp[["selected.mdls"]][grep("AIC_", mcmp[["selected.mdls"]]$sel.cri),"w.AIC"]
  }
  if(sum(grepl("WAAUC", mod.sel.crit))>0) {
    wv.wa <- mcmp[["selected.mdls"]][grep("WAAUC_", mcmp[["selected.mdls"]]$sel.cri),"avg.test.AUC"]
  }
  if(sum(grepl("EBPM", mod.sel.crit))>0) {
    wv.bp <- rep(1, length(grep("EBPM", mcmp[["selected.mdls"]]$sel.cri)))
  }

  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

  thrshld.path <- paste(path.mdls, outpt, "Mdls.thrshld", sep='/')

  if(dir.exists(thrshld.path)==FALSE) {dir.create(thrshld.path)}

  # thrshld.nms <- c("fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd")[thrshld.i]
  thrshld.nms <- tnm[thrshld.i]
  thrshld.crit <- rownames(mxnt.mdls[[1]]@results)[grepl("Cloglog", rownames(mxnt.mdls[[1]]@results)) & # TODO use "outpt" variable
                                                     grepl("threshold", rownames(mxnt.mdls[[1]]@results))][thrshld.i]

  ## extract values of threshold from each model and criteria # sapply(mxnt.mdls[1:length(args)]
  thrshld.crit.v <- as.data.frame(matrix(data=sapply(thrshld.crit, function(y) sapply(mxnt.mdls[1:length(mxnt.mdls)], function(x) x@results[rownames(x@results) == y]) ),
                                         ncol = length(thrshld.i)))
  colnames(thrshld.crit.v) <- thrshld.nms
  rownames(thrshld.crit.v) <- mod.nms
  thrshld.crit.v



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
    thrshld.crit.v)

  ### subset thrshld.mod.crt
  s.nms <- c("LowAIC", "ORmtp", "OR10", "AUCmtp", "AUC10", "^AvgAIC$", "^EBPM$", "^WAAUC$")
  thrshld.mod.crt <- subset(thrshld.mod.crt, grepl(paste0(s.nms, collapse = "|"), rownames(thrshld.mod.crt)))

  brick.nms.t <- paste0("mxnt.pred.", scn.nm, ".", thrshld.nms)
  brick.nms.t.b <- paste0("mxnt.pred.", scn.nm, ".", thrshld.nms, ".b")


  mt.lst <- vector("list", length = length(thrshld.nms))
  names(mt.lst) <- thrshld.nms
  mt <- list(continuous=mt.lst, binary=mt.lst)

  lapply(base::seq_along(thrshld.crit),
         function(t, mod.sel.crit, scn.nm, thrshld.nms, thrshld.path, thrshld.mod.crt,
                  pred.r, brick.nms.t, brick.nms.t.b, mt){


           mod.sel.crit.t <- paste(paste0(mod.sel.crit, ".", scn.nm), thrshld.nms[t], sep=".")
           mod.sel.crit.b <- paste(paste0(mod.sel.crit, ".", scn.nm, ".b"), thrshld.nms[t], sep=".")

           pred.t <- pred.r
           pred.t <- raster::stack(lapply(seq_along(mod.sel.crit), function(m, pred.t, thrshld.mod.crt, t) {
             pred.t[[m]][pred.t[[m]] < thrshld.mod.crt[m,t]] <- 0
             return(pred.t[[m]])
           }, pred.t, thrshld.mod.crt, t))

           names(pred.t) <- mod.sel.crit.t

           mt$continuous[[thrshld.nms[t]]] <<- raster::writeRaster(x = pred.t,
                                                                   filename = paste(thrshld.path, paste0("mxnt.pred", gsub(".mxnt.pred","", paste0(".",scn.nm)), ".", thrshld.nms[t], ".grd"), sep='/'),
                                                                   format = "raster", overwrite = T) #)

           # create presence only raster
           pred.t <- raster::stack(lapply(seq_along(mod.sel.crit), function(m) {
             pred.t[[m]][pred.t[[m]] >= thrshld.mod.crt[m,t]] <- 1
             return(pred.t[[m]])
           }))
           mt$binary[[thrshld.nms[t]]] <<- raster::writeRaster(x = pred.t,
                                                               filename = paste(thrshld.path, paste0("mxnt.pred", gsub(".mxnt.pred","", paste0(".",scn.nm)), ".", thrshld.nms[t], ".b", ".grd"), sep='/'),
                                                               format = "raster", overwrite = T)

           return(thrshld.nms[t])

         }, mod.sel.crit, scn.nm, thrshld.nms, thrshld.path, thrshld.mod.crt,
         pred.r, brick.nms.t, brick.nms.t.b, mt)

  return(mt)
}

# mods.thrshld <- thr(mcmp, thrshld.i = 4:6, pred.args, path.mdls)
# plot(mods.thrshld[[1]][[2]]) # continuous
# plot(mods.thrshld[[2]][[2]]) # binary

#### 4.8.5 threshold for past and future pred
#' Apply threshold for MaxEnt projections for multiple species
#'
#' This function will apply the selected threshold criterias to MaxEnt model projection(s) of a 'mcmp.l' object
#' and save on the folder "3_out.MaxEnt/Mdls.[species name]/Mdls.thrshld". For each projection (species and climatic
#' scenario), two layers will be generated, one with suitability above the threshold value and another with presence/absence only.
#'
#' @param mcmp.l Object returned by "mxntProjB", containing a list of calibrated models
#' and model projections for each species.
#' @inheritParams thr
#' @inheritParams mxntCalib
#' @seealso \code{\link{thr}}
#' @return List of stack or brick of thresholded predictions
#' @examples
#' mods.thrshld.lst <- thrB(mcmp.l=mxnt.mdls.preds.cf)
#' @export
thrB <- function(mcmp.l, thrshld.i = 4:6, numCores = 1) {
  path.res <- "3_out.MaxEnt"
  if(dir.exists(path.res)==FALSE) dir.create(path.res)
  path.sp.m <- paste0("Mdls.", names(mcmp.l))
  path.mdls <- paste(path.res, path.sp.m, sep="/")

  # thrshld for each species
  mods.thrshld <- vector("list", length(mcmp.l))
  names(mods.thrshld) <- names(mcmp.l)

  for(i in base::seq_along(mcmp.l)){ # species i
    path.mdls.i <- path.mdls[i]
    if(dir.exists(path.mdls.i)==FALSE) dir.create(path.mdls.i)
    mcmp <- mcmp.l[[i]]
    scn.nms <- names(mcmp.l[[i]]$mxnt.preds)

    if(numCores>1){
      cl <- parallel::makeCluster(numCores)
      parallel::clusterExport(cl, list("thr"))

      mods.thrshld.spi <- parallel::clusterApply(cl, base::seq_along(scn.nms),
                                                 function(j, mcmp, scn.nms, thrshld.i, path.mdls.i){

                                                   resu <- thr(mcmp, scn.nm = scn.nms[j], thrshld.i, path.mdls = path.mdls.i)
                                                   return(resu)

                                                 }, mcmp, scn.nms, thrshld.i, path.mdls.i)
      parallel::stopCluster(cl)

    }else{
      mods.thrshld.spi <- lapply(base::seq_along(scn.nms),
                                 function(j, mcmp, scn.nms, thrshld.i, path.mdls.i){

                                   resu <- thr(mcmp, scn.nm = scn.nms[j], thrshld.i, path.mdls = path.mdls.i)
                                   return(resu)

                                 }, mcmp, scn.nms, thrshld.i, path.mdls.i)
    }
    names(mods.thrshld.spi) <- scn.nms
    mods.thrshld[[i]] <- append(mods.thrshld[[i]], mods.thrshld.spi)
  }
  return(mods.thrshld)
}


### threshold names to sub
#' @keywords internal
tnm <- c("fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd")
#' @keywords internal
tr <- c("FCV1", "FCV5", "FCV10", "LPT (mtp)", "10P (x10ptp)", "ETSS (etss)", "MTSS", "BTO", "EETD")



