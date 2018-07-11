#### 4.3.3 aplicar threshold
# TODO examples
# name of arg "mxnt.mdls.preds.sp[...]" shortened to "mcmp.spi"

#' Apply threshold for MaxEnt projections of a species
#'
#' This function will apply the selected threshold criterias to MaxEnt model projection(s) of a 'mcmp' object
#' and save on the folder "3_out.MaxEnt/Mdls.[species name]/Mdls.thrshld". For each projection (species and climatic
#' scenario), two layers will be generated, one with suitability above the threshold value and another with presence/absence only.
#'
#' @param mcmp.spi Species "i" of a object returned by "mxnt.p.batch.mscn", containing a list of
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
#' @seealso \code{\link{f.thr.batch}}
#' @return Stack or brick of thresholded predictions
#' @examples
#' mods.thrshld <- f.thr(mcmp.spi=mxnt.mdls.preds, thrshld.i = 4:6, pred.args, path.mdls)
#' plot(mods.thrshld[[1]][[2]]) # continuous
#' plot(mods.thrshld[[2]][[2]]) # binary
#' @keywords internal
# #' @export
f.thr <- function(mcmp.spi, scn.nm = "", thrshld.i = 4:6, path.mdls = NULL) {
  if(is.null(path.mdls)){
    path.mdls <- paste("3_out.MaxEnt", "Mdls.sp", sep = "/")
  }
  if(dir.exists(path.mdls)==FALSE) {dir.create(path.mdls)}
  mxnt.mdls <- mcmp.spi[["mxnt.mdls"]]

  #### TODO use the "slot" to find and loop through predictions
  ##  check here
  pred.r <- mcmp.spi$mxnt.preds[[scn.nm]] # mcmp.spi[[match(pred.nm, names(mcmp.spi))]] # , fixed=TRUE # [pred.i]
  pred.args <- mcmp.spi$pred.args
  mod.sel.crit <- names(pred.r)
  mod.nms  <- mcmp.spi[["selected.mdls"]]$sel.cri

  if(sum(grepl("AvgAICc", mod.sel.crit))>0) {
    wv <- mcmp.spi[["selected.mdls"]][order(mcmp.spi[["selected.mdls"]]$delta.AICc[grep("AICc_", mcmp.spi[["selected.mdls"]]$sel.cri)]),"w.AIC"]
  }

  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

  thrshld.path <- paste(path.mdls, outpt, "Mdls.thrshld", sep='/')

  if(dir.exists(thrshld.path)==FALSE) {dir.create(thrshld.path)}

  thrshld.nms <- c("fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd")[thrshld.i]
  thrshld.crit <- rownames(mxnt.mdls[[1]]@results)[grepl("Cloglog", rownames(mxnt.mdls[[1]]@results)) & # TODO use "outpt" variable
                                                     grepl("threshold", rownames(mxnt.mdls[[1]]@results))][thrshld.i]

  ## extract values of threshold from each model and criteria # sapply(mxnt.mdls[1:length(args)]
  thrshld.crit.v <- as.data.frame(matrix(data=sapply(thrshld.crit, function(y) sapply(mxnt.mdls[1:length(mxnt.mdls)], function(x) x@results[rownames(x@results) == y]) ),
                                         ncol = length(thrshld.i)))
  colnames(thrshld.crit.v) <- thrshld.nms
  rownames(thrshld.crit.v) <- mod.nms
  thrshld.crit.v



  thrshld.mod.crt <- rbind(if(sum(grepl("AICc_", mod.nms))>1){
    matrix(apply(data.frame(thrshld.crit.v[grep("AICc_", mcmp.spi[["selected.mdls"]]$sel.cri),]), 2, function(x, wv) {
      stats::weighted.mean(x, wv)
    }, wv), nrow = 1, dimnames = list("AvgAICc", thrshld.nms) )
  } else {thrshld.crit <- thrshld.crit.v},
  thrshld.crit.v)
  thrshld.mod.crt <- subset(thrshld.mod.crt, grepl(paste0(mcmp.spi$mSel, collapse = "|"), rownames(thrshld.mod.crt)))

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

# mods.thrshld <- f.thr(mcmp.spi, thrshld.i = 4:6, pred.args, path.mdls)
# plot(mods.thrshld[[1]][[2]]) # continuous
# plot(mods.thrshld[[2]][[2]]) # binary

#### 4.8.5 threshold for past and future pred
#' Apply threshold for MaxEnt projections across all species
#'
#' This function will apply the selected threshold criterias to MaxEnt model projection(s) of a 'mcmp.l' object
#' and save on the folder "3_out.MaxEnt/Mdls.[species name]/Mdls.thrshld". For each projection (species and climatic
#' scenario), two layers will be generated, one with suitability above the threshold value and another with presence/absence only.
#'
#' @param mcmp.l Object returned by "mxnt.p.batch.mscn", containing a list of calibrated models
#' and model projections for each species.
#' @inheritParams f.thr
#' @inheritParams mxnt.c
#' @seealso \code{\link{f.thr}}
#' @return List of stack or brick of thresholded predictions
#' @examples
#' mods.thrshld.lst <- f.thr.batch(mcmp.l=mxnt.mdls.preds.cf)
#' @export
f.thr.batch <- function(mcmp.l, thrshld.i = 4:6, numCores = 1) {
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
    mcmp.spi <- mcmp.l[[i]]
    scn.nms <- names(mcmp.l[[i]]$mxnt.preds)

    if(numCores>1){
      cl <- parallel::makeCluster(numCores)
      parallel::clusterExport(cl, list("f.thr"))

      mods.thrshld.spi <- parallel::clusterApply(cl, base::seq_along(scn.nms),
                                                 function(j, mcmp.spi, scn.nms, thrshld.i, path.mdls.i){

                                                   resu <- f.thr(mcmp.spi, scn.nm = scn.nms[j], thrshld.i, path.mdls = path.mdls.i)
                                                   return(resu)

                                                 }, mcmp.spi, scn.nms, thrshld.i, path.mdls.i)
      parallel::stopCluster(cl)

    }else{
      mods.thrshld.spi <- lapply(base::seq_along(scn.nms),
                                 function(j, mcmp.spi, scn.nms, thrshld.i, path.mdls.i){

                                   resu <- f.thr(mcmp.spi, scn.nm = scn.nms[j], thrshld.i, path.mdls = path.mdls.i)
                                   return(resu)

                                 }, mcmp.spi, scn.nms, thrshld.i, path.mdls.i)
    }
    names(mods.thrshld.spi) <- scn.nms
    mods.thrshld[[i]] <- append(mods.thrshld[[i]], mods.thrshld.spi)
  }
  return(mods.thrshld)
}




