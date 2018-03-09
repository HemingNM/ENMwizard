#### 4.3.3 aplicar threshold
# TODO examples
# name of arg "mxnt.mdls.preds.sp[...]" shortened to "mcmp.spi"
#' Apply threshold for a prediction
#'
#' This function will apply the selected threshold criterias to MaxEnt model projection(s) of a 'mcmp' object
#' and save on the folder "4_ENMeval.results/Mdls.[species name]/Mdls.thrshld". For each projection (species and climatic
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
#' mods.thrshld <- f.thr(mxnt.mdls.preds, thrshld.i = 4:6, pred.args, path.mdls)
#' plot(mods.thrshld[[1]][[2]]) # continuous
#' plot(mods.thrshld[[2]][[2]]) # binary
#' @export
f.thr <- function(mcmp.spi, scn.nm = "", thrshld.i = 4:6, path.mdls = NULL) {
  if(is.null(path.mdls)){
    path.mdls <- paste("4_ENMeval.results", "Mdls.sp", sep = "/")
  }
  if(dir.exists(path.mdls)==F) dir.create(path.mdls)
  # args <- 1:length(mcmp.spi[["mxnt.mdls"]])
  mxnt.mdls <- mcmp.spi[["mxnt.mdls"]]

  # pred.nm <- ifelse(pred.nm == "", pred.nm, paste0(".", pred.nm))
  # m.pred.n <- ifelse(pred.nm == "", "mxnt.pred", paste0("mxnt.pred", pred.nm))
  # m.pred.n <- ifelse(pred.nm == "", "mxnt.pred", paste0(pred.nm))
  # print(paste("pred.nm is ", pred.nm))
  # print(paste("m.pred.n is ", m.pred.n))

  #### TODO use the "slot" to find and loop through predictions
  ##  check here
  pred.r <- mcmp.spi$mxnt.preds[[scn.nm]] # mcmp.spi[[match(pred.nm, names(mcmp.spi))]] # , fixed=T # [pred.i]
  pred.args <- mcmp.spi$pred.args
  mod.sel.crit <- names(pred.r)

  if(sum(grepl("AvgAICc", mod.sel.crit))>0) {
    # args.aicc <- 1:(length(args)-4)
    wv <- mcmp.spi[["selected.mdls"]][order(mcmp.spi[["selected.mdls"]]$delta.AICc[grep("Mod.AICc", mcmp.spi[["selected.mdls"]]$sel.cri)]),"w.AIC"]
  } # else { wv <- 1 } # rep(1, 1+length(grep("Mod.AICc", mcmp.spi[["selected.mdls"]]$sel.cri, invert = T))) } # else {wv <- rep(0, max(grep("Mod.AICc", mcmp.spi[[1]]$sel.cri))) }


  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

  thrshld.path <- paste(path.mdls, outpt, "Mdls.thrshld", sep='/')

  if(dir.exists(thrshld.path)==F) dir.create(thrshld.path)

  thrshld.nms <- c("fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd")[thrshld.i]
  thrshld.crit <- rownames(mxnt.mdls[[1]]@results)[grepl("Cloglog", rownames(mxnt.mdls[[1]]@results)) & # TODO use "outpt" variable
                                                     grepl("threshold", rownames(mxnt.mdls[[1]]@results))][thrshld.i]

  ## extract values of threshold from each model and criteria # sapply(mxnt.mdls[1:length(args)]
  thrshld.crit.v <- as.data.frame(matrix(data=sapply(thrshld.crit, function(y) sapply(mxnt.mdls[1:length(mxnt.mdls)], function(x) x@results[rownames(x@results) == y]) ),
                                         ncol = length(thrshld.i)))

  # thrshld.mod.crt <- thrshld.crit.v[1]
  # grep("LowAICc", names(pred.r))

  thrshld.mod.crt <- data.frame(rbind(
    if(sum(grepl("AvgAICc", names(pred.r)))>0){ # # 1:length(args.aicc)
      apply(data.frame(thrshld.crit.v[grep("Mod.AICc", mcmp.spi[["selected.mdls"]]$sel.cri),]), 2, function(x, wv) {
        # if(length(x) != length(w) & length(w) == 1) {w <- rep(1, length(x))}
        return(stats::weighted.mean(x, wv))
        }, wv) ### check if is raster
    } else {NA}, # compute avg.thrshld from each criteria weighted by model importance (AICc W)
    if(sum(grepl("LowAICc", names(pred.r)))>0){
      thrshld.crit.v[grep("Mod.AICc_1$", mcmp.spi[["selected.mdls"]]$sel.cri),]
    } else {NA},
    if(sum(grepl("Mean.ORmin", names(pred.r)))>0){
      thrshld.crit.v[grep("Mean.ORmin", mcmp.spi[["selected.mdls"]]$sel.cri),]
    }else {NA},
    if(sum(grepl("Mean.OR10", names(pred.r)))>0){
      thrshld.crit.v[grep("Mean.OR10", mcmp.spi[["selected.mdls"]]$sel.cri),]
    } else {NA},
    if(sum(grepl("Mean.AUCmin", names(pred.r)))>0){
      thrshld.crit.v[grep("Mean.AUCmin", mcmp.spi[["selected.mdls"]]$sel.cri),]
    } else {NA},
    if(sum(grepl("Mean.AUC10", names(pred.r)))>0){
      thrshld.crit.v[grep("Mean.AUC10", mcmp.spi[["selected.mdls"]]$sel.cri),]
    } else {NA} ) )

  # thrshld.mod.crt <- cbind(thrshld.mod.crt, thrshld.mod.crt)
  thrshld.mod.crt <- as.data.frame(thrshld.mod.crt[!is.na(thrshld.mod.crt[,1]),])
  row.names(thrshld.mod.crt) <- mod.sel.crit
  colnames(thrshld.mod.crt) <- paste0("thrshld.", thrshld.nms)

  brick.nms.t <- paste0("mxnt.pred.", scn.nm, ".", thrshld.nms)
  brick.nms.t.b <- paste0("mxnt.pred.", scn.nm, ".", thrshld.nms, ".b")


  mt.lst <- vector("list", length = length(thrshld.nms))
  names(mt.lst) <- thrshld.nms
  mt <- list(continuous=mt.lst, binary=mt.lst)

  # mt2 <-
  lapply(base::seq_along(thrshld.crit),
         function(t, mod.sel.crit, scn.nm, thrshld.nms, thrshld.path, thrshld.mod.crt,
                  pred.r, brick.nms.t, brick.nms.t.b, mt){


           # for(t in base::seq_along(thrshld.crit)){
           mod.sel.crit.t <- paste(paste0(mod.sel.crit, ".", scn.nm), thrshld.nms[t], sep=".")
           mod.sel.crit.b <- paste(paste0(mod.sel.crit, ".", scn.nm, ".b"), thrshld.nms[t], sep=".")

           pred.t <- pred.r

           pred.t <- raster::stack(lapply(seq_along(mod.sel.crit), function(m) {
             pred.t[[m]][pred.t[[m]] < thrshld.mod.crt[m,t]] <- 0
             return(pred.t[[m]])
           }))
           # for(m in base::seq_along(mod.sel.crit)){
           #   pred.t[[m]][pred.t[[m]] < thrshld.mod.crt[m,t]] <- 0
           # }
           names(pred.t) <- mod.sel.crit.t
           # assign(brick.nms.t[t],
           mt$continuous[[thrshld.nms[t]]] <<- #continuous <-
             raster::writeRaster(x = pred.t,
                                 filename = paste(thrshld.path, paste0("mxnt.pred", gsub(".mxnt.pred","", paste0(".",scn.nm)), ".", thrshld.nms[t], ".grd"), sep='/'),
                                 format = "raster", overwrite = T) #)

           # create presence only raster
           pred.t <- raster::stack(lapply(seq_along(mod.sel.crit), function(m) {
             pred.t[[m]][pred.t[[m]] >= thrshld.mod.crt[m,t]] <- 1
             return(pred.t[[m]])
           }))
           # for(m in base::seq_along(mod.sel.crit)){
           #   pred.t[[m]][pred.t[[m]] >= thrshld.mod.crt[m,t]] <- 1
           # }
           # assign(brick.nms.t.b[t],
           mt$binary[[thrshld.nms[t]]] <<- #binary <-
             raster::writeRaster(x = pred.t,
                                 filename = paste(thrshld.path, paste0("mxnt.pred", gsub(".mxnt.pred","", paste0(".",scn.nm)), ".", thrshld.nms[t], ".b", ".grd"), sep='/'),
                                 format = "raster", overwrite = T) #)
           # }
           # return(list(continuous=continuous, binary=binary))
           # return(mt)
           return(thrshld.nms[t])

         }, mod.sel.crit, scn.nm, thrshld.nms, thrshld.path, thrshld.mod.crt,
         pred.r, brick.nms.t, brick.nms.t.b, mt)

  # mods.t <- lapply(brick.nms.t, function(x) get(x))
  # mods.t.b <- lapply(brick.nms.t.b, function(x) get(x))
  # names(mods.t) <- thrshld.nms
  # names(mods.t.b) <- thrshld.nms
  # return(list(continuous=mods.t, binary=mods.t.b))

  return(mt)
}

# mods.thrshld <- f.thr(mcmp.spi, thrshld.i = 4:6, pred.args, path.mdls)
# plot(mods.thrshld[[1]][[2]]) # continuous
# plot(mods.thrshld[[2]][[2]]) # binary

#### 4.8.5 threshold for past and future pred
#' Apply threshold for all predictions
#'
#' This function will apply the selected threshold criterias to MaxEnt model projection(s) of a 'mcmp.l' object
#' and save on the folder "4_ENMeval.results/Mdls.[species name]/Mdls.thrshld". For each projection (species and climatic
#' scenario), two layers will be generated, one with suitability above the threshold value and another with presence/absence only.
#'
#' @param mcmp.l Object returned by "mxnt.p.batch.mscn", containing a list of calibrated models
#' and model projections for each species.
#' @inheritParams f.thr
#' @inheritParams mxnt.c
#' @seealso \code{\link{f.thr}}
#' @return List of stack or brick of thresholded predictions
#' @examples
#' mods.thrshld.lst <- f.thr.batch(mxnt.mdls.preds.pf)
#' @export
f.thr.batch <- function(mcmp.l, thrshld.i = 4:6, numCores = 1) {
  path.res <- "4_ENMeval.results"
  if(dir.exists(path.res)==F) dir.create(path.res)
  path.sp.m <- paste0("Mdls.", names(mcmp.l))
  path.mdls <- paste(path.res, path.sp.m, sep="/")

  # thrshld for each species
  mods.thrshld <- vector("list", length(mcmp.l))
  names(mods.thrshld) <- names(mcmp.l)
  # pred.nm <- ifelse(pred.nm == "", "mxnt.pred", pred.nm)

  for(i in base::seq_along(mcmp.l)){ # species i
    path.mdls.i <- path.mdls[i]
    if(dir.exists(path.mdls.i)==F) dir.create(path.mdls.i)
    #scn.ind <- grep(n.pred.nm, names(mcmp.l[[i]]))
    mcmp.spi <- mcmp.l[[i]]
    # scn.ind <- grep(pred.nm, names(mcmp.l[[i]]))
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

      # mods.thrshld.spi <- stats::setNames(vector("list", length(scn.nms)), scn.nms)
      # for(j in base::seq_along(scn.nms)){ # climatic scenario
      #   mods.thrshld.spi[[j]] <- f.thr(mcmp.spi, scn.nm = scn.nms[j], thrshld.i, path.mdls = path.mdls.i)
      # }
    }
    names(mods.thrshld.spi) <- scn.nms

    mods.thrshld[[i]] <- append(mods.thrshld[[i]], mods.thrshld.spi)
  }
  return(mods.thrshld)
}




