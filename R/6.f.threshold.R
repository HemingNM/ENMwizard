#### 4.3.3 aplicar threshold
# TODO documentation
# TODO shorten name of arg "mxnt.mdls.preds.sp[...]"
#' Function Title (short description)
#'
#' General function description. A short paragraph (or more) describing what the function does.
#' @param mmp.spi Stack or brick of predictions to apply the threshold
#' @param pred.nm name of prediction to be appended to the final name. Usually "pres", "past" or "fut".
#' @param thrshld.i List of threshold criteria to be applied
#' @return Stack or brick of thresholded predictions
#' @examples
#' plot(mxnt.mdls.preds.lst[[1]][[4]]) # MaxEnt predictions, based on the model selection criteria
#' @export
f.thr <- function(mmp.spi, pred.nm = "", thrshld.i = 4:6) {
  if(is.null(path.mdls)){
    path.mdls <- "4_ENMeval.results"
    if(dir.exists(path.mdls)==F) dir.create(path.mdls)
  }
  # args <- 1:length(mmp.spi[["mxnt.mdls"]])
  mxnt.mdls <- mmp.spi[["mxnt.mdls"]]

  # pred.nm <- ifelse(pred.nm == "", pred.nm, paste0(".", pred.nm))
  # m.pred.n <- ifelse(pred.nm == "", "mxnt.pred", paste0("mxnt.pred", pred.nm))
  m.pred.n <- ifelse(pred.nm == "", "mxnt.pred", paste0(pred.nm))
  print(paste("pred.nm is ", pred.nm))
  print(paste("m.pred.n is ", m.pred.n))

  ##  wrong here
  pred.r <- mmp.spi[[match(m.pred.n, names(mmp.spi))]] # , fixed=T # [pred.i]
  pred.args <- mmp.spi$pred.args
  mod.pred.nms <- names(pred.r)

  if(sum(grepl("AvgAICc", mod.pred.nms))>0) {
    # args.aicc <- 1:(length(args)-4)
    wv <- mmp.spi[[1]][order(mmp.spi[[1]]$delta.AICc[grep("Mod.AICc", mmp.spi[[1]]$sel.cri)]),"w.AIC"]
  } # else {wv <- rep(0, max(grep("Mod.AICc", mmp.spi[[1]]$sel.cri))) }


  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

  thrshld.path <- paste(path.mdls, outpt, "Mdls.thrshld", sep='/')

  if(dir.exists(thrshld.path)==F) dir.create(thrshld.path)

  thrshld.nms <- c("fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd")[thrshld.i]
  thrshld.crit <- rownames(mxnt.mdls[[1]]@results)[grepl("Cloglog", rownames(mxnt.mdls[[1]]@results)) &
                                                     grepl("threshold", rownames(mxnt.mdls[[1]]@results))][thrshld.i]

  ## extract values of threshold from each model and criteria # sapply(mxnt.mdls[1:length(args)]
  thrshld.crit.v <- as.data.frame(matrix(data=sapply(thrshld.crit, function(y) sapply(mxnt.mdls[1:length(mxnt.mdls)], function(x) x@results[rownames(x@results) == y]) ),
                                         ncol = length(thrshld.i)))

  # thrshld.mod.crt <- thrshld.crit.v[1]
  # grep("LowAICc", names(pred.r))

  ## TO DO - change order of all stacks to c("Mod.AvgAICc", "Mod.LowAICc", "Mod.Mean.ORmin", "Mod.Mean.OR10", "Mod.Mean.AUCmin", "Mod.Mean.AUC10")
  thrshld.mod.crt <- data.frame(rbind(
    if(sum(grepl("AvgAICc", names(pred.r)))>0){ # # 1:length(args.aicc)
      apply(thrshld.crit.v[grep("Mod.AICc", mmp.spi[[1]]$sel.cri),], 2, function(x) stats::weighted.mean(x, wv)) ### check if is raster
    } else {NA}, # compute avg.thrshld from each criteria weighted by model importance (AICc W)
    if(sum(grepl("LowAICc", names(pred.r)))>0){
      thrshld.crit.v[grep("Mod.AICc_1$", mmp.spi[[1]]$sel.cri),]
    } else {NA},
    if(sum(grepl("Mean.ORmin", names(pred.r)))>0){
      thrshld.crit.v[grep("Mean.ORmin", mmp.spi[[1]]$sel.cri),]
    }else {NA},
    if(sum(grepl("Mean.OR10", names(pred.r)))>0){
      thrshld.crit.v[grep("Mean.OR10", mmp.spi[[1]]$sel.cri),]
    } else {NA},
    if(sum(grepl("Mean.AUCmin", names(pred.r)))>0){
      thrshld.crit.v[grep("Mean.AUCmin", mmp.spi[[1]]$sel.cri),]
    } else {NA},
    if(sum(grepl("Mean.AUC10", names(pred.r)))>0){
      thrshld.crit.v[grep("Mean.AUC10", mmp.spi[[1]]$sel.cri),]
    } else {NA} ) )

  # thrshld.mod.crt <- cbind(thrshld.mod.crt, thrshld.mod.crt)
  thrshld.mod.crt <- as.data.frame(thrshld.mod.crt[!is.na(thrshld.mod.crt[,1]),])
  row.names(thrshld.mod.crt) <- mod.pred.nms
  # col.names(thrshld.mod.crt) <- thrshld.nms

  brick.nms.t <- paste0("mxnt.pred", pred.nm, ".", thrshld.nms)
  brick.nms.t.b <- paste0("mxnt.pred", pred.nm, ".b.", thrshld.nms)

  for(t in 1:length(thrshld.crit)){
    mod.pred.nms.t <- paste(paste0(mod.pred.nms, ".", pred.nm), thrshld.nms[t], sep=".")
    mod.pred.nms.b <- paste(paste0(mod.pred.nms, ".", pred.nm, ".b"), thrshld.nms[t], sep=".")

    pred.t <- pred.r

    for(m in 1:length(mod.pred.nms)){
      pred.t[[m]][pred.t[[m]] < thrshld.mod.crt[m,t]] <- 0
    }
    names(pred.t) <- mod.pred.nms.t
    assign(brick.nms.t[t],
           raster::writeRaster(x = pred.t,
                       filename = paste(thrshld.path, paste0("mxnt.pred", gsub(".mxnt.pred","", paste0(".",pred.nm)), ".", thrshld.nms[t], ".grd"), sep='/'),
                       format = "raster", overwrite = T))
    # create presence only raster
    for(m in 1:length(mod.pred.nms)){
      pred.t[[m]][pred.t[[m]] >= thrshld.mod.crt[m,t]] <- 1
    }
    assign(brick.nms.t.b[t],
           raster::writeRaster(x = pred.t,
                       filename = paste(thrshld.path, paste0("mxnt.pred", gsub(".mxnt.pred","", paste0(".",pred.nm)), ".b.", thrshld.nms[t], ".grd"), sep='/'),
                       format = "raster", overwrite = T))
  }
  mods.t <- lapply(brick.nms.t, function(x) get(x))
  mods.t.b <- lapply(brick.nms.t.b, function(x) get(x))
  names(mods.t) <- thrshld.nms
  names(mods.t.b) <- thrshld.nms
  return(list(continuous=mods.t, binary=mods.t.b))
}

# mods.thrshld <- f.thr(mmp.spi, thrshld.i = 4:6, pred.args, path.mdls)
# plot(mods.thrshld[[1]][[2]]) # continuous
# plot(mods.thrshld[[2]][[2]]) # binary

#' Function Title (short description)
#'
#' General function description. A short paragraph (or more) describing what the function does.
#' @param mmp.spl List of stack or brick of predictions to apply the threshold
#' @inheritParams f.thr
#' @return List of stack or brick of thresholded predictions
#' @examples
#' plot(mxnt.mdls.preds.lst[[1]][[4]]) # MaxEnt predictions, based on the model selection criteria
#' @export
f.thr.batch <- function(mmp.spl, pred.nm="", thrshld.i = 4:6){
  path.res <- "4_ENMeval.results"
  if(dir.exists(path.res)==F) dir.create(path.res)
  path.sp.m <- paste0("Mdls.", names(mmp.spl))
  path.mdls <- paste(path.res, path.sp.m, sep="/")

  # thrshld for each species
  mods.thrshld <- vector("list", length(mmp.spl))
  names(mods.thrshld) <- names(mmp.spl)
  pred.nm <- ifelse(pred.nm == "", "mxnt.pred", pred.nm)

  for(i in 1:length(mmp.spl)){ # species i
    if(dir.exists(path.mdls[i])==F) dir.create(path.mdls[i])
    #scn.ind <- grep(n.pred.nm, names(mmp.spl[[i]]))
    scn.ind <- grep(pred.nm, names(mmp.spl[[i]]))
    scn.nms <- names(mmp.spl[[i]])[scn.ind]
    mods.thrshld.spi <- stats::setNames(vector("list", length(scn.nms)), scn.nms)
    for(j in seq_along(scn.ind)){ # climatic scenario
      mods.thrshld.spi[[j]] <- f.thr(mmp.spl[[i]], pred.nm = scn.nms[j], thrshld.i) # , path.mdls[i]
    }
    mods.thrshld[[i]] <- append(mods.thrshld[[i]], mods.thrshld.spi)
  }
  return(mods.thrshld)
}


#### 4.8.5 threshold for past and future pred
# mods.thrshld.lst <- f.thr.batch(mxnt.mdls.preds.pf, pred.nm="", thrshld.i = 4:6)


