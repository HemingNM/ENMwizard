### 4.8 predictions for future and past
### functions to predict areas based on fitted models

#' Project calibrated MaxEnt models
#'
#' This function will read an object returned by "calib_mdl", read the calibrated models and project into
#' new areas/climatic scenarios. These new projections will be returned together with (appended to)
#' the original object.
#'
#' @param mcm Objects returned by "calib_mdl", containing calibrated models.
#' @param pred.nm Character. Prefix to add to projection name (e.g. "fut" or "past")
#' @param a.proj A Raster* object or a data.frame where models will be projected. Argument 'x' of dismo::predict
# #' @param numCores Number of cores to use for parallelization. If set to 1, no paralellization is performed
#' @inheritParams calib_mdl
#' @inheritParams raster::writeRaster
#' @seealso \code{\link{mod_sel}}, \code{\link{calib_mdl}}, \code{\link{calib_mdl_b}}, \code{\link[dismo]{maxent}},
#' \code{\link[ENMeval]{ENMevaluate}}, \code{\link{ENMevaluate_b}}, \code{\link{proj_mdl_b}}
#' @return A list containing all the items returned from function "calib_mdl", plus the projection specified in a.proj.
#' Each projection is a raster stack containing model projections ('mxnt.preds'), where each layer is a projection based on
#' a specific model selection criteria (i.e. AvgAIC, LowAIC, avg.test.orMTP, avg.test.or10pct, avg.test.AUC.MTP, avg.test.AUC10pct)
# #' @examples
#' @export
proj_mdl <- function(mcm, sp.nm="species", pred.nm="fut", a.proj, format = "raster",
                   numCores = 1, parallelTunning = TRUE){ # , #, ENMeval.occ.results, occ.b.env, occ.locs,
  path.res <- "3_out.MaxEnt"
  if(dir.exists(path.res)==FALSE) dir.create(path.res)
  path.mdls <- paste(path.res, paste0("Mdls.", sp.nm), sep="/")
  if(dir.exists(path.mdls)==FALSE) dir.create(path.mdls)
  pred.args <- mcm$pred.args

  xsel.mdls <- mcm$selected.mdls # [order(mcm$selected.mdls$delta.AICc),] # mdl.arg[[2]]
  args.all <- mcm$mxnt.args

  mod.nms <- paste0("Mod_", format(xsel.mdls[, "rm"], digits=2), "_", xsel.mdls[, "features"]) #
  # mod.nms <- paste0("Mod.", xsel.mdls$sel.cri)
  # mod.nms <- paste0("Mod.", xsel.mdls[, "settings"]) #
  mod.nms2 <- paste0("Mod.", xsel.mdls$sel.cri)
  # # mod.nms2 <- gsub(paste0("AIC_", 1:length(xsel.mdls$sel.cri), "." , collapse = "|"), "", mod.nms)
  ens.meth <- c("AIC_", "WAAUC_", "EBPM_", "ESOR_")
  ens2gsub <- paste0(ens.meth, rep(1:length(xsel.mdls$sel.cri), each= length(ens.meth)), collapse = "|")
  mod.nms2 <- gsub(ens2gsub, "", mod.nms2)
  mod.nms2 <- gsub(paste0("\\.{", 1:length(xsel.mdls$sel.cri), "}" , collapse = "|"), ".", mod.nms2)


  mod.preds <- raster::stack() #vector("list", length(mod.pred.nms))

  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

  if(dir.exists(paste(path.mdls, outpt, sep='/'))==FALSE) dir.create(paste(path.mdls, outpt, sep='/'))

  mxnt.mdls <- mcm$mxnt.mdls

  #### 4.3.2 predictions

  # # AIC AVG model
  # if(length(grep("AvgAIC", mcm$mSel))>0) {
  #   avg.m.path <- paste(path.mdls, outpt, "Mod.AvgAIC", sep='/') # paste0("3_out.MaxEnt/selected.models/cloglog/", mod.pred.nms[2])
  #   if(dir.exists(avg.m.path)==FALSE) dir.create(avg.m.path)
  #   filename.aicc <- paste(avg.m.path, paste0("Mod.AvgAIC", ".", pred.nm,".grd"), sep='/')#[1:length(args.aicc)]
  # }
  filename.or.auc.laic <- paste(path.mdls, outpt, mod.nms,
                    paste0(mod.nms, ".", pred.nm,".grd"), sep='/')

  ##### list of models to PREDICT
  # mod.all <- vector("list")

  if(numCores>1 & parallelTunning){

    cl <- parallel::makeCluster(numCores)

    mod.all <- parallel::clusterApply(cl, seq_along(args.all), function(i, mxnt.mdls, a.proj, pred.args, filename.or.auc.laic, format) {
      resu <- dismo::predict(mxnt.mdls[[i]], a.proj, args = pred.args, progress = 'text', file = filename.or.auc.laic[i], format = format, overwrite = TRUE)
      return(resu)}, mxnt.mdls, a.proj, pred.args, filename.or.auc.laic, format) #) ## fecha for or lapply

    parallel::stopCluster(cl)

  } else {

    mod.all <- lapply(seq_along(args.all), function(i) {
      resu <- dismo::predict(mxnt.mdls[[i]], a.proj, args = pred.args, progress = 'text', file = filename.or.auc.laic[i], format = format, overwrite = TRUE)
      return(resu)}) #) ## fecha for or lapply
  }
  mod.all <- stats::setNames(mod.all, names(mxnt.mdls))

  #### Ensemble models (AvgAIC, WAAUC, EBPM)
  mod.preds <- ensemble_projs(mcm, filename.or.auc.laic, a.proj, path.mdls, outpt, pred.nm, format)

  #### AUC OR & LowAIC models
  args.or.auc.laic <- grep("All|LowAIC|AUC10|AUCmtp|OR10|ORmtp", mod.nms2)
  # args.or.auc.laic <- grep("All|LowAIC|AUC10|AUCmtp|OR10|ORmtp", xsel.mdls$sel.cri)
  for(i in args.or.auc.laic){
    mod.preds <- raster::addLayer(mod.preds, mod.all[[i]] )
    names(mod.preds)[raster::nlayers(mod.preds)] <- mod.nms2[i]# gsub(paste0("AIC_", 1:length(xsel.mdls$sel.cri), "." , collapse = "|"), "", mod.nms[i])
    }

  mcm$mxnt.preds <- append(mcm$mxnt.preds, stats::setNames(list(mod.preds) , paste0(pred.nm)))
  return(mcm)
}


#' Project calibrated MaxEnt models for several species onto multiple environmental scenarios
#'
#' This function will read an object returned by "calib_mdl_b", read the calibrated models and project into
#' several environmental (areas/climatic) scenarios (specified in a.proj.l). These new projections will be returned together with (appended to)
#' each element (species) the original object.
#'
#' @param mcm.l A list of objects returned by "calib_mdl", containing calibrated models.
#' @param a.proj.l A list of Raster* objects or data.frames where models will be projected. Argument 'x' of dismo::predict
#' @inheritParams proj_mdl
#' @inheritParams calib_mdl_b
#' @seealso \code{\link{mod_sel}}, \code{\link{calib_mdl}}, \code{\link{calib_mdl_b}}, \code{\link[dismo]{maxent}},
#' \code{\link[ENMeval]{ENMevaluate}}, \code{\link{ENMevaluate_b}}, \code{\link{proj_mdl}}
#' @return A 'mcmp.spl' object. A list of objects returned from function "proj_mdl", containing the new (multiple) projections for each element (species) of the list
#' @examples
#' \dontrun{
#' mxnt.mdls.preds.pf <- proj_mdl_b(mcm.l = mxnt.mdls.preds.lst, a.proj.l = area.projection.pf)
#' }
#' @export
proj_mdl_b <- function(mcm.l, a.proj.l, format = "raster", numCores = 1, parallelTunning = TRUE){ #, # cores=2, #, pred.nm="fut", ENMeval.occ.results.lst, occ.b.env.lst, occ.locs.lst,
  mdl.names <- names(mcm.l)

  if(numCores>1 & parallelTunning==FALSE){

    cl<-parallel::makeCluster(numCores)
    parallel::clusterExport(cl,list("proj_mdl"))

    mcmp.l <- parallel::clusterApply(cl, seq_along(mcm.l),

                                     function(i, a.proj.l, mcm.l, format){
                                       mxnt.preds.spi <- list()
                                       mcm <- mcm.l[[i]]
                                       sp.nm <- names(mcm.l)[i]
                                       pred.nm <- names(a.proj.l[[i]])
                                       a.proj = a.proj.l[[i]]

                                       f <- factor(mcm$selected.mdls$features)
                                       beta <- mcm$selected.mdls$rm
                                       print(data.frame(features=f, beta, row.names = mcm$selected.mdls$sel.cri))

                                       for(j in 1:length(a.proj)){
                                         mxnt.preds.spi[j] <- proj_mdl(mcm = mcm,
                                                                     sp.nm = sp.nm, pred.nm = pred.nm[j],
                                                                     a.proj = a.proj[[j]],
                                                                     format = format, parallelTunning = parallelTunning,
                                                                     numCores = numCores)$mxnt.preds[length(mcm$mxnt.preds) + 1]
                                       }
                                       names(mxnt.preds.spi) <- paste0(names(a.proj))

                                       mcm.l[[i]]$mxnt.preds <- append(mcm.l[[i]]$mxnt.preds, mxnt.preds.spi)
                                       return(mcm.l[[i]])}, a.proj.l, mcm.l, format)


    parallel::stopCluster(cl)

  }else{

    mcmp.l <- lapply(seq_along(mcm.l),

                     function(i,a.proj.l,mcm.l,format){
                       cat(c("\n", names(mcm.l)[i], "\n"))
                       mxnt.preds.spi <- list()
                       mcm <- mcm.l[[i]]
                       sp.nm <- names(mcm.l)[i]
                       pred.nm <- names(a.proj.l[[i]])
                       a.proj <- a.proj.l[[i]]

                       f <- factor(mcm$selected.mdls$features)
                       beta <- mcm$selected.mdls$rm
                       print(data.frame(features=f, beta, row.names = mcm$selected.mdls$sel.cri))

                       for(j in 1:length(a.proj)){
                         cat(c("\n", "projection", j, "of", length(a.proj), "-",
                               names(a.proj)[j], "\n"))

                         mxnt.preds.spi[j] <- proj_mdl(mcm = mcm,
                                                     sp.nm = sp.nm, pred.nm = pred.nm[j],
                                                     a.proj = a.proj[[j]],
                                                     format = format, parallelTunning = parallelTunning,
                                                     numCores = numCores)$mxnt.preds[length(mcm$mxnt.preds) + 1]
                       }
                       names(mxnt.preds.spi) <- paste0(names(a.proj))

                       mcm.l[[i]]$mxnt.preds <- append(mcm.l[[i]]$mxnt.preds, mxnt.preds.spi)
                       return(mcm.l[[i]])}, a.proj.l, mcm.l, format)

  }

  names(mcmp.l) <- mdl.names
  return(mcmp.l)
}




#' Ensemble models (AvgAIC, WAAUC, EBPM, ESOR)
#'
#' This function will read an object returned by "calib_mdl", read the calibrated models and project into
#' new areas/climatic scenarios. These new projections will be returned together with (appended to)
#' the original object.
#'
#' @param mcm Objects returned by "calib_mdl", containing calibrated models.
#' @param pred.nm Character. Prefix to add to projection name (e.g. "current", "fut", or "past")
#' @param a.proj A Raster* object or a data.frame where models will be projected. Argument 'x' of dismo::predict
#' @seealso \code{\link{mod_sel}}, \code{\link{calib_mdl}}, \code{\link{calib_mdl_b}}, \code{\link[dismo]{maxent}},
#' \code{\link[ENMeval]{ENMevaluate}}, \code{\link{proj_mdl_b}}
#' @return A list containing all the items returned from function "calib_mdl", plus the projection specified in a.proj.
#' Each projection is a raster stack containing model projections ('mxnt.preds'), where each layer is a projection based on
#' a specific model selection criteria (i.e. AvgAIC, LowAIC, avg.test.orMTP, avg.test.or10pct, avg.test.AUC.MTP, avg.test.AUC10pct)
# #' @examples
#' @keywords internal
# #' @export
ensemble_projs <- function(mcm, filename.or.auc.laic, a.proj, path.mdls, outpt, pred.nm, format){
  mod.preds <- raster::stack()
  xsel.mdls <- mcm$selected.mdls

  ens <- c("AvgAIC", "WAAUC", "EBPM", "ESOR")
  ens.i <- grepl(paste0("^", mcm$mSel, collapse = "|^"), ens)
  if(sum(ens.i)>0){
    ens <- ens[ens.i]
    for(EM in ens){
      # print(EM)
      # filenames
      ens.m.path <- paste(path.mdls, outpt, paste0("Mod.", EM), sep='/') # paste0("3_out.MaxEnt/selected.models/cloglog/", mod.pred.nms[2])
      if(dir.exists(ens.m.path)==FALSE) dir.create(ens.m.path)
      filename.ens <- paste(ens.m.path, paste0("Mod.", EM, ".", pred.nm,".grd"), sep='/')#[1:length(args.aicc)]

      #### 4.3.2.1.2 create model averaged prediction (models*weights, according to model selection)
      argsEns <- grep(EM, xsel.mdls$sel.cri)
      # argsEns <- grep(EM, mcm$mSel)

      # create vector of model weights
      if(EM == "AvgAIC"){ #
        argsEns <- grep("AIC", xsel.mdls$sel.cri)
        wv <- xsel.mdls[argsEns,"w.AIC"]
      } else if(EM == "WAAUC"){ # WAAUC
        wv <- xsel.mdls[argsEns,"avg.test.AUC"]
      } else if(EM == "EBPM"){ # EBPM
        wv <- rep(1, length(argsEns))
      } else if(EM == "ESOR"){ # Cobos et al 2019
        wv <- rep(1, length(argsEns))
      }

      ### stack prediction rasters (to create Average Model prediction)
      filename.avg.stk <- filename.or.auc.laic[argsEns]
      Mod.ens.stack <- raster::stack(filename.avg.stk)
      # plot(Mod.ens.stack[[1]])
      # print(names(Mod.ens.stack))
      # plot(a.proj[[1]])
      # create averaged prediction map
      if(length(argsEns)>1){
        mod.preds <- raster::addLayer(mod.preds, raster::writeRaster(raster::mask((sum(Mod.ens.stack*wv, na.rm = T)/sum(wv)), a.proj[[1]]),
                                                                     filename = filename.ens,
                                                                     format = format, overwrite = T) )
      } else {
        mod.preds <- raster::addLayer(mod.preds, raster::writeRaster(raster::mask(Mod.ens.stack, a.proj[[1]]),
                                                                     filename = filename.ens,
                                                                     format = format, overwrite = T) )
      }
      names(mod.preds)[raster::nlayers(mod.preds)] <- paste0("Mod.", EM)
    }
    return(mod.preds)
  } else {
    # print("No models selected for ensembling. Returning empty stack.")
    return(mod.preds)
  }
}

# mod.preds <- ensemble_projs(mod.preds, mcm, filename.or.auc.laic, pred.nm)

