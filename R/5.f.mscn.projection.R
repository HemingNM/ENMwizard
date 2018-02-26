### 4.8 predictions for future and past
### functions to predict areas based on fitted models
## TODO remove arg "pred.args", replace by mcm$pred.args
## TODO remove args "wAICsum", "randomseed", "responsecurves", "arg1", and "arg2", by storing it on
# replace by mcm$pred.args

#' Projecting Calibrated MaxEnt Models
#'
#' This function will read an object returned by "mxnt.cp", read the calibrated models and project into
#' new areas/climatic scenarios. These new projections will be returned together with (appended to)
#' the original object
#' @param mcm Objects returned by "mxnt.cp", containing calibrated models.
#' @param pred.nm Character. Prefix to add to projection name (e.g. "fut" or "past")
#' @param a.proj A Raster* object or a data.frame where models will be projected. Argument 'x' of dismo::predict
# #' @param numCores Number of cores to use for parallelization. If set to 1, no paralellization is performed
#' @inheritParams mxnt.cp
#' @return A list containing all the items returned from function "mxnt.cp", plus the projection specified in a.proj.
#' Each projection is a raster stack containing model projections ('mxnt.preds'), where each layer is a projection based on
#' a specific model selection criteria (i.e. AvgAICc, LowAICc, Mean.ORmin, Mean.OR10, Mean.AUCmin, Mean.AUC10)
# #' @examples
#' @export
mxnt.p <- function(mcm, sp.nm, pred.nm="fut", a.proj, formt = "raster",numCores=1,parallelTunning=TRUE){ # , #, ENMeval.occ.results, occ.b.env, occ.locs,
                        # pred.args = c("outputformat=cloglog", "doclamp=true", "pictures=true"),
                        # wAICsum=0.99, randomseed=F, responsecurves=T, arg1='noaddsamplestobackground', arg2='noautofeature'){ # wAICsum=0.99,

  path.res <- "4_ENMeval.results"
  if(dir.exists(path.res)==FALSE) dir.create(path.res)
  path.mdls <- paste(path.res, paste0("Mdls.", sp.nm), sep="/")
  if(dir.exists(path.mdls)==FALSE) dir.create(path.mdls)
  pred.args <- mcm$pred.args

  xsel.mdls <- mcm$selected.mdls # mdl.arg[[2]]
  f <- factor(xsel.mdls$features)
  beta <- xsel.mdls$rm
  args.all <- mcm$mxnt.args
  args.aicc <- args.all[grep("Mod.AIC", xsel.mdls$sel.cri)] #[1:2]
  print(data.frame(features=f, beta, row.names = xsel.mdls$sel.cri))

  mod.nms <- paste(xsel.mdls[,"sel.cri"]) # paste0("Mod.", c(1:length(args.aicc), "Mean.ORmin", "Mean.OR10", "Mean.AUCmin", "Mean.AUC10"))

  ## TO DO - change order of all stacks to c("Mod.AvgAICc", "Mod.LowAICc", "Mod.Mean.ORmin", "Mod.Mean.OR10", "Mod.Mean.AUCmin", "Mod.Mean.AUC10")
  mod.pred.nms <- c("Mod.AvgAICc", "Mod.LowAICc", mod.nms[(length(args.aicc)+1):length(args.all)])
  mod.preds <- raster::stack() #vector("list", length(mod.pred.nms))

  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

  if(dir.exists(paste(path.mdls, outpt, sep='/'))==FALSE) dir.create(paste(path.mdls, outpt, sep='/'))

  mxnt.mdls <- mcm$mxnt.mdls

  # if(pred.nm != "") {pred.nm <- paste0(".", pred.nm)}

  #### 4.3.2 predictions

  # AIC AVG model
  avg.m.path <- paste(path.mdls, outpt, mod.pred.nms[1], sep='/') # paste0("4_ENMeval.results/selected.models/cloglog/", mod.pred.nms[2])
  if(dir.exists(avg.m.path)==FALSE) dir.create(avg.m.path)

  ##### list of models to PREDICT
  mod.all <- vector("list")
  filename.aicc <- paste(avg.m.path, mod.nms, paste0(mod.nms, ".", pred.nm,".grd"), sep='/')[1:length(args.aicc)]

  # path2file <-paste(path.mdls, outpt, mod.nms[(length(args.aicc)+1):length(args.all)], sep='/')
  filename.au.om <- paste(path.mdls, outpt, mod.nms[(length(args.aicc)+1):length(args.all)],
                        paste0(mod.nms[(length(args.aicc)+1):length(args.all)], ".", pred.nm,".grd"), sep='/')

  filename <- c(filename.aicc, filename.au.om)

  if(numCores>1&parallelTunning){

    cl<-parallel::makeCluster(numCores)

    mod.all <- parallel::clusterApply(cl, seq_along(args.all), function(i,mxnt.mdls,a.proj,pred.args,filename,formt) {
      resu <- dismo::predict(mxnt.mdls[[i]], a.proj, args = pred.args, progress = 'text',file = filename[i], format = formt, overwrite = TRUE)
      return(resu)}, mxnt.mdls, a.proj, pred.args, filename, formt) #) ## fecha for or lapply

    parallel::stopCluster(cl)

  }else{

    mod.all <- lapply(seq_along(args.all), function(i) {
      resu <- dismo::predict(mxnt.mdls[[i]], a.proj, args = pred.args, progress = 'text',file = filename[i], format = formt, overwrite = TRUE)
      return(resu)}) #) ## fecha for or lapply
  }

  mod.preds <- raster::stack() #vector("list", length(mod.pred.nms))
  #### AIC AVG model
  {
    #### 4.3.2.1.2 create model averaged prediction (models*weights, according to model selection)
    # create vector of model weights
    wv <- xsel.mdls[order(xsel.mdls$delta.AICc),"w.AIC"][seq_along(args.aicc)]

    ### stack prediction rasters (to create Average Model prediction)
    path2stk <- paste(avg.m.path, mod.nms[seq_along(args.aicc)], sep='/')
    filename <- paste(path2stk, paste0(mod.nms[seq_along(args.aicc)], ".", pred.nm,".grd"), sep='/')
    Mod.AICc.stack <- raster::stack(filename)

    # create averaged prediction map
    print(mod.pred.nms[1])
    mod.preds <- raster::addLayer(mod.preds, raster::writeRaster(raster::mask((sum(Mod.AICc.stack*wv, na.rm = T)/sum(wv)), a.proj[[1]]),
                                                 filename = paste(avg.m.path, paste0(mod.pred.nms[1], ".", pred.nm,".grd"), sep='/'),
                                                 format = formt, overwrite = T) )
    names(mod.preds)[raster::nlayers(mod.preds)] <- mod.pred.nms[1]
  }

  #### Low AIC
  # if(i == 1) # usar if(low = T) pra escolher o low aic ou if(grep("low", Mod.pred))
  {
    path2file <- paste(path.mdls, outpt, mod.pred.nms[2], sep='/')
    filename <- paste(path2file, paste0(mod.pred.nms[2], ".", pred.nm,".grd"), sep='/')
    if(dir.exists(path2file) == FALSE) dir.create(path2file)

    #### 4.3.2.1.1 create Low AIC model prediction on a specific path
    print(mod.pred.nms[2])
    mod.preds <- raster::addLayer(mod.preds, raster::writeRaster(mod.all[[1]],
                                                 filename = filename,
                                                 format = formt, overwrite = T) )
    names(mod.preds)[raster::nlayers(mod.preds)] <- mod.pred.nms[2]
  }


  #### AUC OmR models
  if(length(args.all) > length(args.aicc)){

    for(i in (length(args.aicc)+1):length(args.all)){
      # path2file <-paste(path.mdls, outpt, mod.nms[i], sep='/')
      # filename <- paste(path2file, paste0(mod.nms[i], ".", pred.nm,".grd"), sep='/')
      # if(dir.exists(path2file) == FALSE) dir.create(path2file)
      print(mod.nms[i])
      mod.preds <- raster::addLayer(mod.preds, mod.all[[i]] )
      names(mod.preds)[raster::nlayers(mod.preds)] <- mod.nms[i]
    }
  }

  # also changed line 51
  # mcm <- append(mcm, stats::setNames(list(mod.preds) , paste0("mxnt.preds", pred.nm)))

  # if(is.null(mcm$mxnt.preds)){
  #   mcm$mxnt.preds <- stats::setNames(list(mod.preds) , paste0(pred.nm))
  # } else {
    mcm$mxnt.preds <- append(mcm$mxnt.preds, stats::setNames(list(mod.preds) , paste0(pred.nm)))
  # }

  return(mcm)
}

# #### for several species
# #' Projecting calibrated MaxEnt models for several species
# #'
# #' This function will read an object returned by "mxnt.cp.batch", read the calibrated models and project into
# #' new areas/climatic scenarios. These new projections will be returned together with (appended to)
# #' each element (species) the original object
# #' @param mcm.l A list of objects returned by "mxnt.cp", containing calibrated models.
# #' @param a.proj.l A list of Raster* objects or data.frames where models will be projected. Argument 'x' of dismo::predict
# #' @inheritParams mxnt.p
# #' @inheritParams mxnt.cp.batch
# #' @return A list of objects returned from function "mxnt.p"
# #' @examples
# #' mxnt.mdls.preds <- mxnt.p.batch(mcm.l = mxnt.mdls.preds.lst[1],
# #                                     pred.nm ="fut", a.proj.l = areas.projection)
# #' @export
# mxnt.p.batch <- function(mcm.l, pred.nm="fut", a.proj.l, formt = "raster",numCores=1){ #, #, ENMeval.occ.results.lst, occ.b.env.lst, occ.locs.lst,
#                               # pred.args = c("outputformat=cloglog", "doclamp=true", "pictures=true"),
#                               # wAICsum=0.99, randomseed=F, responsecurves=T, arg1='noaddsamplestobackground', arg2='noautofeature'){ #wAICsum=0.99,
#'
#   # path.res <- "4_ENMeval.results"
#   # if(dir.exists(path.res)==F) dir.create(path.res)
#   # path.mdls <- paste(path.res, paste0("Mdls.", names(mcm.l)), sep="/")
#   mcmp.l <- vector("list", length(mcm.l))
#   names(mcmp.l) <- names(mcm.l)
#   for(i in 1:length(mcmp.l)){
#     # if(dir.exists(path.mdls[i])==F) dir.create(path.mdls[i])
#     cat(c(names(mcm.l)[i], "\n"))
#     # print(paste(names(mcm.l)[i]))
#     # compute final models and predictions
#     mcmp.l[[i]] <- mxnt.p(mcm = mcm.l[[i]], #ENMeval.occ.results = ENMeval.occ.results.lst[[i]],
#                                   sp.nm = names(mcm.l)[i], pred.nm = pred.nm, a.proj = a.proj.l[[i]],
#                                        # occ.b.env = occ.b.env.lst[[i]], occ.locs = occ.locs.lst[[i]],
#                                   formt = formt,
#                                   numCores=numCores) # , #pred.args = pred.args,
#                                        # wAICsum = wAICsum,
#                                        # randomseed = randomseed, responsecurves = responsecurves, arg1 = arg1, arg2 = arg2)
#     # mcm.l[[i]]$pred.args <- pred.args
#   }
#   return(mcmp.l)
#   }



#' Projecting calibrated MaxEnt models for several species onto multiple environmental scenarios
#'
#' This function will read an object returned by "mxnt.cp.batch", read the calibrated models and project into
#' several environmental (areas/climatic) scenarios (specified in a.proj.l). These new projections will be returned together with (appended to)
#' each element (species) the original object.
#' @param mcm.l A list of objects returned by "mxnt.cp", containing calibrated models.
#' @param a.proj.l A list of Raster* objects or data.frames where models will be projected. Argument 'x' of dismo::predict
#' @inheritParams mxnt.p
#' @inheritParams mxnt.cp.batch
#' @return A 'mcmp.spl' object. A list of objects returned from function "mxnt.p", containing the new (multiple) projections for each element (species) of the list
#' @examples
#' mxnt.mdls.preds.pf <- mxnt.p.batch.Mscn(mxnt.mdls.preds.lst, a.proj.l = area.projection.pf)
#' @export

mxnt.p.batch.mscn <- function(mcm.l, a.proj.l, formt = "raster", numCores=1, parallelTunning=TRUE){ #, # cores=2, #, pred.nm="fut", ENMeval.occ.results.lst, occ.b.env.lst, occ.locs.lst,
  mdl.names <- names(mcm.l)

  if(numCores>1 & parallelTunning==FALSE){

    cl<-parallel::makeCluster(numCores)
    parallel::clusterExport(cl,list("mxnt.p"))

    mcmp.l <- parallel::clusterApply(cl, seq_along(mcm.l),

                                    function(i, a.proj.l, mcm.l, formt){
                                      mxnt.preds.spi <- list()
                                      mcm <- mcm.l[[i]]
                                      sp.nm <- names(mcm.l)[i]
                                      pred.nm <- names(a.proj.l[[i]])
                                      a.proj = a.proj.l[[i]]

                                       for(j in 1:length(a.proj)){
                                         mxnt.preds.spi[j] <- mxnt.p(mcm = mcm,
                                                                     sp.nm = sp.nm, pred.nm = pred.nm[j],
                                                                     a.proj = a.proj[[j]],
                                                                     formt = formt,parallelTunning=parallelTunning,
                                                                     numCores = numCores)$mxnt.preds[length(mcm$mxnt.preds) + 1]
                                       }
                                      names(mxnt.preds.spi) <- paste0(names(a.proj))

                   # resu <- append(mcm.l[[i]], mxnt.preds.spi)
                   # mcm.l[[i]]$mxnt.preds <- mxnt.preds.spi
                   mcm.l[[i]]$mxnt.preds <- append(mcm.l[[i]]$mxnt.preds, mxnt.preds.spi)
                   # resu <- mcm.l[[i]]
                   return(mcm.l[[i]])}, a.proj.l, mcm.l, formt)


    parallel::stopCluster(cl)

  }else{

    mcmp.l <- lapply(seq_along(mcm.l),

                              function(i,a.proj.l,mcm.l,formt){
                                cat(c(names(mcm.l)[i], "\n"))
                                mxnt.preds.spi <- list()
                                mcm <- mcm.l[[i]]
                                sp.nm <- names(mcm.l)[i]
                                pred.nm <- names(a.proj.l[[i]])
                                a.proj <- a.proj.l[[i]]

                                for(j in 1:length(a.proj)){
                                  cat(c("\n", paste0("mxnt.pred.", names(a.proj)[j]), "\n",
                                        "projection ", j, " of ", length(a.proj), "\n"))

                                  mxnt.preds.spi[j] <- mxnt.p(mcm = mcm,
                                                              sp.nm = sp.nm, pred.nm = pred.nm[j],
                                                              a.proj = a.proj[[j]],
                                                              formt = formt,parallelTunning=parallelTunning,
                                                              numCores = numCores)$mxnt.preds[length(mcm$mxnt.preds) + 1]
                                }
                                names(mxnt.preds.spi) <- paste0(names(a.proj))

                                # resu <- append(mcm.l[[i]], mxnt.preds.spi)
                                # mcm.l[[i]]$mxnt.preds <- mxnt.preds.spi
                                mcm.l[[i]]$mxnt.preds <- append(mcm.l[[i]]$mxnt.preds, mxnt.preds.spi)
                                # resu <- mcm.l[[i]]
                                return(mcm.l[[i]])}, a.proj.l, mcm.l, formt)

  }

  # names(mcm.l) <- mdl.names
  names(mcmp.l) <- mdl.names

  # return(mcm.l)
  return(mcmp.l)
}
