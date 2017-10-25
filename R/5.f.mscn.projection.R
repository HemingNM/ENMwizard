### 4.8 predictions for future and past
### functions to predict areas based on fitted models
## TODO remove arg "pred.args", replace by mxnt.c.mdls$pred.args
## TODO remove args "wAICsum", "randomseed", "responsecurves", "arg1", and "arg2", by storing it on
# replace by mxnt.c.mdls$pred.args

#' Projecting calibrated MaxEnt models
#'
#' This function will read an object returned by "mxnt.cp", read the calibrated models and project into
#' new areas/climatic scenarios. These new projections will be returned together with (appended to)
#' the original object
#' @param mxnt.c.mdls Objects returned by "mxnt.cp", containing calibrated models.
#' @param pred.nm Character. Prefix to add to projection name (e.g. "fut" or "past")
#' @inheritParams mxnt.cp
#' @return A list containing all the items returned from function "mxnt.cp", plus the projection specified in a_proj.
# #' @examples
#' @export
mxnt.p <- function(mxnt.c.mdls, sp.nm, pred.nm="fut", a_proj, formt = "raster"){ # , #, ENMeval_occ_results, occ_b_env, occ_locs,
                        # pred.args = c("outputformat=cloglog", "doclamp=true", "pictures=true"),
                        # wAICsum=0.99, randomseed=F, responsecurves=T, arg1='noaddsamplestobackground', arg2='noautofeature'){ # wAICsum=0.99,

  # { #library(xlsx)
  #   library(dismo)
  #   library(rJava)
  #   library(raster)}

  path.res <- "4_ENMeval_results"
  if(dir.exists(path.res)==F) dir.create(path.res)
  path.mdls <- paste(path.res, paste0("Mdls_", sp.nm), sep="/")
  if(dir.exists(path.mdls)==F) dir.create(path.mdls)
  pred.args <- mxnt.c.mdls$pred.args

  # x <- ENMeval_occ_results
  # x <- mxnt.c.mdls$selected.mdls#[grep("AICc", mxnt.c.mdls$selected.mdls$sel.cri),]
  xsel.mdls <- mxnt.c.mdls$selected.mdls # mdl.arg[[2]]
  f <- factor(xsel.mdls$features)
  beta <- xsel.mdls$rm
  # mdl.arg <- f.args(x=x, wAICsum=wAICsum, randomseed=randomseed, responsecurves=responsecurves, arg1=arg1, arg2=arg2) # sum(x$w.AIC)
  # args.all <- mdl.arg[[1]]
  args.all <- mxnt.c.mdls$mxnt.args
  # args.all <- paste(paste0(arg1), paste0(arg2),
  #                   ifelse(grepl("H", f), paste("hinge"), paste("nohinge")),
  #                   ifelse(grepl("L", f), paste("linear"), paste("nolinear")),
  #                   ifelse(grepl("Q", f), paste("quadratic"), paste("noquadratic")),
  #                   ifelse(grepl("P", f), paste("product"), paste("noproduct")),
  #                   ifelse(grepl("T", f), paste("threshold"), paste("nothreshold")),
  #                   paste0("betamultiplier=", beta),
  #                   paste0("responsecurves=", ifelse(responsecurves==T, "true", "false")),
  #                   paste0("randomseed=", ifelse(randomseed==T, "true", "false")),
  #                   sep = ",")
  # # args.aicc <- args.all[1:(length(args.all)-4)]
  # # print(data.frame(features=f, beta, row.names = c(paste("AICc", 1:length(args.aicc)),
  # #                                                 "Mean.ORmin", "Mean.OR10", "Mean.AUC_min", "Mean.AUC_10")))
  args.aicc <- args.all[grep("Mod_AIC", xsel.mdls$sel.cri)] #[1:2]
  print(data.frame(features=f, beta, row.names = xsel.mdls$sel.cri))


  # # exportar planilha de resultados
  # write.xlsx(xsel.mdls, paste0(path.mdls,"/sel.mdls.xlsx"))
  # res.tbl <- xsel.mdls[,c("sel.cri", "features","rm","AICc", "w.AIC", "nparam", "rankAICc", "Mean.OR10", "Mean.ORmin", "Mean.AUC")]
  # colnames(res.tbl) <- c("Optimality criteria", "FC", "RM", "AICc", "wAICc", "NP", "Rank", "OR10", "ORLPT", "AUC")
  # write.xlsx(res.tbl, paste0(path.mdls,"/sel.mdls.smmr.xlsx"))

  mod.nms <- paste(xsel.mdls[,"sel.cri"]) # paste0("Mod_", c(1:length(args.aicc), "Mean.ORmin", "Mean.OR10", "Mean.AUC_min", "Mean.AUC_10"))

  ## TO DO - change order of all stacks to c("Mod_AvgAICc", "Mod_LowAICc", "Mod_Mean.ORmin", "Mod_Mean.OR10", "Mod_Mean.AUC_min", "Mod_Mean.AUC_10")
  mod.pred.nms <- c("Mod_AvgAICc", "Mod_LowAICc", mod.nms[(length(args.aicc)+1):length(args.all)])
  mod.preds <- raster::stack() #vector("list", length(mod.pred.nms))

  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

  if(dir.exists(paste(path.mdls, outpt, sep='/'))==F) dir.create(paste(path.mdls, outpt, sep='/'))

  # mxnt.mdls <- vector("list", length(args.all))
  mxnt.mdls <- mxnt.c.mdls$mxnt.mdls

  if(pred.nm != "") {pred.nm <- paste0(".", pred.nm)}

  #### 4.3.2 predictions

  # AIC AVG model
  avg.m.path <- paste(path.mdls, outpt, mod.pred.nms[1], sep='/') # paste0("4_ENMeval_results/selected_models/cloglog/", mod.pred.nms[2])
  if(dir.exists(avg.m.path)==F) dir.create(avg.m.path)

  ##### list of models to average
  mod.avg.i <- vector("list")
  filename <- paste(avg.m.path, mod.nms, paste0(mod.nms, pred.nm, ".grd"), sep='/')
  lapply(seq_along(args.aicc), function(i) {
    mod.avg.i[[i]] <<- dismo::predict(mxnt.mdls[[i]], a_proj, args = pred.args, progress = 'text',
                               file = filename[i], format = formt, overwrite = T)
  }) #) ## fecha for or lapply

  mod.preds <- raster::stack() #vector("list", length(mod.pred.nms))
  {
    #### 4.3.2.1.2 create model averaged prediction (models*weights, according to model selection)
    # create vector of model weights
    wv <- xsel.mdls[order(xsel.mdls$delta.AICc),"w.AIC"][seq_along(args.aicc)]

    ### stack prediction rasters (to create Average Model prediction)
    path2stk <- paste(avg.m.path, mod.nms[seq_along(args.aicc)], sep='/')
    filename <- paste(path2stk, paste0(mod.nms[seq_along(args.aicc)], pred.nm, ".grd"), sep='/')
    Mod_AICc.stack <- raster::stack(filename)

    # create averaged prediction map
    print(mod.pred.nms[1])
    mod.preds <- raster::addLayer(mod.preds, raster::writeRaster(raster::mask((sum(Mod_AICc.stack*wv, na.rm = T)/sum(wv)), a_proj[[1]]),
                                                 filename = paste(avg.m.path, paste0(mod.pred.nms[1], pred.nm, ".grd"), sep='/'),
                                                 format = formt, overwrite = T) )
    names(mod.preds)[raster::nlayers(mod.preds)] <- mod.pred.nms[1]
  }

  #### Low AIC
  # if(i == 1) # usar if(low = T) pra escolher o low aic ou if(grep("low", Mod_pred))
  {
    path2file <- paste(path.mdls, outpt, mod.pred.nms[2], sep='/')
    filename <- paste(path2file, paste0(mod.pred.nms[2], pred.nm, ".grd"), sep='/')
    if(dir.exists(path2file) == F) dir.create(path2file)

    #### 4.3.2.1.1 create Low AIC model prediction on a specific path
    print(mod.pred.nms[2])
    mod.preds <- raster::addLayer(mod.preds, raster::writeRaster(mod.avg.i[[1]],
                                                 filename = filename,
                                                 format = formt, overwrite = T) )
    names(mod.preds)[raster::nlayers(mod.preds)] <- mod.pred.nms[2]
  }



  if(length(args.all) > length(args.aicc)){

    ##### TODO # run predictions for each xsel.mdls$sel.cri [i > length(args.aicc)]
    ## TODO # find sel.cri and run prediction; eliminate 'for', 'if's, and 'else's

    for(i in (length(args.aicc)+1):length(args.all)){
      path2file <-paste(path.mdls, outpt, mod.nms[i], sep='/')
      filename <- paste(path2file, paste0(mod.nms[i], pred.nm, ".grd"), sep='/')
      if(dir.exists(path2file) == F) dir.create(path2file)
      # grep("Mod_Mean.ORmin", xsel.mdls$sel.cri)
      # m <-  grep(mod.nms[i], xsel.mdls$sel.cri)
      print(mod.nms[i])
      mod.preds <- raster::addLayer(mod.preds, dismo::predict(mxnt.mdls[[i]], a_proj, args=pred.args, progress = 'text',
                                               file = filename,
                                               format = formt, overwrite = T) )
      names(mod.preds)[raster::nlayers(mod.preds)] <- mod.nms[i]
    }
  }

  mxnt.c.mdls <- append(mxnt.c.mdls, stats::setNames(list(mod.preds) , paste0("mxnt.preds", pred.nm)))
  return(mxnt.c.mdls)
}

#### for several species
#' Projecting calibrated MaxEnt models for several species
#'
#' This function will read an object returned by "mxnt.cpb", read the calibrated models and project into
#' new areas/climatic scenarios. These new projections will be returned together with (appended to)
#' each element (species) the original object
#' @param mxnt.c.mdls.lst A list of objects returned by "mxnt.cp", containing calibrated models.
#' @inheritParams mxnt.p
#' @inheritParams mxnt.cpb
#' @return A list of objects returned from function "mxnt.p"
#' @examples
#' mxnt.mdls.preds <- mxnt.pb(mxnt.c.mdls.lst = mxnt.mdls.preds.lst[1],
#                                     pred.nm ="fut", a_proj.lst = areas_projection)
#' @export
mxnt.pb <- function(mxnt.c.mdls.lst, pred.nm="fut", a_proj.lst, formt = "raster"){ #, #, ENMeval_occ_results.lst, occ_b_env.lst, occ_locs.lst,
                              # pred.args = c("outputformat=cloglog", "doclamp=true", "pictures=true"),
                              # wAICsum=0.99, randomseed=F, responsecurves=T, arg1='noaddsamplestobackground', arg2='noautofeature'){ #wAICsum=0.99,

  # path.res <- "4_ENMeval_results"
  # if(dir.exists(path.res)==F) dir.create(path.res)
  # path.mdls <- paste(path.res, paste0("Mdls_", names(mxnt.c.mdls.lst)), sep="/")
  mxnt.preds.lst <- vector("list", length(mxnt.c.mdls.lst))
  names(mxnt.preds.lst) <- names(mxnt.c.mdls.lst)
  for(i in 1:length(mxnt.preds.lst)){
    # if(dir.exists(path.mdls[i])==F) dir.create(path.mdls[i])
    cat(c(names(mxnt.c.mdls.lst)[i], "\n"))
    # print(paste(names(mxnt.c.mdls.lst)[i]))
    # compute final models and predictions
    mxnt.preds.lst[[i]] <- mxnt.p(mxnt.c.mdls = mxnt.c.mdls.lst[[i]], #ENMeval_occ_results = ENMeval_occ_results.lst[[i]],
                                  sp.nm = names(mxnt.c.mdls.lst)[i], pred.nm = pred.nm, a_proj = a_proj.lst[[i]],
                                       # occ_b_env = occ_b_env.lst[[i]], occ_locs = occ_locs.lst[[i]],
                                  formt = formt) # , #pred.args = pred.args,
                                       # wAICsum = wAICsum,
                                       # randomseed = randomseed, responsecurves = responsecurves, arg1 = arg1, arg2 = arg2)
    # mxnt.c.mdls.lst[[i]]$pred.args <- pred.args
  }
  return(mxnt.preds.lst)
}



#' Projecting calibrated MaxEnt models for several species onto multiple environmental scenarios
#'
#' This function will read an object returned by "mxnt.cpb", read the calibrated models and project into
#' several environmental (areas/climatic) scenarios (specified in a_proj.lst). These new projections will be returned together with (appended to)
#' each element (species) the original object.
#' @inheritParams mxnt.pb
#' @return A list of objects returned from function "mxnt.p", containing the new (multiple) projections for each element (species) of the list
#' @examples
#' mxnt.mdls.preds.pf <- mxnt.pb.Mscn(mxnt.mdls.preds.lst, a_proj.lst = area_projection.pf)
#' @export
mxnt.pb.mscn <- function(mxnt.c.mdls.lst, a_proj.lst, formt = "raster"){ #, # cores=2, #, pred.nm="fut", ENMeval_occ_results.lst, occ_b_env.lst, occ_locs.lst,
                                   # pred.args = c("outputformat=cloglog", "doclamp=true", "pictures=true"),
                                   # wAICsum=0.99, randomseed=F, responsecurves=T, arg1='noaddsamplestobackground', arg2='noautofeature'){ #wAICsum=0.99,

  # path.res <- "4_ENMeval_results"
  # if(dir.exists(path.res)==F) dir.create(path.res)
  # path.mdls <- paste(path.res, paste0("Mdls_", names(mxnt.c.mdls.lst)), sep="/")

  for(i in seq_along(mxnt.c.mdls.lst)){ # loop for each species
    # if(dir.exists(path.mdls[i])==F) dir.create(path.mdls[i])
    cat(c(names(mxnt.c.mdls.lst)[i], "\n"))
    mxnt.preds.spi <- vector("list", length(a_proj.lst[[i]]))
    # compute final models and predictions
    for(j in seq_along(mxnt.preds.spi)){
      cat(c("\n", paste0("mxnt.pred.", names(a_proj.lst[[i]])[j]), "\n",
            "projection ", j, " of ", length(mxnt.preds.spi), "\n"))
      mxnt.preds.spi[j] <- mxnt.p(mxnt.c.mdls = mxnt.c.mdls.lst[[i]],
                                  sp.nm = names(mxnt.c.mdls.lst)[i], pred.nm = names(a_proj.lst[[i]])[j],
                                  a_proj = a_proj.lst[[i]][[j]],
                                  formt = formt)[length(mxnt.c.mdls.lst[[i]]) + 1] #, # pred.args = pred.args,
                                       # wAICsum = wAICsum,
                                       # randomseed = randomseed, responsecurves = responsecurves,
                                       # arg1 = arg1, arg2 = arg2)[length(mxnt.c.mdls.lst[[i]]) + 1]
      names(mxnt.preds.spi)[j] <- paste0("mxnt.pred.", names(a_proj.lst[[i]])[j])
    }
    mxnt.c.mdls.lst[[i]] <- append(mxnt.c.mdls.lst[[i]], mxnt.preds.spi)
  }
  return(mxnt.c.mdls.lst)
}


