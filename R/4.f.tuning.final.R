## 4.2 Generating final models for occ
#' Creating arguments for selected models
#'
#' This function will read an object of class ENMevaluation (See ?ENMeval::ENMevaluate for details) and
#' return the necessary arguments for final model calibration and predictions.
#' @param x object of class ENMevaluation
#' @param wAICsum cumulative sum of top ranked models for which arguments will be created
#' @param save should save args only ("A"), selected models only ("M") or both ("B")?
#' @param randomseed logical. Args to be passed to dismo::maxent. See ?dismo::maxent and the MaxEnt help for more information.
#' @param responsecurves logical. Args to be passed to dismo::maxent. See ?dismo::maxent and the MaxEnt help for more information.
#' @param arg1 charater. Args to be passed to dismo::maxent. See ?dismo::maxent and the MaxEnt help for more information.
#' @param arg2 charater. Args to be passed to dismo::maxent. See ?dismo::maxent and the MaxEnt help for more information.
#' @examples
#' ENMeval_res.lst <- ENMevaluate.batch(occ_locs, occ_b_env, parallel = T , numCores = 7)
#' f.args(ENMeval_res.lst[[1]]@results)
#' @return A vector of args (if save="A"), data.frame of selected models (if save="M") or
#' a list with both, args and selected models, (if save="B")
#' @export
f.args <- function(x, wAICsum=0.99, save = "B", randomseed=F, responsecurves=T, arg1='noaddsamplestobackground', arg2= 'noautofeature'){ # , seq=T
  # AICc
  x <- x[order(x$delta.AICc),]
  if(is.null(x$rankAICc)) {x$rankAICc <- 1:nrow(x)} # if something goes wrong with function f.mxnt.mdl.pred, check this line
  # seq
  x.m <- x[order(x$Mean.ORmin, -x$Mean.AUC),][1,]
  # seqOr10
  x.10 <- x[order(x$Mean.OR10, -x$Mean.AUC),][1,]
  #AUCmin
  x.AUCmin <- x[order(-x$Mean.AUC, x$Mean.ORmin),][1,]
  # AUC10
  x.AUC10 <- x[order(-x$Mean.AUC, x$Mean.OR10),][1,]

  wsum <- 1:(which(cumsum(x$w.AIC) >= wAICsum)[1])
  cat("\n", paste(length(wsum)), "of",nrow(x), "models selected using AICc")# from a total of", "models")
  cat("\n", "Total AIC weight (sum of Ws) of selected models is", round(sum(x$w.AIC[wsum]), 4), "of 1")
  x <- x[wsum,]
  # if(seq==T){
  f <- factor(c(x$features, x.m$features, x.10$features, x.AUCmin$features, x.AUC10$features),
              levels=1:nlevels(x$features), labels=levels(x$features))

  beta <- c(x$rm, x.m$rm, x.10$rm, x.AUCmin$rm, x.AUC10$rm)
  cat("\n", "arguments used for building models", "\n")
  x.mdls <- data.frame(rbind(x, x.m, x.10, x.AUCmin, x.AUC10))
  x.mdls$sel.cri <- paste0("Mod_",c(paste0("AICc_", 1:length(wsum)),
                                    "Mean.ORmin", "Mean.OR10", "Mean.AUC_min", "Mean.AUC_10"))

  print(data.frame(features=f, beta, row.names = c(paste("AICc", 1:length(wsum)),
                                                   "Mean.ORmin", "Mean.OR10", "Mean.AUC_min", "Mean.AUC_10")))
  cat("\n")
  args <- paste(paste0(arg1), paste0(arg2),
                ifelse(grepl("H", f), paste("hinge"), paste("nohinge")),
                ifelse(grepl("L", f), paste("linear"), paste("nolinear")),
                ifelse(grepl("Q", f), paste("quadratic"), paste("noquadratic")),
                ifelse(grepl("P", f), paste("product"), paste("noproduct")),
                ifelse(grepl("T", f), paste("threshold"), paste("nothreshold")),
                paste0("betamultiplier=", beta),
                paste0("responsecurves=", ifelse(responsecurves==T, "true", "false")),
                paste0("randomseed=", ifelse(randomseed==T, "true", "false")),
                sep = ",")

  if(save == "A"){
    return(c(strsplit(args, ",")))
  } else if(save == "M"){
    return(x.mdls)
  } else if (save == "B"){
    return(list(c(strsplit(args, ",")),x.mdls))
  }
}

#### 4.3 Run top corresponding models and save predictions
#### 4.3.1 save maxent best models and predictions for each model
#' Calibrating and predicting selected models
#'
#' This function will read an object of class ENMevaluation (See ?ENMeval::ENMevaluate for details) and
#' return selected maxent model calibrations and predictions.
#' @param x Slot "results" of object of class ENMevaluation
#' @param sp.nm Species name. Used to name the output folder
#' @param area_projection A Raster* object or a data.frame where models will be projected. Argument 'x' of dismo::predict
#' @param occ_b_env Predictors. Used in model calibration. Argument 'x' of dismo::maxent. Raster* object or SpatialGridDataFrame, containing grids with
#' predictor variables. These will be used to extract values from for the point locations. Can
#' also be a data.frame, in which case each column should be a predictor variable and each row
#' a presence or background record.
#' @param occ_locs Occurrence data. Argument 'p' of dismo::maxent. This can be a data.frame, matrix,
#' SpatialPoints* object, or a vector. If p is a data.frame or matrix it represents a set of point
#' locations; and it must have two columns with the first being the x-coordinate (longitude) and
#' the second the y-coordinate (latitude). Coordinates can also be specified with a SpatialPoints*
#' object
#' If x is a data.frame, p should be a vector with a length equal to nrow(x) and contain 0
#' (background) and 1 (presence) values, to indicate which records (rows) in data.frame x are
#' presence records, and which are background records
#' @param pred.args Charater. Argument 'args' of dismo::maxent. Additional argument that can be
#' passed to MaxEnt. See the MaxEnt help for more information. The R maxent function only uses the
#' arguments relevant to model fitting. There is no point in using args='outputformat=raw' when
#' *fitting* the model; but you can use arguments relevant for *prediction* when using the predict
#' function. Some other arguments do not apply at all to the R implementation. An example is
#' 'outputfiletype', because the 'predict' function has its own 'filename' argument for that.
#' @inheritParams f.args
#' @return A list containing the models ('selected.mdls') used for model calibration and prediction,
#' calibrated maxent models ('mxnt.mdls'), arguments used for prediction/calibration ('pred.args'), and
#' a raster stack containing model projections ('mxnt.preds'), where each layer is a projection based on
#' a specific model selection criteria (i.e. AvgAICc, LowAICc, Mean.ORmin, Mean.OR10, Mean.AUC_min, Mean.AUC_10)
#' @export
f.mxnt.mdl.pred <- function(x, sp.nm, area_projection, occ_b_env, occ_locs, formt = "raster",
                            pred.args = c("outputformat=cloglog", "doclamp=true", "pictures=true"),
                            wAICsum=0.99, randomseed=F, responsecurves=T, arg1='noaddsamplestobackground', arg2='noautofeature'){

  # {library(xlsx)
  #   library(dismo)
  #   library(rJava)
  #   library(raster)}

  path.res <- "4_ENMeval_results"
  if(dir.exists(path.res)==F) dir.create(path.res)
  path.mdls <- paste(path.res, paste0("Mdls_", sp.nm), sep="/")
  if(dir.exists(path.mdls)==F) dir.create(path.mdls)

  mdl.arg <- f.args(x, wAICsum=wAICsum, randomseed=randomseed, responsecurves=responsecurves, arg1=arg1, arg2=arg2)
  xsel.mdls <- mdl.arg[[2]]

  args.all <- mdl.arg[[1]]
  args.aicc <- args.all[grep("Mod_AIC", xsel.mdls$sel.cri)]

  # exportar planilha de resultados
  # write.xlsx(xsel.mdls, paste0(path.mdls,"/sel.mdls.xlsx"))
  xlsx::write.xlsx(xsel.mdls, paste0(path.mdls,"/sel.mdls.", gsub("4_ENMeval_results/Mdls_", "", path.mdls), ".xlsx"))
  res.tbl <- xsel.mdls[,c("sel.cri", "features","rm","AICc", "w.AIC", "nparam", "rankAICc", "Mean.OR10", "Mean.ORmin", "Mean.AUC")]
  colnames(res.tbl) <- c("Optimality criteria", "FC", "RM", "AICc", "wAICc", "NP", "Rank", "OR10", "ORLPT", "AUC")
  xlsx::write.xlsx(res.tbl, paste0(path.mdls,"/sel.mdls.smmr.", gsub("4_ENMeval_results/Mdls_", "", path.mdls), ".xlsx"))


  mod.nms <- paste(xsel.mdls[,"sel.cri"]) # paste0("Mod_", c(1:length(args.aicc), "Mean.ORmin", "Mean.OR10", "Mean.AUC_min", "Mean.AUC_10"))
  mod.pred.nms <- c("Mod_AvgAICc", "Mod_LowAICc", mod.nms[(length(args.aicc)+1):length(args.all)])
  mod.preds <- raster::stack() #vector("list", length(mod.pred.nms))

  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

  if(dir.exists(paste(path.mdls, outpt, sep='/'))==F) dir.create(paste(path.mdls, outpt, sep='/'))

  mxnt.mdls <- vector("list", length(args.all))

  #### AIC AVG model
  {
    avg.m.path <- paste(path.mdls, outpt, mod.pred.nms[1], sep='/') # paste0("4_ENMeval_results/selected_models/cloglog/", mod.pred.nms[2])
    if(dir.exists(avg.m.path)==F) dir.create(avg.m.path)

    ##### list of models to average
    mod.avg.i <- vector("list", length(args.aicc))
    # filename <- paste(avg.m.path, mod.nms, paste0(mod.nms, ".grd"), sep='/')
    lapply(seq_along(args.aicc), function(i) {
      path2file <- paste(avg.m.path, mod.nms[i], sep='/')
      filename <- paste(path2file, paste0(mod.nms[i], ".grd"), sep='/')
      # maxent models
      set.seed(1)
      mxnt.mdls[[i]] <<- dismo::maxent(occ_b_env, occ_locs, path=path2file, args=args.all[[i]]) # final model fitting/calibration

      mod.avg.i[[i]] <<- dismo::predict(mxnt.mdls[[i]], area_projection, args=pred.args, progress='text',
                                 file = filename, format = formt, overwrite=T)
    }) #) ## fecha for or lapply

    mod.preds <- raster::stack() #vector("list", length(mod.pred.nms))

    #### 4.3.2.1.2 create model averaged prediction (models*weights, according to model selection)
    # create vector of model weights
    wv <- xsel.mdls[order(xsel.mdls$delta.AICc),"w.AIC"][seq_along(args.aicc)]

    ### stack prediction rasters (to create Average Model prediction)
    path2stk <- paste(avg.m.path, mod.nms[seq_along(args.aicc)], sep='/')
    filename <- paste(path2stk, paste0(mod.nms[seq_along(args.aicc)], pred.nm, ".grd"), sep='/')
    Mod_AICc.stack <- raster::stack(filename)

    # create averaged prediction map
    cat(c("\n",paste(mod.pred.nms[1]), "\n"))
    mod.preds <- raster::addLayer(mod.preds, raster::writeRaster(mask((sum(Mod_AICc.stack*wv, na.rm=T)/sum(wv)), area_projection[[1]]),
                                                 filename = paste(avg.m.path, paste0(mod.pred.nms[1], ".grd"), sep='/'),
                                                 format = formt, overwrite = T) )
    names(mod.preds)[nlayers(mod.preds)] <- mod.pred.nms[1]
  }

  #### Low AIC
  # if(i == 1) # usar if(low = T) pra escolher o low aic ou if(grep("low", Mod_pred))
  {
    path2file <- paste(path.mdls, outpt, mod.pred.nms[2], sep='/')
    filename <- paste(path2file, paste0(mod.pred.nms[2], ".grd"), sep='/')
    if(dir.exists(path2file)==F) dir.create(path2file)
    # print(mod.pred.nms[1])

    #### 4.3.2.1.1 create Low AIC model prediction on a specific path
    cat(c(paste(mod.pred.nms[2]), "\n"))
    mod.preds <- raster::addLayer(mod.preds, raster::writeRaster(mod.avg.i[[1]],
                                                 filename = filename,
                                                 format = formt, overwrite = T) )
    names(mod.preds)[nlayers(mod.preds)] <- mod.pred.nms[2]
  }


  #### other models
  if(length(args.all) > length(args.aicc)){

    for(i in (length(args.aicc)+1):length(args.all)){
      path2file <- paste(path.mdls, outpt, mod.nms[i], sep='/')
      filename  <- paste(path2file, paste0(mod.nms[i], ".grd"), sep='/')
      if(dir.exists(path2file)==F) dir.create(path2file)
      # maxent models
      set.seed(1)
      mxnt.mdls[[i]] <- dismo::maxent(occ_b_env, occ_locs, path=path2file, args=args.all[[i]])

      # grep("Mod_Mean.ORmin", xsel.mdls$sel.cri)
      # m <-  grep(mod.nms[i], xsel.mdls$sel.cri)
      cat(c(paste(mod.nms[i]), "\n"))
      mod.preds <- raster::addLayer(mod.preds, dismo::predict(mxnt.mdls[[i]], area_projection, args=pred.args, progress='text',
                                               file = filename,
                                               format = formt, overwrite = T) )
      names(mod.preds)[nlayers(mod.preds)] <- mod.nms[i]
    }
  }
  return(list(selected.mdls = xsel.mdls, mxnt.mdls=mxnt.mdls, mxnt.args = args.all, pred.args = pred.args, mxnt.preds = mod.preds))
}

#' Calibrating and predicting selected models for several species
#'
#' This function will read a list of objects of class ENMevaluation (See ?ENMeval::ENMevaluate for details) and
#' return selected maxent model calibrations and predictions. Each element on the list is usually a species.
#' area_proj.lst, occ_b_env.lst, occ_locs.lst are lists with occurence data, projection and calibration/predictor data.
#' Species in these lists must all be in the same order of species in ENMeval_res.
#' @param ENMeval_res List of objects of class ENMevaluation
#' @param area_proj.lst List of projection areas. See argument "area_projection" in f.mxnt.mdl.pred.
#' @param occ_b_env.lst List of predictor areas.
#' @param occ_locs.lst List of occurence data. See argument "occ_locs" in f.mxnt.mdl.pred.
#' @param formt Character. Output file type. Argument 'format' of raster::writeRaster
#' @inheritParams f.mxnt.mdl.pred
#' @return A list of objects returned from f.mxnt.mdl.pred
#' @examples
#' mxnt.mdls.preds.lst <- f.mxnt.mdl.pred.batch(ENMeval_res=ENMeval_res.lst,  area_proj.lst=area_projection,
#' occ_b_env.lst=occ_b_env, occ_locs.lst=occ_locs, wAICsum=0.99)
#' mxnt.mdls.preds.lst[[1]][[1]] # models [ENM]evaluated and selected using sum of wAICc
#' mxnt.mdls.preds.lst[[1]][[2]] # MaxEnt models
#' mxnt.mdls.preds.lst[[1]][[3]] # used prediction arguments
#' plot(mxnt.mdls.preds.lst[[1]][[4]]) # MaxEnt predictions, based on the model selection criteria
#' @export
f.mxnt.mdl.pred.batch <- function(ENMeval_res, area_proj.lst, occ_b_env.lst, occ_locs.lst, formt = "raster",
                                  pred.args = c("outputformat=cloglog", "doclamp=true", "pictures=true"),
                                  wAICsum=0.99, randomseed=F, responsecurves=T, arg1='noaddsamplestobackground', arg2='noautofeature'){

  # path.res <- "4_ENMeval_results"
  # if(dir.exists(path.res)==F) dir.create(path.res)
  # path.mdls <- paste(path.res, paste0("Mdls_", names(ENMeval_res)), sep="/")
  mxnt.mdls.preds.lst <- vector("list", length(ENMeval_res))
  names(mxnt.mdls.preds.lst) <- names(ENMeval_res)
  for(i in seq_along(ENMeval_res)){
    ## TODO - check this, decide if keep other fields before or remove only here (in which use loop to get)
    ENMeval_res[[i]] <- ENMeval_res[[i]]@results
    cat(c(names(mxnt.mdls.preds.lst)[i], "\n"))
    # if(dir.exists(path.mdls[i])==F) dir.create(path.mdls[i])
    # compute final models and predictions
    mxnt.mdls.preds.lst[[i]] <- f.mxnt.mdl.pred(x = ENMeval_res[[i]], sp.nm = names(ENMeval_res[i]), area_projection = area_proj.lst[[i]],
                                                occ_b_env = occ_b_env.lst[[i]], occ_locs = occ_locs.lst[[i]],
                                                formt = formt, pred.args = pred.args, wAICsum = wAICsum,
                                                randomseed = randomseed, responsecurves = responsecurves, arg1 = arg1, arg2 = arg2)
    # mxnt.mdls.preds.lst[[i]]$pred.args <- pred.args
  }
  return(mxnt.mdls.preds.lst)
}




