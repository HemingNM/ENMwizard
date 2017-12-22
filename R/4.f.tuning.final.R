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
#' #' @examples
#' ENMeval.res.lst <- ENMevaluate.batch(occ.locs, occ.b.env, parallel = T , numCores = 7)
#' f.args(ENMeval.res.lst[[1]]@results)
#' @return A vector of args (if save="A"), data.frame of selected models (if save="M") or
#' a list with both, args and selected models, (if save="B")
#' @export
f.args <- function(x, wAICsum=0.99, save = "B", randomseed=F, responsecurves=T, arg1='noaddsamplestobackground', arg2= 'noautofeature'){ # , seq=T
  # AICc
  x <- x[order(x$delta.AICc),]
  if(is.null(x$rankAICc)) {x$rankAICc <- 1:nrow(x)} # if something goes wrong with function mxnt.cp, check this line
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
  x.mdls$sel.cri <- paste0("Mod.",c(paste0("AICc_", 1:length(wsum)),
                                    "Mean.ORmin", "Mean.OR10", "Mean.AUCmin", "Mean.AUC10"))

  print(data.frame(features=f, beta, row.names = c(paste("AICc", 1:length(wsum)),
                                                   "Mean.ORmin", "Mean.OR10", "Mean.AUCmin", "Mean.AUC10")))
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
# "f.mxnt.mdl.pred" renamed to "mxnt.cp"
#' Calibrating and predicting selected models
#'
#' This function will read an object of class ENMevaluation (See ?ENMeval::ENMevaluate for details) and
#' return selected maxent model calibrations and predictions.
#' @param x Slot "results" of object of class ENMevaluation
#' @param sp.nm Species name. Used to name the output folder
#' @param a.calib Predictors (cropped environmental variables) for model tuning. Used in model calibration. Argument 'x' of dismo::maxent. Raster* object or SpatialGridDataFrame, containing grids with
#' predictor variables. These will be used to extract values from for the point locations. Can
#' also be a data.frame, in which case each column should be a predictor variable and each row
#' a presence or background record.
# #' @param a.proj A Raster* object or a data.frame where models will be projected. Argument 'x' of dismo::predict
#' @param occ Occurrence data. Argument 'p' of dismo::maxent. This can be a data.frame, matrix,
#' SpatialPoints* object, or a vector. If p is a data.frame or matrix it represents a set of point
#' locations; and it must have two columns with the first being the x-coordinate (longitude) and
#' the second the y-coordinate (latitude). Coordinates can also be specified with a SpatialPoints*
#' object
#' If x is a data.frame, p should be a vector with a length equal to nrow(x) and contain 0
#' (background) and 1 (presence) values, to indicate which records (rows) in data.frame x are
#' presence records, and which are background records
#' @param formt Character. Output file type. Argument 'format' of raster::writeRaster
#' @param pred.args Charater. Argument 'args' of dismo::maxent. Additional argument that can be
#' passed to MaxEnt. See the MaxEnt help for more information. The R maxent function only uses the
#' arguments relevant to model fitting. There is no point in using args='outputformat=raw' when
#' *fitting* the model; but you can use arguments relevant for *prediction* when using the predict
#' function. Some other arguments do not apply at all to the R implementation. An example is
#' 'outputfiletype', because the 'predict' function has its own 'filename' argument for that.
#' @param numCores Number of cores to use for parallelization. If set to 1, no paralellization is performed
#' @param parallelTunning Should parallelize within species (parallelTunning=TRUE) or between species (parallelTunning=FALSE)
#' @inheritParams f.args
#' @return A list containing the models ('selected.mdls') used for model calibration and prediction,
#' calibrated maxent models ('mxnt.mdls'), arguments used for prediction/calibration ('pred.args'), and
#' a raster stack containing model projections ('mxnt.preds'), where each layer is a projection based on
#' a specific model selection criteria (i.e. AvgAICc, LowAICc, Mean.ORmin, Mean.OR10, Mean.AUCmin, Mean.AUC10)
#' @export
mxnt.cp <- function(x, sp.nm, a.calib, occ, formt = "raster", # , a.proj
                    pred.args = c("outputformat=cloglog", "doclamp=true", "pictures=true"),
                    wAICsum=0.99, randomseed=F, responsecurves=T, arg1='noaddsamplestobackground', arg2='noautofeature',
                    numCores = 1, parallelTunning = TRUE){

  path.res <- "4_ENMeval.results"
  if(dir.exists(path.res)==FALSE) dir.create(path.res)
  path.mdls <- paste(path.res, paste0("Mdls.", sp.nm), sep="/")
  if(dir.exists(path.mdls)==FALSE) dir.create(path.mdls)

  mdl.arg <- f.args(x, wAICsum=wAICsum, randomseed=randomseed, responsecurves=responsecurves, arg1=arg1, arg2=arg2)
  xsel.mdls <- mdl.arg[[2]]

  args.all <- mdl.arg[[1]]
  args.aicc <- args.all[grep("Mod.AIC", xsel.mdls$sel.cri)]

  # exportar planilha de resultados
  # write.xlsx(xsel.mdls, paste0(path.mdls,"/sel.mdls.xlsx"))
  xlsx::write.xlsx(xsel.mdls, paste0(path.mdls,"/sel.mdls.", gsub("4_ENMeval.results/Mdls.", "", path.mdls), ".xlsx"))
  res.tbl <- xsel.mdls[,c("sel.cri", "features","rm","AICc", "w.AIC", "nparam", "rankAICc", "Mean.OR10", "Mean.ORmin", "Mean.AUC")]
  colnames(res.tbl) <- c("Optimality criteria", "FC", "RM", "AICc", "wAICc", "NP", "Rank", "OR10", "ORLPT", "AUC")
  xlsx::write.xlsx(res.tbl, paste0(path.mdls,"/sel.mdls.smmr.", gsub("4_ENMeval.results/Mdls.", "", path.mdls), ".xlsx"))


  mod.nms <- paste(xsel.mdls[,"sel.cri"]) # paste0("Mod.", c(1:length(args.aicc), "Mean.ORmin", "Mean.OR10", "Mean.AUCmin", "Mean.AUC10"))
  mod.pred.nms <- c("Mod.AvgAICc", "Mod.LowAICc", mod.nms[(length(args.aicc)+1):length(args.all)])
  mod.preds <- raster::stack() #vector("list", length(mod.pred.nms))

  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

  if(dir.exists(paste(path.mdls, outpt, sep='/'))==F) dir.create(paste(path.mdls, outpt, sep='/'))

  mxnt.mdls <- vector("list", length(args.all))

  #### AIC AVG model
  {
    avg.m.path <- paste(path.mdls, outpt, mod.pred.nms[1], sep='/') # paste0("4_ENMeval.results/selected.models/cloglog/", mod.pred.nms[2])
    if(dir.exists(avg.m.path)==F) dir.create(avg.m.path)

    ##### list of models to average
    mod.avg.i <- vector("list", length(args.aicc))
    # filename <- paste(avg.m.path, mod.nms, paste0(mod.nms, ".grd"), sep='/')

    if(numCores>1 & parallelTunning){

      cl<-parallel::makeCluster(numCores)

      mxnt.mdls <- parallel::clusterApply(cl, seq_along(args.all), function(i, args.all, mod.nms, a.calib, occ) {

        if(i<=length(args.aicc)){
          path2file <- paste(getwd(), avg.m.path, mod.nms[i], sep='/')
          }else{
          path2file <- paste(path.mdls, outpt, mod.nms[i], sep='/')
          }


        filename <- paste(path2file, paste0(mod.nms[i], ".grd"), sep='/')
        # maxent models
        set.seed(1)
        resu <- dismo::maxent(a.calib, occ, path=path2file, args=args.all[[i]]) # final model fitting/calibration
        return(resu)
        # mod.avg.i[[i]] <<- dismo::predict(mxnt.mdls[[i]], a.proj, args=pred.args, progress='text',
        #                            file = filename, format = formt, overwrite=T)
      }, args.all, mod.nms, a.calib, occ) #)

      parallel::stopCluster(cl)



    }else{
      mxnt.mdls<-lapply(seq_along(args.all), function(i) {

        if(i<=length(args.aicc)){
          path2file <- paste(getwd(), avg.m.path, mod.nms[i], sep='/')
        }else{
          path2file <- paste(path.mdls, outpt, mod.nms[i], sep='/')
        }

      filename <- paste(path2file, paste0(mod.nms[i], ".grd"), sep='/')
      # maxent models
      set.seed(1)
      resu <- dismo::maxent(a.calib, occ, path=path2file, args=args.all[[i]]) # final model fitting/calibration
      return(resu)
      # mod.avg.i[[i]] <<- dismo::predict(mxnt.mdls[[i]], a.proj, args=pred.args, progress='text',
      #                            file = filename, format = formt, overwrite=T)
    }) #) ## fecha for or lapply

    } # closes else


}
  return(list(selected.mdls = xsel.mdls, mxnt.mdls=mxnt.mdls, mxnt.args = args.all, pred.args = pred.args)) #, mxnt.preds = mod.preds))
}

# "f.mxnt.mdl.pred.batch" renamed to "mxnt.cp.batch"
#' Calibrating and predicting selected models for several species
#'
#' This function will read a list of objects of class ENMevaluation (See ?ENMeval::ENMevaluate for details) and
#' return selected maxent model calibrations and predictions. Each element on the list is usually a species.
#' a.proj.l, a.calib.l, occ.l are lists with occurence data, projection and calibration/predictor data.
#' Species in these lists must all be in the same order of species in ENMeval.res.
#' @param ENMeval.res List of objects of class ENMevaluation
#' @param a.calib.l List of predictors (cropped environmental variables) for model tuning. Used in model calibration. Argument 'x' of dismo::maxent. Raster* object or SpatialGridDataFrame, containing grids with
#' predictor variables. These will be used to extract values from for the point locations. Can
#' also be a data.frame, in which case each column should be a predictor variable and each row
#' a presence or background record..
# #' @param a.proj.l List of projection areas. See argument "a.proj" in mxnt.cp.
#' @param occ.l List of occurence data. See argument "occ" in mxnt.cp.
#' @param numCores Number of cores to use for parallelization. If set to 1, no paralellization is performed
#' @inheritParams mxnt.cp
#' @return A list of objects returned from function "mxnt.cp"
#' @examples
#' mxnt.mdls.preds.lst <- mxnt.cp.batch(ENMeval.res=ENMeval.res.lst,
#' a.calib.l=occ.b.env, a.proj.l=areas.projection, occ.l=occ, wAICsum=0.99)
#' mxnt.mdls.preds.lst[[1]][[1]] # models [ENM]evaluated and selected using sum of wAICc
#' mxnt.mdls.preds.lst[[1]][[2]] # MaxEnt models
#' mxnt.mdls.preds.lst[[1]][[3]] # used prediction arguments
#' plot(mxnt.mdls.preds.lst[[1]][[4]]) # MaxEnt predictions, based on the model selection criteria
#' @export
mxnt.cp.batch <- function(ENMeval.res, a.calib.l, occ.l, formt = "raster", # , a.proj.l
                          pred.args = c("outputformat=cloglog", "doclamp=true", "pictures=true"),
                          wAICsum=0.99, randomseed=F, responsecurves=T, arg1='noaddsamplestobackground', arg2='noautofeature',
                          numCores = 1, parallelTunning = TRUE){

  # path.res <- "4_ENMeval.results"
  # if(dir.exists(path.res)==F) dir.create(path.res)
  # path.mdls <- paste(path.res, paste0("Mdls.", names(ENMeval.res)), sep="/")

  if(numCores>1 & parallelTunning==FALSE){

    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl,list("mxnt.cp","f.args"))

  mxnt.mdls.preds.lst <- parallel::clusterApply(cl, base::seq_along(ENMeval.res), function(i, ENMeval.res, a.calib.l, occ.l, formt, pred.args, wAICsum, randomseed, responsecurves, arg1, arg2, numCores, parallelTunning){
      ## TODO - check this, decide if keep other fields before or remove only here (in which use loop to get)
      ENMeval.res[[i]] <- ENMeval.res[[i]]@results
      cat(c(names(ENMeval.res[i]), "\n"))
      # if(dir.exists(path.mdls[i])==F) dir.create(path.mdls[i])
      # compute final models and predictions
     resu <- mxnt.cp(x = ENMeval.res[[i]], sp.nm = names(ENMeval.res[i]),
                                          a.calib = a.calib.l[[i]], # a.proj = a.proj.l[[i]],
                                          occ = occ.l[[i]], formt = formt,
                                          pred.args = pred.args, wAICsum = wAICsum,
                                          randomseed = randomseed, responsecurves = responsecurves, arg1 = arg1, arg2 = arg2,numCores=numCores,parallelTunning=parallelTunning)

    return(resu)
  }, ENMeval.res, a.calib.l, occ.l, formt, pred.args, wAICsum, randomseed, responsecurves, arg1, arg2, numCores, parallelTunning)

  parallel::stopCluster(cl)

  }else{

    mxnt.mdls.preds.lst <- lapply(base::seq_along(ENMeval.res), function(i, ENMeval.res, a.calib.l, occ.l, formt, pred.args, wAICsum, randomseed, responsecurves, arg1, arg2, numCores, parallelTunning){
      ## TODO - check this, decide if keep other fields before or remove only here (in which use loop to get)
      ENMeval.res[[i]] <- ENMeval.res[[i]]@results
      cat(c(names(ENMeval.res[i]), "\n"))
      # if(dir.exists(path.mdls[i])==F) dir.create(path.mdls[i])
      # compute final models and predictions
      resu <- mxnt.cp(x = ENMeval.res[[i]], sp.nm = names(ENMeval.res[i]),
                      a.calib = a.calib.l[[i]], # a.proj = a.proj.l[[i]],
                      occ = occ.l[[i]], formt = formt,
                      pred.args = pred.args, wAICsum = wAICsum,
                      randomseed = randomseed, responsecurves = responsecurves, arg1 = arg1, arg2 = arg2,numCores=numCores,parallelTunning=parallelTunning)

      return(resu)
    }, ENMeval.res, a.calib.l, occ.l, formt, pred.args, wAICsum, randomseed, responsecurves, arg1, arg2, numCores, parallelTunning)


  }

  names(mxnt.mdls.preds.lst) <- names(ENMeval.res)

  # mxnt.mdls.preds.lst <- vector("list", length(ENMeval.res))
  # names(mxnt.mdls.preds.lst) <- names(ENMeval.res)

  # for(i in base::seq_along(ENMeval.res)){
  #   ## TODO - check this, decide if keep other fields before or remove only here (in which use loop to get)
  #   ENMeval.res[[i]] <- ENMeval.res[[i]]@results
  #   cat(c(names(mxnt.mdls.preds.lst)[i], "\n"))
  #   # if(dir.exists(path.mdls[i])==F) dir.create(path.mdls[i])
  #   # compute final models and predictions
  #   mxnt.mdls.preds.lst[[i]] <- mxnt.cp(x = ENMeval.res[[i]], sp.nm = names(ENMeval.res[i]),
  #                                       a.calib = a.calib.l[[i]], # a.proj = a.proj.l[[i]],
  #                                       occ = occ.l[[i]], formt = formt,
  #                                       pred.args = pred.args, wAICsum = wAICsum,
  #                                       randomseed = randomseed, responsecurves = responsecurves, arg1 = arg1, arg2 = arg2,numCores=numCores,parallelTunning=parallelTunning)
  #   # mxnt.mdls.preds.lst[[i]]$pred.args <- pred.args
  # }

  return(mxnt.mdls.preds.lst)
}




