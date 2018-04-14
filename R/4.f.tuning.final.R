## 4.2 Generating final models for occ
#' Creating arguments for selected models
#'
#' This function will read an object of class ENMevaluation (See ?ENMeval::ENMevaluate for details) and
#' return the necessary arguments for final model calibration and predictions.
#'
#' @param x Slot "results" of object of class ENMevaluation
#' @param mSel character vector. Which criteria to use when selecting model(s). Currently implemented:
#' "AICavg", "LowAIC", "OR", "AUC"
#' @param wAICsum cumulative sum of top ranked models for which arguments will be created
#' @param save should save args only ("A"), selected models only ("M") or both ("B")?
#' @inheritParams dismo::maxent
#' @param randomseed logical. Args to be passed to dismo::maxent. See ?dismo::maxent and the MaxEnt help for more information.
#' @param responsecurves logical. Args to be passed to dismo::maxent. See ?dismo::maxent and the MaxEnt help for more information.
#' @param arg1 charater. Args to be passed to dismo::maxent. See ?dismo::maxent and the MaxEnt help for more information.
#' @param arg2 charater. Args to be passed to dismo::maxent. See ?dismo::maxent and the MaxEnt help for more information.
#' @seealso \code{\link{mxnt.c}}, \code{\link{mxnt.c.batch}}, \code{\link[dismo]{maxent}}, \code{\link[ENMeval]{ENMevaluate}},
#' \code{\link{mxnt.p}}, \code{\link{mxnt.p.batch.mscn}}
#' #' @examples
#' ENMeval.res.lst <- ENMevaluate.batch(occ.locs, occ.b.env, parallel = T , numCores = 7)
#' f.args(ENMeval.res.lst[[1]]@results)
#' @return A vector of args (if save="A"), data.frame of selected models (if save="M") or
#' a list with both, args and selected models, (if save="B")
#' @export
f.args <- function(x, mSel=c("AICavg", "LowAIC", "OR", "AUC"), wAICsum=0.99, save="B", randomseed=FALSE, responsecurves=TRUE, arg1='noaddsamplestobackground', arg2='noautofeature'){ # , seq=TRUE

  x.a.i <- order(x$delta.AICc)
  x <- x[x.a.i,]
  if(is.null(x$rankAICc)) {x$rankAICc <- 1:nrow(x)} # if something goes wrong with function mxnt.c, check this line
  x$ID <- as.numeric(rownames(x))

  x$opt.crit <- ""
  # AICcAvg
  if("AICavg" %in% mSel){
    if(any(cumsum(x$w.AIC) >= wAICsum)){
      wsum <- 1:(which(cumsum(x$w.AIC) >= wAICsum)[1])
    } else {
      wsum <- 1:length(x$w.AIC)
    }

    x.a.i <- x.a.i[wsum]
    x$opt.crit[wsum] <- paste0(x$opt.crit[wsum], paste0("AICc_", wsum))

    cat("\n", paste(length(wsum)), "of", nrow(x), "models selected using AICc")# from a total of", "models")
    cat("\n", "Total AIC weight (sum of Ws) of selected models is", round(sum(x$w.AIC[wsum]), 4), "of 1")

  } else {
    x.a.i <- NULL
  }
  x <- x[order(x$ID),]

  # LowAIC
  if("LowAIC" %in% mSel){
    x.la.i <- order(x$delta.AICc)[1]
    x$opt.crit[x.la.i] <- sub("^\\.", "", paste(x$opt.crit[x.la.i], "LowAICc", sep = "."))
  } else {
    x.la.i <- NULL
  }

  # OR
  if("OR" %in% mSel){
    # seq
    xORm.i <- order(x$Mean.ORmin, -x$Mean.AUC)[1]
    x$opt.crit[xORm.i] <- sub("^\\.", "", paste(x$opt.crit[xORm.i], "ORmin", sep = "."))
    # seqOr10
    xOR10.i <- order(x$Mean.OR10, -x$Mean.AUC)[1]
    x$opt.crit[xOR10.i] <- sub("^\\.", "", paste(x$opt.crit[xOR10.i], "OR10", sep = "."))
  } else {
    xORm.i <- NULL
    xOR10.i <- NULL
  }

  # AUC
  if("AUC" %in% mSel){
    #AUCmin
    xAUCmin.i <- order(-x$Mean.AUC, x$Mean.ORmin)[1]
    x$opt.crit[xAUCmin.i] <- sub("^\\.", "", paste(x$opt.crit[xAUCmin.i], "AUCmin", sep = "."))
    # AUC10
    xAUC10.i <- order(-x$Mean.AUC, x$Mean.OR10)[1]
    x$opt.crit[xAUC10.i] <- sub("^\\.", "", paste(x$opt.crit[xAUC10.i], "AUC10", sep = "."))
  } else {
    xAUCmin.i <- NULL
    xAUC10.i <- NULL
  }


  xsel.mdls <- x#[(unique(c(x.a.i, x.la.i, xORm.i, xOR10.i, xAUCmin.i, xAUC10.i))),]
  xsel.mdls$ID <- NULL

  f <- factor(xsel.mdls$features)
  beta <- c(xsel.mdls$rm)
  cat("\n", "arguments used for building models", "\n")
  print(data.frame(optimality.criteria = xsel.mdls$opt.crit, features=xsel.mdls$features, beta=xsel.mdls$rm))

  cat("\n")
  args <- paste(paste0(arg1), paste0(arg2),
                ifelse(grepl("H", f), paste("hinge"), paste("nohinge")),
                ifelse(grepl("L", f), paste("linear"), paste("nolinear")),
                ifelse(grepl("Q", f), paste("quadratic"), paste("noquadratic")),
                ifelse(grepl("P", f), paste("product"), paste("noproduct")),
                ifelse(grepl("T", f), paste("threshold"), paste("nothreshold")),
                paste0("betamultiplier=", beta),
                paste0("responsecurves=", ifelse(responsecurves==TRUE, "true", "false")),
                paste0("randomseed=", ifelse(randomseed==TRUE, "true", "false")),
                sep = ",")

  if(save == "A"){
    return(c(strsplit(args, ",")))
  } else if(save == "M"){
    return(xsel.mdls)
  } else if (save == "B"){
    return(list(args=c(strsplit(args, ",")), mdls=xsel.mdls))
  }
}

#### 4.3 Run top corresponding models and save predictions
#### 4.3.1 save maxent best models and predictions for each model
# "f.mxnt.mdl.pred" renamed to "mxnt.c"
#' Calibrating selected models
#'
#' This function will read an object of class ENMevaluation (See ?ENMeval::ENMevaluate for details) and
#' return selected maxent models calibrated.
#'
#' @param ENMeval.o Object of class ENMevaluation
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
#' @param use.ENMeval.bgpts Logical. Use background points from ENMeval or sample new ones?
#' @inheritParams f.args
#' @inheritParams dismo::maxent
#' @seealso \code{\link{f.args}}, \code{\link{mxnt.c.batch}}, \code{\link{ENMevaluate.batch}}, \code{\link[ENMeval]{ENMevaluate}},
#' \code{\link[dismo]{maxent}}, \code{\link{mxnt.p}}, \code{\link{mxnt.p.batch.mscn}}
#' @return A 'mcm' (mxnt.c.mdls, Maxent Calibrated Models). A list containing the models ('selected.mdls') used for model calibration,
#' calibrated maxent models ('mxnt.mdls'), and arguments used for calibration ('pred.args').
#' @export
mxnt.c <- function(ENMeval.o, sp.nm, a.calib, occ = NULL, use.ENMeval.bgpts = TRUE, nbg=10000, formt = "raster", # , a.proj
                   pred.args = c("outputformat=cloglog", "doclamp=true", "pictures=true"),
                   mSel=c("AICavg", "LowAIC", "OR", "AUC"), wAICsum = 0.99, randomseed = FALSE,
                   responsecurves = TRUE, arg1 = 'noaddsamplestobackground', arg2 = 'noautofeature',
                   numCores = 1, parallelTunning = TRUE){

  path.res <- "3_out.MaxEnt"
  if(dir.exists(path.res)==FALSE) dir.create(path.res)
  path.mdls <- paste(path.res, paste0("Mdls.", sp.nm), sep="/")
  if(dir.exists(path.mdls)==FALSE) dir.create(path.mdls)

  if(is.null(occ)){ occ <- ENMeval.o@occ.pts }

  if(use.ENMeval.bgpts){
    a <- ENMeval.o@bg.pts
  } else {
    a <- dismo::randomPoints(a.calib, nbg, occ)
  }

  ENMeval.r <- ENMeval.o@results
  # cat(c(names(ENMeval.r[i]), "\n"))

  mdl.arg <- f.args(x=ENMeval.r, mSel=mSel, wAICsum=wAICsum, randomseed=randomseed, responsecurves=responsecurves, arg1=arg1, arg2=arg2)
  xsel.mdls <- mdl.arg[[2]]
  ENMeval.r <- cbind(ENMeval.r, opt.crit=xsel.mdls$opt.crit)

  xsel.mdls <- xsel.mdls[xsel.mdls$opt.crit!="",]


  args.all <- mdl.arg[[1]]
  args.aicc <- grep("AIC", xsel.mdls$opt.crit)

  # exportar planilha de resultados
  # write.xlsx(xsel.mdls, paste0(path.mdls,"/sel.mdls.xlsx"))
  # xlsx::write.xlsx(xsel.mdls, paste0(path.mdls,"/sel.mdls.", gsub("3_out.MaxEnt/Mdls.", "", path.mdls), ".xlsx"))
  utils::write.csv(ENMeval.r, paste0(path.mdls,"/sel.mdls.", gsub("3_out.MaxEnt/Mdls.", "", path.mdls), ".csv"))
  res.tbl <- ENMeval.r[,c("opt.crit", "features","rm","AICc", "w.AIC", "nparam", "rankAICc", "Mean.OR10", "Mean.ORmin", "Mean.AUC")]
  colnames(res.tbl) <- c("Optimality criteria", "FC", "RM", "AICc", "wAICc", "NP", "Rank", "OR10", "ORLPT", "AUC")
  # xlsx::write.xlsx(res.tbl, paste0(path.mdls,"/sel.mdls.smmr.", gsub("3_out.MaxEnt/Mdls.", "", path.mdls), ".xlsx"))
  utils::write.csv(res.tbl, paste0(path.mdls,"/sel.mdls.smmr.", gsub("3_out.MaxEnt/Mdls.", "", path.mdls), ".csv"))

  mod.nms <- paste0("Mod.", xsel.mdls[, "opt.crit"]) # mod.nms <- paste(xsel.mdls[, "opt.crit"]) # paste0("Mod.", c(1:length(args.aicc), "Mean.ORmin", "Mean.OR10", "Mean.AUCmin", "Mean.AUC10"))
  # mod.pred.nms <- c("Mod.AvgAICc", "Mod.LowAICc", mod.nms[(length(args.aicc)+1):length(args.all)])
  # mod.pred.nms <- c(if(length(args.aicc)>0){
  #   c("Mod.AvgAICc", "Mod.LowAICc")
  # }, mod.nms[(length(args.aicc)+1):length(args.all)])
  mod.pred.nms <- c(if(length(args.aicc)>1){"Mod.AvgAICc"}, # if(length(grep("LowAIC", xsel.mdls$opt.crit))>0){"Mod.LowAICc"},
                    paste0("Mod.", mod.nms[1:length(args.all)]))

  mod.preds <- raster::stack() #vector("list", length(mod.pred.nms))

  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

  if(dir.exists(paste(path.mdls, outpt, sep='/'))==FALSE) dir.create(paste(path.mdls, outpt, sep='/'))

  mxnt.mdls <- vector("list", length(args.all))

  #### AIC AVG model
  {
    # if(length(args.aicc)>0) {
    #   avg.m.path <- paste(path.mdls, outpt, mod.pred.nms[1], sep='/') # paste0("3_out.MaxEnt/selected.models/cloglog/", mod.pred.nms[2])
    #   if(dir.exists(avg.m.path)==FALSE) dir.create(avg.m.path)
    # }

    ##### list of models to average
    # mod.avg.i <- vector("list", length(args.aicc))
    # filename <- paste(avg.m.path, mod.nms, paste0(mod.nms, ".grd"), sep='/')

    if(numCores>1 & parallelTunning){

      cl <- parallel::makeCluster(numCores)

      mxnt.mdls <- parallel::clusterApply(cl, seq_along(args.all), function(i, args.all, mod.nms, a.calib, occ, a) {

        # if(i<=length(args.aicc)){
        #   path2file <- paste(avg.m.path, mod.nms[i], sep='/')
        # }else{
          path2file <- paste(path.mdls, outpt, mod.nms[i], sep='/')
        # }

        filename <- paste(path2file, paste0(mod.nms[i], ".grd"), sep='/')
        # maxent models
        set.seed(1)
        resu <- dismo::maxent(a.calib, occ, a, path=path2file, args=args.all[[i]]) # final model fitting/calibration
        return(resu)
        # mod.avg.i[[i]] <<- dismo::predict(mxnt.mdls[[i]], a.proj, args=pred.args, progress='text',
        #                            file = filename, format = formt, overwrite=TRUE)
      }, args.all, mod.nms, a.calib, occ, a) #)

      parallel::stopCluster(cl)



    } else {
      mxnt.mdls <- lapply(seq_along(args.all), function(i, args.all, mod.nms, a.calib, occ, a) {

        # if(i<=length(args.aicc)){
        #   path2file <- paste(avg.m.path, mod.nms[i], sep='/')
        # }else{
          path2file <- paste(path.mdls, outpt, mod.nms[i], sep='/')
        # }

        filename <- paste(path2file, paste0(mod.nms[i], ".grd"), sep='/')
        # maxent models
        set.seed(1)
        resu <- dismo::maxent(a.calib, occ, a, path=path2file, args=args.all[[i]]) # final model fitting/calibration
        return(resu)
        # mod.avg.i[[i]] <<- dismo::predict(mxnt.mdls[[i]], a.proj, args=pred.args, progress='text',
        #                            file = filename, format = formt, overwrite=TRUE)
      }, args.all, mod.nms, a.calib, occ, a) #) ## fecha for or lapply

    } # closes else


  }
  return(list(ENMeval.results = ENMeval.r, mxnt.mdls = mxnt.mdls, selected.mdls = xsel.mdls,
              occ.pts = occ, bg.pts = a,
              mxnt.args = args.all, pred.args = pred.args)) #, mxnt.preds = mod.preds))
}

# "f.mxnt.mdl.pred.batch" renamed to "mxnt.c.batch"
#' Calibrating and predicting selected models for several species
#'
#' This function will read a list of objects of class ENMevaluation (See ?ENMeval::ENMevaluate for details) and
#' return selected maxent model calibrations and predictions. Each element on the list is usually a species.
#' a.proj.l, a.calib.l, occ.l are lists with occurence data, projection and calibration/predictor data.
#' Species in these lists must all be in the same order of species in ENMeval.o.
#' @param ENMeval.o.l List of objects of class ENMevaluation
#' @param a.calib.l List of predictors (cropped environmental variables) for model tuning. Used in model calibration. Argument 'x' of dismo::maxent. Raster* object or SpatialGridDataFrame, containing grids with
#' predictor variables. These will be used to extract values from for the point locations. Can
#' also be a data.frame, in which case each column should be a predictor variable and each row
#' a presence or background record..
# #' @param a.proj.l List of projection areas. See argument "a.proj" in mxnt.c.
#' @param occ.l List of occurence data. See argument "occ" in mxnt.c.
# #' @param numCores Number of cores to use for parallelization. If set to 1, no paralellization is performed
#' @inheritParams mxnt.c
#' @inheritParams poly.c.batch
#' @seealso \code{\link{f.args}}, \code{\link{mxnt.c}}, \code{\link{ENMevaluate.batch}}, \code{\link[ENMeval]{ENMevaluate}},
#' \code{\link[dismo]{maxent}}, \code{\link{mxnt.p}}, \code{\link{mxnt.p.batch.mscn}}
#' @return A 'mcm.l' object. A list of 'mcm' (mxnt.c.mdls, Maxent Calibrated Models), returned from function "mxnt.c"
#' @examples
#' mxnt.mdls.preds.lst <- mxnt.c.batch(ENMeval.o=ENMeval.res.lst,
#' a.calib.l=occ.b.env, a.proj.l=areas.projection, occ.l=occ, wAICsum=0.99)
#' mxnt.mdls.preds.lst[[1]][[1]] # models [ENM]evaluated and selected using sum of wAICc
#' mxnt.mdls.preds.lst[[1]][[2]] # MaxEnt models
#' mxnt.mdls.preds.lst[[1]][[3]] # used prediction arguments
#' plot(mxnt.mdls.preds.lst[[1]][[4]]) # MaxEnt predictions, based on the model selection criteria
#' @export
mxnt.c.batch <- function(ENMeval.o.l, a.calib.l, occ.l = NULL, use.ENMeval.bgpts = TRUE, formt = "raster", # , a.proj.l
                         pred.args = c("outputformat=cloglog", "doclamp=true", "pictures=true"),
                         mSel=c("AICavg", "LowAIC", "OR", "AUC"), wAICsum = 0.99, randomseed = FALSE,
                         responsecurves = TRUE, arg1 = 'noaddsamplestobackground', arg2 = 'noautofeature',
                         numCores = 1, parallelTunning = TRUE){

  # path.res <- "3_out.MaxEnt"
  # if(dir.exists(path.res)==FALSE) dir.create(path.res)
  # path.mdls <- paste(path.res, paste0("Mdls.", names(ENMeval.res)), sep="/")

  if(numCores>1 & parallelTunning==FALSE){

    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl,list("mxnt.c","f.args"))

    mxnt.mdls.preds.lst <- parallel::clusterApply(cl, base::seq_along(ENMeval.o.l), function(i, ENMeval.o.l, a.calib.l, occ.l, use.ENMeval.bgpts, formt, pred.args, wAICsum, randomseed, responsecurves, arg1, arg2, numCores, parallelTunning){
      ## TODO - check this, decide if keep other fields before or remove only here (in which use loop to get)
      # ENMeval.o.l[[i]] <- ENMeval.o.l[[i]]@results

      # a <- NULL
      # if(use.ENMeval.bgpts){ a <- ENMeval.o.l[[i]]@bg.pts }

      cat(c(names(ENMeval.o.l[i]), "\n"))
      # if(dir.exists(path.mdls[i])==FALSE) dir.create(path.mdls[i])
      # compute final models and predictions
      resu <- mxnt.c(ENMeval.o = ENMeval.o.l[[i]], sp.nm = names(ENMeval.o.l[i]),
                     a.calib = a.calib.l[[i]], # a.proj = a.proj.l[[i]],
                     occ = occ.l[[i]], use.ENMeval.bgpts = use.ENMeval.bgpts, # a=ENMeval.o.l[[i]]@bg.pts,
                     formt = formt,
                     pred.args = pred.args, mSel = mSel, wAICsum = wAICsum,
                     randomseed = randomseed, responsecurves = responsecurves, arg1 = arg1, arg2 = arg2,
                     numCores = numCores, parallelTunning = parallelTunning)

      return(resu)
    }, ENMeval.o.l, a.calib.l, occ.l, use.ENMeval.bgpts, formt, pred.args, mSel, wAICsum, randomseed, responsecurves, arg1, arg2, numCores, parallelTunning)

    parallel::stopCluster(cl)

  } else {

    mxnt.mdls.preds.lst <- lapply(base::seq_along(ENMeval.o.l), function(i, ENMeval.o.l, a.calib.l, occ.l, use.ENMeval.bgpts, formt, pred.args, wAICsum, randomseed, responsecurves, arg1, arg2, numCores, parallelTunning){
      ## TODO - check this, decide if keep other fields before or remove only here (in which use loop to get)
      # ENMeval.o.l[[i]] <- ENMeval.o.l[[i]]@results
      cat(c(names(ENMeval.o.l[i]), "\n"))
      # if(dir.exists(path.mdls[i])==FALSE) dir.create(path.mdls[i])
      # compute final models and predictions
      resu <- mxnt.c(ENMeval.o = ENMeval.o.l[[i]], sp.nm = names(ENMeval.o.l[i]),
                     a.calib = a.calib.l[[i]], # a.proj = a.proj.l[[i]],
                     occ = occ.l[[i]], use.ENMeval.bgpts = use.ENMeval.bgpts, formt = formt,
                     pred.args = pred.args, mSel = mSel, wAICsum = wAICsum,
                     randomseed = randomseed, responsecurves = responsecurves, arg1 = arg1, arg2 = arg2,numCores=numCores,parallelTunning=parallelTunning)

      return(resu)
    }, ENMeval.o.l, a.calib.l, occ.l, use.ENMeval.bgpts, formt, pred.args, wAICsum, randomseed, responsecurves, arg1, arg2, numCores, parallelTunning)


  }

  names(mxnt.mdls.preds.lst) <- names(ENMeval.o.l)

  # mxnt.mdls.preds.lst <- vector("list", length(ENMeval.o.l))
  # names(mxnt.mdls.preds.lst) <- names(ENMeval.o.l)

  # for(i in base::seq_along(ENMeval.o.l)){
  #   ## TODO - check this, decide if keep other fields before or remove only here (in which use loop to get)
  #   ENMeval.o.l[[i]] <- ENMeval.o.l[[i]]@results
  #   cat(c(names(mxnt.mdls.preds.lst)[i], "\n"))
  #   # if(dir.exists(path.mdls[i])==FALSE) dir.create(path.mdls[i])
  #   # compute final models and predictions
  #   mxnt.mdls.preds.lst[[i]] <- mxnt.c(x = ENMeval.o.l[[i]], sp.nm = names(ENMeval.o.l[i]),
  #                                       a.calib = a.calib.l[[i]], # a.proj = a.proj.l[[i]],
  #                                       occ = occ.l[[i]], formt = formt,
  #                                       pred.args = pred.args, wAICsum = wAICsum,
  #                                       randomseed = randomseed, responsecurves = responsecurves, arg1 = arg1, arg2 = arg2,numCores=numCores,parallelTunning=parallelTunning)
  #   # mxnt.mdls.preds.lst[[i]]$pred.args <- pred.args
  # }

  return(mxnt.mdls.preds.lst)
}




