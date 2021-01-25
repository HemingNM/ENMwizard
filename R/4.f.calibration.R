## 4.2 Generate final models for occ
#' Model selection and creation of MaxEnt arguments for selected models
#'
#' This function will read an object of class ENMevaluation (See ?ENMeval::ENMevaluate for details) and
#' return the results table with models selected by the chosen criteria. It can also return the
#' necessary arguments for final model calibration and predictions.
#'
#' @param x Object of class ENMevaluation
#' @param mSel Character vector. Which criteria to use when selecting model(s). Currently implemented:
#' "AvgAIC", "LowAIC", "OR", "AUC"
#' @param wAICsum Cumulative sum of top ranked models for which arguments will be created
#' @param dAICc Maximum delta AICc of models to be selected.
#' @param save Should save args only ("A"), selected models only ("M") or both ("B")?
# #' @inheritParams dismo::maxent
#' @param randomseed logical. Args to be passed to dismo::maxent. See ?dismo::maxent and the MaxEnt help for more information.
#' @param responsecurves logical. Args to be passed to dismo::maxent. See ?dismo::maxent and the MaxEnt help for more information.
#' @param arg1 charater. Args to be passed to dismo::maxent. See ?dismo::maxent and the MaxEnt help for more information.
#' @param arg2 charater. Args to be passed to dismo::maxent. See ?dismo::maxent and the MaxEnt help for more information.
#' @seealso \code{\link{calib_mdl}}, \code{\link{calib_mdl_b}}, \code{\link[dismo]{maxent}}, \code{\link[ENMeval]{ENMevaluate}},
#' \code{\link{proj_mdl}}, \code{\link{proj_mdl_b}}
#' @examples
#' \dontrun{
#' ENMeval.res.lst <- ENMevaluate(occ.locs, occ.b.env, parallel = F , numCores = 1)
#' mod_sel(ENMeval.res.lst[[1]])
#' }
#' @return A vector of args (if save="A"), data.frame of selected models (if save="M") or
#' a list with both, args and selected models, (if save="B")
# #' @keywords internal
#' @export
mod_sel <- function(x, mSel=c("AvgAIC", "EBPM", "WAAUC", "ESORIC", "LowAIC", "OR", "AUC"), wAICsum=0.99, dAICc=2, randomseed=FALSE, responsecurves=TRUE, arg1='noaddsamplestobackground', arg2='noautofeature', save="M"){ # , seq=TRUE
  x <- x@results
  x$sel.cri <- ""
  x$ID <- as.numeric(rownames(x))

  if ("All" %in% mSel) {
    x$sel.cri <- sub("^\\.", "",
                     paste(x$sel.cri, paste0("All_", x$ID), sep = "."))
  }

  x.a.i <- order(x$delta.AICc)
  x <- x[x.a.i,]
  if(is.null(x$rankAIC)) {x$rankAIC <- 1:nrow(x)} # if something goes wrong with function calib_mdl, check this line

  # AvgAIC
  if("AvgAIC" %in% mSel){
    if(any(cumsum(x$w.AIC) >= wAICsum)){
      wsum <- 1:(which(cumsum(x$w.AIC) >= wAICsum)[1])
    } else {
      wsum <- 1:length(x$w.AIC)
    }

    x$sel.cri[wsum] <- paste0(x$sel.cri[wsum], paste0("AIC_", wsum))

    cat("\n", paste(length(wsum)), "of", nrow(x), "models selected using AIC")# from a total of", "models")
    cat("\n", "Total AIC weight (sum of Ws) of selected models is", round(sum(x$w.AIC[wsum]), 4), "of 1")

  }
  x <- x[order(x$ID),]

  # WA consensus through AUC - Marmion et al 2009
  # get half of top AUC models
  if("WAAUC" %in% mSel){
    xAUCmtp.WA <- order(-x$avg.test.AUC, x$avg.test.orMTP)
    xAUCmtp.WA <- xAUCmtp.WA[1:(round(length(xAUCmtp.WA)/2))]
    x$sel.cri[xAUCmtp.WA] <- sub("^\\.", "", paste(x$sel.cri[xAUCmtp.WA], paste0("WAAUC_", 1:length(xAUCmtp.WA)), sep = "."))
    # mji2 <- mji[1:(round(length(xAUCmtp.WA)/2))]
    # sum(auc2*mji2)/sum(auc2)
    # x$sel.cri[xAUCmtp.i] <- sub("^\\.", "", paste(x$sel.cri[xAUCmtp.i], "AUCmtp", sep = "."))
  }

  # EBPM (Ensemble of best-performing models) - Boria et al 2016
  # 10% top-performing models based on the sequential criteria:
  # Lowest average omission rate (OR) and, subsequently, the highest average AUCevaluation
  if("EBPM" %in% mSel){
    xORm <- order(x$avg.test.orMTP, -x$avg.test.AUC) # order(-x$avg.test.AUC, x$avg.test.orMTP)
    EBPM <- xORm[ifelse(round(length(xORm)*.1)==0, 1, 1:round(length(xORm)*.1))]
    # if(length(EBPM)==1){
    #   print("only one model selected, Not performing EBPM")
    # } else {
    x$sel.cri[EBPM] <- sub("^\\.", "", paste(x$sel.cri[EBPM], paste0("EBPM_", 1:length(EBPM)), sep = "."))
    # }
  }

  # ESOR (Ensemble Significant pROC, low Omission Rate) - Cobos et al 2019
  # Selects models according to: pROC<=0.05; OR<=ORspecified; AICc<=2
  if(any(c("ESOR", "ESORIC") %in% mSel)){
    mod.cobos <- x #[!is.na(x$avg.pROC.p),]
    mod.cobos[is.na(mod.cobos$avg.pROC.p),"avg.pROC.p"] <- 1
    ESOR.pROC <- mod.cobos$avg.pROC.p<=0.05
    if(sum(ESOR.pROC, na.rm = T)==0) {
      ESOR.pROC <- mod.cobos$avg.pROC.p == min(mod.cobos$avg.pROC.p, na.rm = T)
      warning("ESORIC: No model with pROC p value <= 0.05. Using model with lowest pROC.p")
    }
    ESOR.or <- mod.cobos$avg.test.or10pct<=0.1
    if(sum(ESOR.or, na.rm = T)==0) {
      ESOR.or <- mod.cobos$avg.test.or10pct == min(mod.cobos$avg.test.or10pct, na.rm = T)
      warning("ESORIC: No model with OR <= OR criteria. Using model with lowest OR")
    }
    ESOR <- ESOR.pROC & ESOR.or
    if(sum(ESOR, na.rm = T)==0) {
      ESOR <- mod.cobos$avg.test.or10pct == min(mod.cobos$avg.test.or10pct[ESOR.pROC], na.rm = T) #mod.cobos$avg.pROC.p<=0.05 &
      warning("ESORIC: No model match pROC criteria. Using model with lowest OR from models with lowest pROC.p")
    }
    delta <- ifelse((mod.cobos$AICc - min(mod.cobos$AICc[ESOR], na.rm = T)<0), NA, (mod.cobos$AICc - min(mod.cobos$AICc[ESOR], na.rm = T)))
    ESOR <- which(delta<=dAICc & delta>=0 & ESOR)
    ESOR <- ESOR[order(delta[ESOR])]
    x$sel.cri[ESOR] <- sub("^\\.", "", paste(x$sel.cri[ESOR], paste0("ESORIC_", 1:length(ESOR)), sep = "."))
  }

  # # Mean ALL - Marmion et al 2009
  # mean(mji)
  # # Median ALL - Marmion et al 2009
  # median(mji)
  #
  # # Median(PCA) - Marmion et al 2009
  # # Runs PCA on suitability values:
  # # Selects half models for which the variance of predictions along PC1 is the greatest


  # LowAIC
  if("LowAIC" %in% mSel){
    x.la.i <- x.a.i[1]
    x$sel.cri[x.la.i] <- sub("^\\.", "", paste(x$sel.cri[x.la.i], "LowAIC", sep = "."))
  } else {
    x.la.i <- NULL
  }

  # OR
  if("OR" %in% mSel){
    # seq
    xORm.i <- order(x$avg.test.orMTP, -x$avg.test.AUC)[1]
    x$sel.cri[xORm.i] <- sub("^\\.", "", paste(x$sel.cri[xORm.i], "ORmtp", sep = "."))
    # seqOr10
    xOR10.i <- order(x$avg.test.or10pct, -x$avg.test.AUC)[1]
    x$sel.cri[xOR10.i] <- sub("^\\.", "", paste(x$sel.cri[xOR10.i], "OR10", sep = "."))
  } else {
    xORm.i <- NULL
    xOR10.i <- NULL
  }

  # AUC
  if("AUC" %in% mSel){
    #AUCmtp
    xAUCmtp.i <- order(-x$avg.test.AUC, x$avg.test.orMTP)[1]
    x$sel.cri[xAUCmtp.i] <- sub("^\\.", "", paste(x$sel.cri[xAUCmtp.i], "AUCmtp", sep = "."))
    # AUC10
    xAUC10.i <- order(-x$avg.test.AUC, x$avg.test.or10pct)[1]
    x$sel.cri[xAUC10.i] <- sub("^\\.", "", paste(x$sel.cri[xAUC10.i], "AUC10", sep = "."))
  } else {
    xAUCmtp.i <- NULL
    xAUC10.i <- NULL
  }

  xsel.mdls <- x#[x$sel.cri!="",]

  f <- factor(xsel.mdls$features)
  beta <- c(xsel.mdls$rm)
  cat("\n", "arguments used for building models", "\n")
  print(data.frame(ID=xsel.mdls$ID, optimality.criteria = xsel.mdls$sel.cri, settings=xsel.mdls$settings))
  xsel.mdls$ID <- NULL

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

  mod.nms <-  paste0("Mod_", format(xsel.mdls[, "rm"], digits=2), "_", xsel.mdls[, "features"]) #
  args <- stats::setNames(c(strsplit(args, ",")), mod.nms)
  # args <- stats::setNames(c(strsplit(args, ",")), xsel.mdls$settings)

  if(save == "A"){
    return(args)
  } else if(save == "M"){
    return(xsel.mdls)
  } else if (save == "B"){
    return(list(args=args, mdls=xsel.mdls))
  }
}

#### 4.3 Run top corresponding models and save predictions
#### 4.3.1 save maxent best models and predictions for each model
# "f.mxnt.mdl.pred" renamed to "calib_mdl" and now to "calib_mdl"
#' Calibrate MaxEnt models based on model selection criteria
#'
#' This function will read an object of class ENMevaluation (See ?ENMeval::ENMevaluate for details) and
#' calibrate the selected maxent models.
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
#' @param format Character. Output file type. Argument 'format' of raster::writeRaster
#' @param pred.args Charater. Argument 'args' of dismo::maxent. Additional argument that can be
#' passed to MaxEnt. See the MaxEnt help for more information. The R maxent function only uses the
#' arguments relevant to model fitting. There is no point in using args='outputformat=raw' when
#' *fitting* the model; but you can use arguments relevant for *prediction* when using the predict
#' function. Some other arguments do not apply at all to the R implementation. An example is
#' 'outputfiletype', because the 'predict' function has its own 'filename' argument for that.
#' @param numCores Number of cores to use for parallelization. If set to 1, no paralellization is performed
#' @param parallelTunning Should parallelize within species (parallelTunning=TRUE) or between species (parallelTunning=FALSE)
#' @param use.ENMeval.bgpts Logical. Use background points from ENMeval or sample new ones?
#' @inheritParams mod_sel
#' @inheritParams dismo::maxent
#' @seealso \code{\link{mod_sel}}, \code{\link{calib_mdl_b}}, \code{\link{ENMevaluate_b}}, \code{\link[ENMeval]{ENMevaluate}},
#' \code{\link[dismo]{maxent}}, \code{\link{proj_mdl}}, \code{\link{proj_mdl_b}}
#' @return A 'mcm' (calib.mdls, Maxent Calibrated Models). A list containing the models ('selected.mdls') used for model calibration,
#' calibrated maxent models ('mxnt.mdls'), and arguments used for calibration ('pred.args').
#' @export
calib_mdl <- function(ENMeval.o, sp.nm = "species", a.calib, occ = NULL, use.ENMeval.bgpts = TRUE, nbg=10000, format = "raster", # , a.proj
                   pred.args = c("outputformat=cloglog", "doclamp=true", "pictures=true"),
                   mSel = c("AvgAIC", "LowAIC", "OR", "AUC"), wAICsum = 0.99, dAICc=2, randomseed = FALSE,
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

  # ENMeval.r <- ENMeval.o@results
  algorithm <- ENMeval.o@algorithm

  mdl.arg <- mod_sel(x=ENMeval.o, mSel=mSel, wAICsum=wAICsum, dAICc=dAICc, randomseed=randomseed, responsecurves=responsecurves, arg1=arg1, arg2=arg2, save="B")
  xsel.mdls <- mdl.arg[[2]]
  ENMeval.r <- xsel.mdls[order(as.numeric(rownames(xsel.mdls))),]
  mdls.keep <- xsel.mdls$sel.cri!=""
  xsel.mdls <- xsel.mdls[mdls.keep,]
  args.all <- mdl.arg[[1]][mdls.keep]

  # exportar planilha de resultados
  utils::write.csv(ENMeval.r, paste0(path.mdls,"/sel.mdls.", gsub("3_out.MaxEnt/Mdls.", "", path.mdls), ".csv"))
  if("avg.pROC.p" %in% colnames(xsel.mdls)){
    res.tbl <- xsel.mdls[,c("sel.cri", "features","rm","AICc", "delta.AICc", "w.AIC", "parameters", "rankAIC", "avg.test.pROC.ratio", "avg.pROC.p",
                            "avg.test.or10pct", "avg.test.orMTP", "avg.test.AUC")]
    colnames(res.tbl) <- c("Optimality criteria", "FC", "RM", "AICc", "delta AICc", "wAIC", "parameters", "Rank AIC", "avg pROC ratio", "avg pROC pValue",
                           "avg OR10", "avg ORLPT", "avg AUC")
  } else {
    res.tbl <- xsel.mdls[,c("sel.cri", "features","rm","AICc", "delta.AICc", "w.AIC", "parameters", "rankAIC", # "avg.test.pROC.ratio", "avg.pROC.p",
                            "avg.test.or10pct", "avg.test.orMTP", "avg.test.AUC")]
    colnames(res.tbl) <- c("Optimality criteria", "FC", "RM", "AICc", "delta AICc", "wAIC", "parameters", "Rank AIC", # "avg pROC ratio", "avg pROC pValue",
                           "avg OR10", "avg ORLPT", "avg AUC")

  }
  utils::write.csv(res.tbl, paste0(path.mdls,"/sel.mdls.smmr.", gsub("3_out.MaxEnt/Mdls.", "", path.mdls), ".csv"))

  mod.nms <- paste0("Mod_", format(xsel.mdls[, "rm"], digits=2), "_", xsel.mdls[, "features"]) #
  # mod.nms <- paste0("Mod.", xsel.mdls[, "settings"])
  # mod.nms <- paste0("Mod.", xsel.mdls[, "sel.cri"])
  mod.preds <- raster::stack()

  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

  if(dir.exists(paste(path.mdls, outpt, sep='/'))==FALSE) dir.create(paste(path.mdls, outpt, sep='/'))

  # mxnt.mdls <- setNames(vector("list", length(args.all)), names(args.all))

  #### Calibrate ALL selected models
  {
    if(numCores>1 & parallelTunning){

      cl <- parallel::makeCluster(numCores)

      mxnt.mdls <- parallel::clusterApply(cl, seq_along(args.all), function(i, args.all, pred.args, mod.nms, a.calib, occ, a) {

        path2file <- paste(path.mdls, outpt, mod.nms[i], sep='/')

        filename <- paste(path2file, paste0(mod.nms[i], ".grd"), sep='/')
        # maxent models
        set.seed(1)
        resu <- dismo::maxent(a.calib, occ, a, path=path2file, args=c(args.all[[i]], pred.args)) # final model fitting/calibration
        return(resu)
      }, args.all, pred.args, mod.nms, a.calib, occ, a) #)

      parallel::stopCluster(cl)



    } else {
      mxnt.mdls <- lapply(seq_along(args.all), function(i, args.all, pred.args, mod.nms, a.calib, occ, a) {

        path2file <- paste(path.mdls, outpt, mod.nms[i], sep='/')

        filename <- paste(path2file, paste0(mod.nms[i], ".grd"), sep='/')
        # maxent models
        set.seed(1)
        resu <- dismo::maxent(a.calib, occ, a, path=path2file, args=c(args.all[[i]], pred.args)) # final model fitting/calibration
        return(resu)
      }, args.all, pred.args, mod.nms, a.calib, occ, a) #) ## fecha for or lapply

    } # closes else
    mxnt.mdls <- stats::setNames(mxnt.mdls, names(args.all))

  }
  return(list(algorithm = algorithm,
              ENMeval.results = ENMeval.r, mxnt.mdls = mxnt.mdls,
              selected.mdls = xsel.mdls, mSel = mSel,
              occ.pts = occ, bg.pts = a,
              occ.grp = ENMeval.o@occ.grp, bg.grp = ENMeval.o@bg.grp,
              mxnt.args = args.all, pred.args = pred.args)) #, mxnt.preds = mod.preds))
}

# "f.mxnt.mdl.pred.batch" renamed to "calib_mdl_b" and now to "calib_mdl_b"


#' Calibrate MaxEnt models based on model selection criteria for several species
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
# #' @param a.proj.l List of projection areas. See argument "a.proj" in calib_mdl.
#' @param occ.l List of occurence data. See argument "occ" in calib_mdl.
# #' @param numCores Number of cores to use for parallelization. If set to 1, no paralellization is performed
#' @inheritParams calib_mdl
#' @inheritParams set_calibarea_b
#' @seealso \code{\link{mod_sel}}, \code{\link{calib_mdl}}, \code{\link{ENMevaluate_b}}, \code{\link[ENMeval]{ENMevaluate}},
#' \code{\link[dismo]{maxent}}, \code{\link{proj_mdl}}, \code{\link{proj_mdl_b}}
#' @return A 'mcm.l' object. A list of 'mcm' (calib.mdls, Maxent Calibrated Models), returned from function "calib_mdl"
#' @examples
#' \dontrun{
#' mxnt.mdls.preds.lst <- calib_mdl_b(ENMeval.o.l=ENMeval.res.lst,
#' a.calib.l=occ.b.env, a.proj.l=areas.projection, occ.l=occ, wAICsum=0.99)
#' mxnt.mdls.preds.lst[[1]][[1]] # models [ENMevaluate]d and selected using sum of wAIC
#' mxnt.mdls.preds.lst[[1]][[2]] # MaxEnt models
#' mxnt.mdls.preds.lst[[1]][[3]] # used prediction arguments
#' plot(mxnt.mdls.preds.lst[[1]][[4]]) # MaxEnt predictions, based on the model selection criteria
#' }
#' @export
calib_mdl_b <- function(ENMeval.o.l, a.calib.l, occ.l = NULL, use.ENMeval.bgpts = TRUE, format = "raster", # , a.proj.l
                         pred.args = c("outputformat=cloglog", "doclamp=true", "pictures=true"),
                         mSel = c("AvgAIC", "LowAIC", "OR", "AUC"), wAICsum = 0.99, dAICc=2, randomseed = FALSE,
                         responsecurves = TRUE, arg1 = 'noaddsamplestobackground', arg2 = 'noautofeature',
                         numCores = 1, parallelTunning = TRUE){

  if(numCores>1 & parallelTunning==FALSE){

    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl, list("calib_mdl","mod_sel"))

    mxnt.m.p.lst <- parallel::clusterApply(cl, base::seq_along(ENMeval.o.l), function(i, ENMeval.o.l, a.calib.l, occ.l,
                                                                                             use.ENMeval.bgpts, format, pred.args,
                                                                                             mSel, wAICsum, dAICc, randomseed, responsecurves,
                                                                                             arg1, arg2, numCores, parallelTunning){
      cat(c(names(ENMeval.o.l[i]), "\n"))
      # compute final models and predictions
      resu <- calib_mdl(ENMeval.o = ENMeval.o.l[[i]], sp.nm = names(ENMeval.o.l[i]),
                     a.calib = a.calib.l[[i]], # a.proj = a.proj.l[[i]],
                     occ = occ.l[[i]], use.ENMeval.bgpts = use.ENMeval.bgpts, # a=ENMeval.o.l[[i]]@bg.pts,
                     format = format,
                     pred.args = pred.args, mSel = mSel, wAICsum = wAICsum, dAICc,
                     randomseed = randomseed, responsecurves = responsecurves, arg1 = arg1, arg2 = arg2,
                     numCores = numCores, parallelTunning = parallelTunning)

      return(resu)
    }, ENMeval.o.l, a.calib.l, occ.l,
    use.ENMeval.bgpts, format, pred.args,
    mSel, wAICsum, dAICc, randomseed, responsecurves,
    arg1, arg2, numCores, parallelTunning)

    parallel::stopCluster(cl)

  } else {

    mxnt.m.p.lst <- lapply(base::seq_along(ENMeval.o.l), function(i, ENMeval.o.l, a.calib.l, occ.l,
                                                                         use.ENMeval.bgpts, format, pred.args,
                                                                         mSel, wAICsum, dAICc, randomseed, responsecurves,
                                                                         arg1, arg2, numCores, parallelTunning){
      cat(c(names(ENMeval.o.l[i]), "\n"))
      # compute final models and predictions
      resu <- calib_mdl(ENMeval.o = ENMeval.o.l[[i]], sp.nm = names(ENMeval.o.l[i]),
                     a.calib = a.calib.l[[i]], # a.proj = a.proj.l[[i]],
                     occ = occ.l[[i]], use.ENMeval.bgpts = use.ENMeval.bgpts, format = format,
                     pred.args = pred.args, mSel = mSel, wAICsum = wAICsum, dAICc,
                     randomseed = randomseed, responsecurves = responsecurves,
                     arg1 = arg1, arg2 = arg2, numCores=numCores, parallelTunning=parallelTunning)

      return(resu)
    }, ENMeval.o.l, a.calib.l, occ.l,
    use.ENMeval.bgpts, format, pred.args,
    mSel, wAICsum, dAICc, randomseed, responsecurves,
    arg1, arg2, numCores, parallelTunning)


  }
  names(mxnt.m.p.lst) <- names(ENMeval.o.l)
  return(mxnt.m.p.lst)
}




