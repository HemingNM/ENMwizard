
#########	Omission Rate for AICc Averaged Model #############
#' Compute Omission Rate for a species' ensembled model
#'
#' This function will compute the omission rate (OR) for a species' ensembled model
#' from a 'mcmp' object, based on the selected threshold value.
#'
#' @param ORt Threshold value to be used to compute OR
#' @inheritParams calib_mdl
#' @inheritParams proj_mdl
#' @inheritParams ENMeval::ENMevaluate
#' @inheritParams ENMeval::modelTune.maxentJar
#' @seealso \code{\link{get_or_ensemble_b}}, \code{\link{get_tsa}}, \code{\link{get_cont_permimport}},
#' \code{\link{get_fpa}}, \code{\link{get_cont_permimport_b}}, \code{\link{get_fpa_b}}, \code{\link{get_tsa_b}}
#' @return Data frame with average and variance of OR values across partition groups of data
# #' @examples
#' @export
get_or_ensemble <- function(mcm, a.calib,
                           ORt=seq(0, 0.2, 0.05), userArgs=NULL, categoricals, sp.nm="species"){

  allMaxentArgs <- c("addsamplestobackground", "addallsamplestobackground", "allowpartialdata",
                     "beta_threshold", "beta_categorical", "beta_lqp", "beta_hinge", "convergencethreshold",
                     "defaultprevalence", "extrapolate", "fadebyclamping", "jackknife", "maximumbackground",
                     "maximumiterations", "removeduplicates")
  if (length(userArgs) == 0) {
    userArgs <- NULL
  } else {
    if (!all(names(userArgs) %in% allMaxentArgs)) {
      stop("The maxent argument given is not implemented in ENMeval or is misspelled.")
    } else {
      userArgs <- paste(names(userArgs), unlist(userArgs), sep='=')
    }
  }

  xsel.mdls <- mcm$selected.mdls

  ens <- c("AvgAIC", "WAAUC", "EBPM", "ESOR")
  ens.i <- grepl(paste0("^", mcm$mSel, collapse = "|^"), ens)
  ens <- ens[ens.i]
  statsTbl2 <- stats::setNames(vector("list", length(ens)), ens)
  for(EM in ens){
    #### 4.3.2.1.2 create model averaged prediction (models*weights, according to model selection)
    argsEns <- grep(EM, xsel.mdls$sel.cri)
    # argsEns <- grep(EM, mcm$mSel)

    # create vector of model weights for ensemble
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
  # wv <- xsel.mdls$w.AIC[grep("AIC", xsel.mdls$sel.cri)]
  # # pred.args <- mcm$pred.args
  # # ensemble model args
    args <- mcm$mxnt.args[argsEns]
  # args <- mcm$mxnt.args[grep("AIC", xsel.mdls$sel.cri)]

  group.data <- list(occ.grp = mcm$occ.grp,
                     bg.grp = mcm$bg.grp)
  pres <- as.data.frame(raster::extract(a.calib, mcm$occ.pts))
  bg <- as.data.frame(raster::extract(a.calib, mcm$bg.pts))

  nk <- length(unique(group.data$bg.grp))

  # set up empty data.frame for stats
  statsTbl <- as.data.frame(matrix(nrow = nk, ncol = length(ORt)))
  colnames(statsTbl) <- paste0("or", ORt)

  # cross-validation on partitions
  for (k in 1:nk) {
    # set up training and testing data groups
    train.val <- pres[group.data$occ.grp != k,, drop = FALSE]
    test.val <- pres[group.data$occ.grp == k,, drop = FALSE]
    bg.val <- bg[group.data$bg.grp != k,, drop = FALSE]
    occ.test <- mcm$occ.pts[group.data$occ.grp == k,, drop = FALSE]

    # redefine x and p for partition groups
    x <- rbind(train.val, bg.val)
    p <- c(rep(1, nrow(train.val)), rep(0, nrow(bg.val)))

    # run the current test model
    p.train <- data.frame(matrix(nrow = length(args), ncol = nrow(train.val)))
    p.test <- data.frame(matrix(nrow = length(args), ncol = nrow(test.val)))

    # predict values for training and testing data for each AICc model
    for(i in seq_along(args)){
      # create mod with train data
      mod <- dismo::maxent(x, p, args = c(args[[i]], userArgs), factors = categoricals)
      # predict values for training and testing data
      p.train[i,] <- dismo::predict(mod, train.val, args = mcm$pred.args)
      p.test[i,] <- dismo::predict(mod, test.val, args = mcm$pred.args)
    }

    p.train.wm <- apply(p.train, 2, stats::weighted.mean, w=wv)
    p.test.wm <- apply(p.test, 2, stats::weighted.mean, w=wv)

    # compute OR for each threshold value
    for(t in seq_along(ORt)) {
      # figure out % of total no. of training records
      if (nrow(train.val) < 10) {
        nt <- floor(nrow(train.val) * (1-ORt[t]))
      } else {
        nt <- ceiling(nrow(train.val) * (1-ORt[t]))
      }
      train.thr <- rev(sort(p.train.wm))[nt]
      statsTbl[k, t] <- mean(p.test.wm < train.thr)
    }
  }

  statsTbl <- data.table::dcast(data.table::melt(cbind(data.frame(
    stat=c("avgOR", "varOR"),
    apply(statsTbl, 2, function(x)c(mean(x), stats::var(x)) )
  )), id.vars = "stat"), variable~stat)

  statsTbl <- cbind(ORt, statsTbl[,-1])
  # statsTbl <- cbind(variable=statsTbl[,1], ORt, statsTbl[,-1])
  # colnames(statsTbl)[c(3,4)] <- c("avgOR", "varOR")
  statsTbl2[[EM]] <- statsTbl
  }
  statsTbl2 <- data.table::rbindlist(statsTbl2, idcol = "ensemble")
  utils::write.csv(statsTbl2, paste0("3_out.MaxEnt/Mdls.", sp.nm, "/metric.or.ensemble.mdl.", sp.nm, ".csv")) # reorder ds
  # utils::write.csv(ensOR.lst.c, paste0("3_out.MaxEnt/metric.or.ensemble.mdl.csv")) # reorder ds
  return(statsTbl2)
}




#' Compute Omission Rate for a list of species' ensembled models
#'
#' This function will compute the omission rate (OR) for each species' ensembled model
#' from a 'mcmp.l' object, based on the selected threshold value.
#'
#' @inheritParams get_or_ensemble
#' @inheritParams calib_mdl_b
#' @inheritParams proj_mdl_b
#' @inheritParams ENMeval::ENMevaluate
#' @inheritParams ENMeval::modelTune.maxentJar
#' @seealso \code{\link{get_or_ensemble}}
#' @return Data frame with average and variance of OR values across partition groups of data
# #' @examples
#' @export
get_or_ensemble_b <- function(mcm.l, a.calib.l,
                             ORt=seq(0, 0.2, 0.05), userArgs=NULL, categoricals=NULL,
                             numCores = 1){
  if(numCores>1){

    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl, list("get_or_ensemble"))

    ensOR.lst <- parallel::clusterApply(cl, base::seq_along(mcm.l),
                                        function(i, mcm.l, a.calib.l,
                                                 ORt, userArgs, categoricals){
                                          cat(c(names(mcm.l[i]), "\n"))
                                          # Compute Omission Rate for a species' AICc Averaged Model
                                          resu <- get_or_ensemble(mcm = mcm.l[[i]], a.calib = a.calib.l[[i]],
                                                           ORt=ORt, userArgs=userArgs, categoricals=categoricals,
                                                           sp.nm=names(mcm.l[i]))

                                          return(resu)
                                        }, mcm.l, a.calib.l,
                                        ORt, userArgs, categoricals)

    parallel::stopCluster(cl)

  } else {

    ensOR.lst <- lapply(base::seq_along(mcm.l), function(i, mcm.l, a.calib.l,
                                                         ORt, userArgs, categoricals){
      cat(c(names(mcm.l[i]), "\n"))
      # Compute Omission Rate for a species' AICc Averaged Model
      resu <- get_or_ensemble(mcm = mcm.l[[i]], a.calib = a.calib.l[[i]],
                       ORt=ORt, userArgs=userArgs, categoricals=categoricals,
                       sp.nm=names(mcm.l[i]))
      return(resu)
    }, mcm.l, a.calib.l,
    ORt, userArgs, categoricals)
   }
  names(ensOR.lst) <- names(mcm.l)
  ensOR.lst.c <- data.table::rbindlist(ensOR.lst, idcol = "sp")
  utils::write.csv(ensOR.lst.c, paste0("3_out.MaxEnt/metric.or.ensemble.mdl.csv")) # reorder ds
  return(ensOR.lst.c)
}


