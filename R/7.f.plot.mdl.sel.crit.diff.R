#### 4.4 plot differences among model selection criteria predictions
#' Format (underscript) selected texts for plotting
#'
#' Format (underscript) selected texts (criteria used to select models) to be used on plotting.
#'
#' @param x list of text to be formatted
#' @seealso \code{\link{plotMdlDiff}}, \code{\link{plotScnDiff}}
#' @return list of formatted text
#' @examples
#' make.underscript(c("AUC (OR10p)", "AUC (ORlpt)", "OR10p (AUC)", "ORlpt (AUC)"))
#' @keywords internal
# #' @export
make.underscript <- function(x) as.expression(lapply(x, function(y) {
  sp <- grepl("ORlpt", y)
  cf <- grepl("OR10p", y)
  y <- unlist(strsplit(y, "OR10p|ORlpt"))
  y1 <- y[1]
  y2 <- y[2]
  if(cf){
    y <- bquote(.(y1)*OR["10P"]*.(y2))
  } else if(sp){
    y <- bquote(.(y1)*OR["LPT"]*.(y2))
  } else {
    y <- bquote(.(y))
  }
  y
}))


#### 4.4 plot differences between predictions of model selection criteria
#' Plot differences between predictions of models selected from several criteria
#'
#' Plot differences among predictions of species model projetions selected using distinct criteria.
#' This function will save the figures on pdf files in the folder "Mdls.thrshld/figs".
#'
#' @inheritParams f.thr.batch
# #' @param pred.nm name of prediction to be appended to the final name. Usually "pres", "past" or "fut".
#' @param mtp.l List of stack or brick of thresholded predictions
#' @param basemap Shapefile to be plotted with. Usually a continent or country shapefile
#' @seealso \code{\link{plotMdlDiff}}, \code{\link{plotScnDiff}}
#' @return won't return any object. Will save pdf's with differences among model predictions
#' @examples
#' f.plot.mxnt.preds(mxnt.mdls.preds.lst, mods.thrshld.lst, basemap=NewWorld)
#' @keywords internal
# #' @export
# f.plot.mxnt.preds <- function(mcmp.l, mtp.l, basemap=NULL){ #, pred.nm=""
#   { path.res <- "3_out.MaxEnt"
#   if(dir.exists(path.res)==FALSE) dir.create(path.res)
#   path.sp.m <- paste0("Mdls.", names(mcmp.l))
#   path.mdls <- paste(path.res, path.sp.m, sep="/")
#   pred.args <- mcmp.l[[1]]$pred.args}
#
#   # unlist(strsplit(names(mods.thrshld$binary[[2]]), "."))
#
#   thrshld.nms <- c("fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd")
#   thrshld.nms.mod <- paste(c("Mod.", paste0(".",thrshld.nms)), collapse ="|")
#   thrshld.nms <- paste(thrshld.nms, collapse = "|")
#
#   # ncol(comb.plots)
#
#   outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
#                   ifelse(grep("logistic", pred.args)==1, 'logistic',
#                          ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))
#   # pred.nm <- ifelse(pred.nm != "", paste0(".", pred.nm), pred.nm)
#   for(i in 1:length(mtp.l)){ # species
#     comb.plots <- utils::combn(raster::nlayers(mtp.l[[i]][[1]]$continuous[[2]]), 2)
#
#     thrshld.path <- paste(path.mdls[i], outpt, "Mdls.thrshld", "figs", sep='/')
#     if(dir.exists(thrshld.path)==FALSE) dir.create(thrshld.path)
#
#     mods.thrshld <- mtp.l[[i]]
#
#     for(l in 1:length(mods.thrshld$binary)){ # threshold criteria
#       thr.crt <- grep(thrshld.nms,  unique(unlist(strsplit(names(mods.thrshld$binary[[l]]), "."))), value=TRUE)
#       grDevices::pdf(paste(thrshld.path, paste0("Mod.diff.bin", ".", thr.crt, ".pdf"), sep='/'), # , pred.nm
#                      width = 20, height = 10)
#       graphics::par(mfrow=c(3,5), mar=c(2,1,2,1), oma = c(1, 1, 4, 1))
#
#       for(j in 1:ncol(comb.plots)){ #ncol(comb.plots)
#         # r.dif <- mods.thrshld$binary[[l]][[comb.plots[1,j]]] - mods.thrshld$binary[[l]][[comb.plots[2,j]]]#, col=c("red", "white", "blue")
#         r.dif <- raster::overlay(mods.thrshld$binary[[l]][[comb.plots[1,j]]], mods.thrshld$binary[[l]][[comb.plots[2,j]]], fun=function(r1,r2) {r1-r2})
#         n1 <- sub("_", ".", sub("min", "Min", sub("Mean", "M", sub("OR", "or", gsub(thrshld.nms.mod, "", names(mods.thrshld$binary[[l]][[comb.plots[1,j]]]))))))
#         n2 <- sub("_", ".", sub("min", "Min", sub("Mean", "M", sub("OR", "or", gsub(thrshld.nms.mod, "", names(mods.thrshld$binary[[l]][[comb.plots[2,j]]]))))))
#         main.nms <- paste0(n1, " vs. ", n2)
#         # r.dif, breaks= c(-1, -.33, .33, 1), col=c("blue", "white", "red")
#         raster::plot(r.dif, breaks= c(-1, -.33, .33, 1), col=c("blue", "white", "red"),# col=c("blue", "blue","white", "red", "red"), # smallplot=c(0.01, .81, 0.01, .99),
#                      main= main.nms, legend=FALSE)
#         if(j==1){
#           graphics::title(paste("Threshold criteria:", thr.crt), line = 2, outer = T, cex.main=2)
#         }
#         raster::plot(r.dif,  legend.only=TRUE, smallplot=c(.78, .79, .2, .8),  #horiz=TRUE,#  c(.79, .80, .2, .8)
#                      breaks= c(-1, -.34, .34, 1), col=c("blue", "white", "red"),
#                      # bg = "white",
#                      #col=c("blue", "blue","white", "red", "red"),
#                      axis.args=list(at=seq(-1, 1),
#                                     labels=c(n2, "equal", n1)))
#         if(!is.null(basemap)) sp::plot(basemap, add= T)
#       }
#
#       grDevices::dev.off()
#     }
#
#   }
#
# }


#### 4.8.6 plot prediction diff between models
#' Plot differences (for multiple climatic scenarios) among model predictions selected from several criteria
#'
#' Plot differences (for multiple climatic scenarios) among model predictions selected from several
#' criteria (e.g. "AvgAIC", "LowAIC", "OR", "AUC")
#' selected using distinct criteria. This function will save the figures on pdf files in the folder "Mdls.thrshld/figs".
#'
# #' @inheritParams f.plot.mxnt.preds
#' @inheritParams f.thr.batch
#' @param mtp.l List of stack or brick of thresholded predictions
#' @param basemap Shapefile to be plotted with. Usually a continent or country shapefile
#' @param save Logical. If TRUE will save plots in pdf.
#' @seealso \code{\link{plotScnDiff}}
#' @return won't return any object. Will save pdf's with differences among model predictions (for multiple climatic scenarios)
#' @examples
#' plotMdlDiff(mcmp.l=mxnt.mdls.preds.lst, mtp.l=mods.thrshld.lst, basemap=NewWorld)
#' plotMdlDiff(mcmp.l=mxnt.mdls.preds.pf[1], mtp.l=mods.thrshld.lst[1], basemap=NewWorld)
#' @export
plotMdlDiff <- function(mcmp.l, mtp.l, basemap=NULL, save=FALSE, numCores=1){
  { path.res <- "3_out.MaxEnt"
  if(dir.exists(path.res)==FALSE) dir.create(path.res)
  path.sp.m <- paste0("Mdls.", names(mcmp.l))
  path.mdls <- paste(path.res, path.sp.m, sep="/")
  pred.args <- mcmp.l[[1]]$pred.args}

  a <- c("avg.test.AUC10pct", "avg.test.AUC.MTP", "avg.test.or10pct", "avg.test.orMTP")
  b <- c("AUC (OR10p)", "AUC (ORlpt)", "OR10p (AUC)", "ORlpt (AUC)")
  sca <- c("cc26bi70", "cc45bi70", "cc60bi70", "cc85bi70", "mp26bi70", "mp45bi70", "mp85bi70", "mr26bi70", "mr45bi70", "mr60bi50", "mr85bi70", "cclgmbi", "ccmidbi", "lig_30s_bio_", "melgmbi", "memidbi", "mrlgmbi", "mrmidbi", "mxnt.preds")
  scb <- c("2070-CCSM4-rcp2.6", "2070-CCSM4-rcp4.5", "2070-CCSM4-rcp6.0", "2070-CCSM4-rcp8.5",
           "2070-MPI-ESM-LR-rcp2.6", "2070-MPI-ESM-LR-rcp4.5", "2070-MPI-ESM-LR-rcp8.5",
           "2070-MIROC-ESM-rcp2.6", "2070-MIROC-ESM-rcp4.5", "2070-MIROC-ESM-rcp6.0", "2070-MIROC-ESM-rcp8.5",
           "LGM-CCSM4", "MH-CCSM4", "LIG-CCSM3", "LGM-MPI-ESM-P", "MH-MPI-ESM-P", "LGM-MIROC-ESM", "MH-MIROC-ESM", "Present")

  t.nms <- c("fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd")
  t.NMS <- c("FCV1", "FCV5", "FCV10", "LPT (mtp)", "10P (x10ptp)", "ETSS (etss)", "MTSS", "BTO", "EETD")
  thrshld.nms.mod <- paste(c("Mod\\.", paste0("\\.",t.nms)), collapse ="|")
  # thrshld.nms.mod <- paste(c("Mod.", paste0(".",t.nms)), collapse ="|")
  thrshld.nms <- paste(t.nms, collapse = "|")

  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

  f.plot <- function(sc, sp, mtp.l, t.NMS, t.nms, thrshld.path,
                     comb.plots, thrshld.nms.mod, basemap, make.underscript) {
    # for(sc in names(mtp.l[[sp]])){ # climatic scenario
    mods.thrshld <- mtp.l[[sp]][[sc]]
    cat(c("\n", "Climatic Scenario: ", sc))
    cat(c("\n", "Threshold: "))

    for(l in names(mods.thrshld$binary)){ # threshold criteria
      thr.crt <- l
      thr.CRT <- t.NMS[which(t.nms %in% l)] #}
      cat(paste0(" - ", thr.CRT))

      n.t <- length(names(mods.thrshld$binary))
      n.scn <- ncol(comb.plots)
      if(save){
        grDevices::pdf(paste(thrshld.path, paste0("Mod.diff.bin.", sc, ".", thr.crt, ".pdf"), sep='/'),
                       width = n.scn*5+2, height = n.t*5)
      }
      graphics::par(mfrow=c(n.t, n.scn), oma = c(3.5, 0, 5.5, 2))

      for(j in 1:ncol(comb.plots)){ #ncol(comb.plots)
        r.dif <- raster::overlay(mods.thrshld$binary[[l]][[comb.plots[1,j]]], mods.thrshld$binary[[l]][[comb.plots[2,j]]], fun=function(r1,r2) {r1-r2})
        n1 <- gsub(paste0("\\.",sc), "", gsub(thrshld.nms.mod, "", names(mods.thrshld$binary[[l]][[comb.plots[1,j]]]) )  )
        n2 <- gsub(paste0("\\.",sc), "", gsub(thrshld.nms.mod, "", names(mods.thrshld$binary[[l]][[comb.plots[2,j]]]) )  )
        if(n1 %in% a) { n1 <- b[which(a %in% n1)] }
        if(n2 %in% a) { n2 <- b[which(a %in% n2)] }

        main.nms <- paste0( make.underscript(n1), " vs. ", make.underscript(n2))
        graphics::par(mar=c(2,4,2,4.7)) # par(mar=c(2,4,2,5))
        raster::plot(mods.thrshld$binary[[l]][[comb.plots[1,j]]], breaks= c(0, .5, 1), col=c("white", "gray90"),
                     legend=FALSE) # main= main.nms,
        if(j==1){
          sc1 <- gsub("mxnt.pred.fut.|mxnt.pred.past.", "", sc)
          if(sc1 %in% sca) { sc1 <-  scb[which(sca %in% sc1)] }
          graphics::title(paste0("Climatic Scenario: ", sub("mxnt.pred.", "", sub("mxnt.preds", "bioclim", sc1))), line = 1.5, outer = T, cex.main=2)
          graphics::title(paste0("Threshold: ", thr.CRT), line = -.5, outer = T, cex.main=2)
        }
        if(!is.null(basemap)) raster::plot(basemap, border="gray50", add= T)
        raster::plot(r.dif, breaks= c(-1, -.33, .33, 1), col=c("blue", grDevices::rgb(0,0,0,0), "red"),
                     legend=FALSE, add=TRUE) # main= main.nms,

        graphics::par(mar=c(2,1,2,6))
        raster::plot(r.dif,  legend.only=TRUE, legend.width=1.75, legend.shrink=.75,
                     xpd = TRUE, zlim=c(0, 1),
                     breaks= c(-1, -.34, .34, 1), col=c("blue", "gray90", "red"),
                     axis.args=list(at=seq(-1, 1), labels=c(make.underscript(n2), "equal", make.underscript(n1) )))
      }
      if(save){
        grDevices::dev.off()
      }
    }
  }

  for(sp in 1:length(mtp.l)){ # species
    N.mdls <- raster::nlayers(mtp.l[[sp]][[1]]$binary[[1]])

    if(N.mdls<2){
      m.nms <- names(mtp.l[[sp]][[1]]$binary[[1]])
      stop("No models to compare. \nJust one model selected: ", m.nms,
           "\nModel selection criteria: ", paste(mcmp.l[[sp]]$mSel, collapse = " "), sep="")
    }

    comb.plots <- utils::combn(N.mdls, 2)
    thrshld.path <- paste(path.mdls[sp], outpt, "Mdls.thrshld", "figs", sep='/')
    if(dir.exists(thrshld.path)==FALSE) dir.create(thrshld.path)
    cat(c("\n", "Species: " , names(mtp.l)[sp]))

    ## TODO use mclapply

    if(numCores>1){

      cl<-parallel::makeCluster(numCores)

      parallel::clusterApply(cl, names(mtp.l[[sp]]), # climatic scenario
                             function(sc, sp, mtp.l, t.NMS, t.nms, thrshld.path,
                                      comb.plots, thrshld.nms.mod, basemap, make.underscript){

                               f.plot(sc, sp, mtp.l, t.NMS, t.nms, thrshld.path,
                                      comb.plots, thrshld.nms.mod, basemap, make.underscript)

                             }, sp, mtp.l, t.NMS, t.nms, thrshld.path,
                             comb.plots, thrshld.nms.mod, basemap, make.underscript)

      parallel::stopCluster(cl)

    }else{
      lapply(names(mtp.l[[sp]]), # climatic scenario
             function(sc, sp, mtp.l, t.NMS, t.nms, thrshld.path,
                      comb.plots, thrshld.nms.mod, basemap, make.underscript){

               f.plot(sc, sp, mtp.l, t.NMS, t.nms, thrshld.path,
                      comb.plots, thrshld.nms.mod, basemap, make.underscript)

             }, sp, mtp.l, t.NMS, t.nms, thrshld.path,
             comb.plots, thrshld.nms.mod, basemap, make.underscript)

    }
  }
}





#' Plot differences between climatic scenarios
#'
#' Plot differences between a selected climatic scenario and all other climatic scenarios for each species.
#' This function will plota and (optionally) save the figures on pdf files in the folder "Mdls.thrshld/figs".
#'
#' @inheritParams plotMdlDiff
#' @inheritParams mxnt.c
#' @param sel.clim.scn Selected climatic scenario to compare with all others. Usually "current" one.
#' @param mSel Name of selection criteria to be compared: AvgAICc, LowAICc, avg.test.AUC10pct, avg.test.AUC.MTP,
#' avg.test.or10pct, avg.test.orMTP
#' @param save Export to pdf or not?
#' @seealso \code{\link{plotMdlDiff}}
#' @return won't return any object. Will save pdf's with differences among model predictions (for multiple climatic scenarios)
#' @examples
#' plotScnDiff(mcmp.l=mxnt.mdls.preds.cf, mtp.l=mods.thrshld.lst)
#' @export
plotScnDiff <- function(mcmp.l, mtp.l, mSel = mcmp.l[[1]]$mSel, sel.clim.scn="current",
                          basemap=NULL, save=FALSE, numCores=1){
  if(is.null(mSel)){
    stop("Need to specify 'mSel'")
  }

  { path.res <- "3_out.MaxEnt"
    if(dir.exists(path.res)==FALSE) dir.create(path.res)
    path.sp.m <- paste0("Mdls.", names(mcmp.l))
    path.mdls <- paste(path.res, path.sp.m, sep="/")
    pred.args <- mcmp.l[[1]]$pred.args}

  a <- c("avg.test.AUC10pct", "avg.test.AUC.MTP", "avg.test.or10pct", "avg.test.orMTP")
  b <- c("AUC (OR10p)", "AUC (ORlpt)", "OR10p (AUC)", "ORlpt (AUC)")
  sca <- c("cc26bi70", "cc45bi70", "cc60bi70", "cc85bi70", "mp26bi70", "mp45bi70", "mp85bi70", "mr26bi70", "mr45bi70", "mr60bi50", "mr85bi70", "cclgmbi", "ccmidbi", "lig_30s_bio_", "melgmbi", "memidbi", "mrlgmbi", "mrmidbi", "mxnt.preds")
  scb <- c("2070-CCSM4-rcp2.6", "2070-CCSM4-rcp4.5", "2070-CCSM4-rcp6.0", "2070-CCSM4-rcp8.5",
           "2070-MPI-ESM-LR-rcp2.6", "2070-MPI-ESM-LR-rcp4.5", "2070-MPI-ESM-LR-rcp8.5",
           "2070-MIROC-ESM-rcp2.6", "2070-MIROC-ESM-rcp4.5", "2070-MIROC-ESM-rcp6.0", "2070-MIROC-ESM-rcp8.5",
           "LGM-CCSM4", "MH-CCSM4", "LIG-CCSM3", "LGM-MPI-ESM-P", "MH-MPI-ESM-P", "LGM-MIROC-ESM", "MH-MIROC-ESM", "Present")

  t.nms <- c("fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd")
  t.NMS <- c("FCV1", "FCV5", "FCV10", "LPT (mtp)", "10P (x10ptp)", "ETSS (etss)", "MTSS", "BTO", "EETD")
  thrshld.nms.mod <- paste(c("Mod.", paste0(".",t.nms)), collapse ="|")
  thrshld.nms <- paste(t.nms, collapse = "|")

  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))


  f.plot <- function(sp, mtp.l, mdl.crit=NULL, t.NMS, t.nms, thrshld.path,
                     thrshld.nms.mod, basemap, make.underscript) {
    mSel <- mdl.crit
    mdl.crit <- grep(mdl.crit, names(mtp.l[[sp]][[1]]$binary[[1]]))

    comb.plots <- utils::combn(length(mtp.l[[sp]]), 2)
    cli.scn.pres <- which(names(mtp.l[[sp]]) == sel.clim.scn)
    sel.col <- apply(comb.plots == cli.scn.pres, 2, sum)==TRUE # rowsum(comb.plots == cli.scn.pres)==TRUE # comb.plots[, ]
    comb.plots <- matrix(comb.plots[, sel.col], nrow = 2)
    comb.plots[, comb.plots[2,] == cli.scn.pres] <- comb.plots[c(2,1), comb.plots[2,] == cli.scn.pres]

    mods.thrshld <- mtp.l[[sp]][[comb.plots[1,1]]]

    n.t <- length(mods.thrshld$binary)
    n.scn <- ncol(comb.plots)
    if(save){
      grDevices::pdf(paste(thrshld.path, paste0("Mod.clim.scn.diff.bin", ".pdf"), sep='/'),
                     width = n.scn*5+2, height = n.t*5)
    }
    graphics::par(mfcol=c(n.t, n.scn), oma = c(3.5, 0, 5.5, 2)) # , mar=c(2,4,2,5)

    for(tc in names(mods.thrshld$binary)){ # threshold criteria
      thr.CRT <- t.NMS[which(t.nms %in% tc)] #}
      cat(paste0(" - ", thr.CRT, ":", mSel))

      for(j in 1:n.scn){ # climatic scenario

        mod.sc1 <- mtp.l[[sp]][[comb.plots[1,j]]]$binary[[tc]][[mdl.crit]]
        mod.sc2 <- mtp.l[[sp]][[comb.plots[2,j]]]$binary[[tc]][[mdl.crit]]

        r.dif <- raster::overlay(mod.sc1, mod.sc2, fun=function(r1, r2) {r1-r2})

        # clim.scn name
        n1 <- names(mtp.l[[sp]])[comb.plots[1,j]]
        n2 <- names(mtp.l[[sp]])[comb.plots[2,j]]
        main.nms <- paste0("Threshold: ", thr.CRT, "\nMod. Sel. Crit: ", mSel) #, "\nClim. scen.: ", n1, " vs. ", n2)

        graphics::par(mar=c(2,4,2,4.7)) # par(mar=c(2,4,2,5))
        raster::plot(mod.sc1, breaks= c(0, .5, 1), col=c("white", "gray90"),
                     main= main.nms, legend=FALSE) #

        if(!is.null(basemap)) raster::plot(basemap, border="gray50", add= T)
        raster::plot(r.dif, breaks= c(-1, -.33, .33, 1), col=c("blue", grDevices::rgb(0,0,0,0), "red"),
                     legend=FALSE, add=TRUE) # main= main.nms,
        #
        graphics::par(mar=c(2,1,2,6)) # par(mar=c(2,3,2,6))
        raster::plot(r.dif,  legend.only=TRUE, legend.width=1.75, legend.shrink=.75,
                     xpd = TRUE, zlim=c(0, 1),
                     breaks= c(-1, -.34, .34, 1), col=c("blue", "gray90", "red"),
                     axis.args=list(at=seq(-1, 1), labels=c(n2, "equal", n1 )))
      }  # climatic scenario

    }
    if(save){
      grDevices::dev.off()
    }
  }

  for(sp in 1:length(mtp.l)){ # species
    thrshld.path <- paste(path.mdls[sp], outpt, "Mdls.thrshld", "figs", sep='/')
    if(dir.exists(thrshld.path)==FALSE) dir.create(thrshld.path)
    cat(c("\n", "Species: " , names(mtp.l)[sp]))
    mSel.crt <- mcmp.l[[sp]]$mSel
    if(all(!mSel %in% mSel.crt)){
      stop("Need to specify 'mSel' correctly. \nOptions available from selected models are: ", paste(mSel.crt, collapse = " "))
    }
    mSel <- mSel[mSel %in% mSel.crt]
    for(m in mSel){
      f.plot(sp, mtp.l, mdl.crit=m, t.NMS, t.nms, thrshld.path,
             thrshld.nms.mod, basemap, make.underscript) # comb.plots,
    }
  }
}


