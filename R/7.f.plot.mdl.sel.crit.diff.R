#### 4.4 plot differences among model selection criteria predictions


#' Format (underscript) selected texts for plotting
#'
#' General function description. A short paragraph (or more) describing what the function does.
#' @param x list of text to be formatted
#' @return list of formatted text
#' @examples
#' make.underscript(c("AUC (OR10p)", "AUC (ORlpt)", "OR10p (AUC)", "ORlpt (AUC)"))
#' @export
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


#### 4.4 plot differences among model selection criteria predictions
#' Plot differences among model predictions selected from several criteria
#'
#' General function description. A short paragraph (or more) describing what the function does.
#' @inheritParams f.thr.batch
#' @param mtp.spl List of stack or brick of thresholded predictions
#' @param basemap Shapefile to be plotted with. Usually a continent or country shapefile
#' @return won't return any object. Will save pdf's with differences among model predictions
#' @examples
#' f.plot.mxnt.preds(mxnt.mdls.preds.lst, mods.thrshld.lst, basemap=NewWorld)
#' @export
f.plot.mxnt.preds <- function(mmp.spl, mtp.spl, basemap=NULL, pred.nm=""){
  { path.res <- "4_ENMeval.results"
  if(dir.exists(path.res)==F) dir.create(path.res)
  path.sp.m <- paste0("Mdls.", names(mmp.spl))
  path.mdls <- paste(path.res, path.sp.m, sep="/")
  pred.args <- mmp.spl[[1]]$pred.args}

  comb.plots <- utils::combn(raster::nlayers(mtp.spl[[1]]$continuous[[2]]), 2)

  # unlist(strsplit(names(mods.thrshld$binary[[2]]), "."))

  thrshld.nms <- c("fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd")
  thrshld.nms.mod <- paste(c("Mod.", paste0(".",thrshld.nms)), collapse ="|")
  thrshld.nms <- paste(thrshld.nms, collapse = "|")

  # ncol(comb.plots)

  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))
  pred.nm <- ifelse(pred.nm != "", paste0(".", pred.nm), pred.nm)
  for(i in 1:length(mtp.spl)){ # species
    thrshld.path <- paste(path.mdls[i], outpt, "Mdls.thrshld", "figs", sep='/')
    if(dir.exists(thrshld.path)==F) dir.create(thrshld.path)

    mods.thrshld <- mtp.spl[[i]]

    for(l in 1:length(mods.thrshld$binary)){ # threshold criteria
      thr.crt <- grep(thrshld.nms,  unique(unlist(strsplit(names(mods.thrshld$binary[[l]]), "."))), value=T)
      grDevices::pdf(paste(thrshld.path, paste0("Mod.diff.bin", pred.nm, ".", thr.crt, ".pdf"), sep='/'),
                     width = 20, height = 10)
      graphics::par(mfrow=c(3,5), mar=c(2,1,2,1), oma = c(1, 1, 4, 1))

      for(j in 1:ncol(comb.plots)){ #ncol(comb.plots)
        # r.dif <- mods.thrshld$binary[[l]][[comb.plots[1,j]]] - mods.thrshld$binary[[l]][[comb.plots[2,j]]]#, col=c("red", "white", "blue")
        r.dif <- raster::overlay(mods.thrshld$binary[[l]][[comb.plots[1,j]]], mods.thrshld$binary[[l]][[comb.plots[2,j]]], fun=function(r1,r2) {r1-r2})
        n1 <- sub("_", ".", sub("min", "Min", sub("Mean", "M", sub("OR", "or", gsub(thrshld.nms.mod, "", names(mods.thrshld$binary[[l]][[comb.plots[1,j]]]))))))
        n2 <- sub("_", ".", sub("min", "Min", sub("Mean", "M", sub("OR", "or", gsub(thrshld.nms.mod, "", names(mods.thrshld$binary[[l]][[comb.plots[2,j]]]))))))
        main.nms <- paste0(n1, " vs. ", n2)
        # r.dif, breaks= c(-1, -.33, .33, 1), col=c("blue", "white", "red")
        raster::plot(r.dif, breaks= c(-1, -.33, .33, 1), col=c("blue", "white", "red"),# col=c("blue", "blue","white", "red", "red"), # smallplot=c(0.01, .81, 0.01, .99),
                     main= main.nms, legend=FALSE)
        if(j==1){
          graphics::title(paste("Threshold criteria:", thr.crt), line = 2, outer = T, cex.main=2)
        }
        raster::plot(r.dif,  legend.only=TRUE, smallplot=c(.78, .79, .2, .8),  #horiz=T,#  c(.79, .80, .2, .8)
                     breaks= c(-1, -.34, .34, 1), col=c("blue", "white", "red"),
                     # bg = "white",
                     #col=c("blue", "blue","white", "red", "red"),
                     axis.args=list(at=seq(-1, 1),
                                    labels=c(n2, "equal", n1)))
        if(!is.null(basemap)) sp::plot(basemap, add= T)
      }

      grDevices::dev.off()
    }

  }

}


#### 4.8.6 plot prediction diff between models
#' Plot differences (for multiple climatic scenarios) among model predictions selected from several criteria
#'
#' Plot differences (for multiple climatic scenarios) among model predictions selected from several criteria
#' @inheritParams f.plot.mxnt.preds
#' @inheritParams mxnt.cp
#' @return won't return any object. Will save pdf's with differences among model predictions (for multiple climatic scenarios)
#' @examples
#' f.plot.mxnt.preds.mscn(mxnt.mdls.preds.lst, mods.thrshld.lst, basemap=NewWorld)
#' f.plot.mxnt.preds.mscn(mxnt.mdls.preds.pf[1], mods.thrshld.lst[1], basemap=NewWorld)
#' @export
f.plot.mxnt.preds.mscn <- function(mmp.spl, mtp.spl, basemap=NULL, numCores=1){
  { path.res <- "4_ENMeval.results"
  if(dir.exists(path.res)==F) dir.create(path.res)
  path.sp.m <- paste0("Mdls.", names(mmp.spl))
  path.mdls <- paste(path.res, path.sp.m, sep="/")
  pred.args <- mmp.spl[[1]]$pred.args}

  a <- c("Mean.AUC10", "Mean.AUCmin", "Mean.OR10", "Mean.ORmin")
  b <- c("AUC (OR10p)", "AUC (ORlpt)", "OR10p (AUC)", "ORlpt (AUC)")
  sca <- c("cc26bi70", "cc45bi70", "cc60bi70", "cc85bi70", "mp26bi70", "mp45bi70", "mp85bi70", "mr26bi70", "mr45bi70", "mr60bi50", "mr85bi70", "cclgmbi", "ccmidbi", "lig_30s_bio_", "melgmbi", "memidbi", "mrlgmbi", "mrmidbi", "mxnt.preds")
  scb <- c("2070-CCSM4-rcp2.6", "2070-CCSM4-rcp4.5", "2070-CCSM4-rcp6.0", "2070-CCSM4-rcp8.5",
           "2070-MPI-ESM-LR-rcp2.6", "2070-MPI-ESM-LR-rcp4.5", "2070-MPI-ESM-LR-rcp8.5",
           "2070-MIROC-ESM-rcp2.6", "2070-MIROC-ESM-rcp4.5", "2070-MIROC-ESM-rcp6.0", "2070-MIROC-ESM-rcp8.5",
           "LGM-CCSM4", "MH-CCSM4", "LIG-CCSM3", "LGM-MPI-ESM-P", "MH-MPI-ESM-P", "LGM-MIROC-ESM", "MH-MIROC-ESM", "Present")

  comb.plots <- utils::combn(raster::nlayers(mtp.spl[[1]][[1]]$binary[[1]]), 2)

  # unlist(strsplit(names(mods.thrshld$binary[[2]]), "."))

  t.nms <- c("fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd")
  t.NMS <- c("FCV1", "FCV5", "FCV10", "LPT (mtp)", "10P (x10ptp)", "ETSS (etss)", "MTSS", "BTO", "EETD")
  thrshld.nms.mod <- paste(c("Mod.", paste0(".",t.nms)), collapse ="|")
  thrshld.nms <- paste(t.nms, collapse = "|")

  # ncol(comb.plots)

  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))
  # pred.nm <- ifelse(pred.nm != "", paste0(".", pred.nm), pred.nm)
  # cat(c("\n"))
  for(sp in 1:length(mtp.spl)){ # species
    thrshld.path <- paste(path.mdls[sp], outpt, "Mdls.thrshld", "figs", sep='/')
    if(dir.exists(thrshld.path)==F) dir.create(thrshld.path)
    cat(c("\n", "Species: " , names(mtp.spl)[sp]))
    # cat(c("\n", "Climatic Scenario:"))
    ## TODO use mclapply

    f.plot <- function(sc, sp, mtp.spl, mods.thrshld, t.NMS, t.nms, thrshld.path,
                       comb.plots, thrshld.nms.mod, basemap, make.underscript) {
      # for(sc in names(mtp.spl[[sp]])){ # climatic scenario
      mods.thrshld <- mtp.spl[[sp]][[sc]]
      cat(c("\n", "Climatic Scenario: ", sc))
      cat(c("\n", "Threshold Criterion: "))

      for(l in names(mods.thrshld$binary)){ # threshold criteria
        thr.crt <- l #sub("x10ptp", "10P", sub("mtp", "LPT", l))  #grep(thrshld.nms,  unique(unlist(strsplit(names(mods.thrshld$binary[[l]]), "."))), value=T)
        #### TODO
        # if(l %in% ta) {
        thr.CRT <- t.NMS[which(t.nms %in% l)] #}
        cat(paste0(" - ", thr.CRT)) # cat(paste0(thr.CRT))

        grDevices::pdf(paste(thrshld.path, paste0("Mod.diff.bin.", sc, ".", thr.crt, ".pdf"), sep='/'),
                       width = 20, height = 10)
        # par(mfrow=c(3,5), mar=c(2,1,2,1), oma = c(1, 1, 4, 1))
        graphics::par(mfrow=c(3,5), mar=c(2,4,2,5), oma = c(3.5, 0, 5.5, 2))

        for(j in 1:ncol(comb.plots)){ #ncol(comb.plots)
          # r.dif <- mods.thrshld$binary[[l]][[comb.plots[1,j]]] - mods.thrshld$binary[[l]][[comb.plots[2,j]]]#, col=c("red", "white", "blue")
          r.dif <- raster::overlay(mods.thrshld$binary[[l]][[comb.plots[1,j]]], mods.thrshld$binary[[l]][[comb.plots[2,j]]], fun=function(r1,r2) {r1-r2})
          # n1 <- gsub(paste0(".",sc), "", sub("_", ".", sub("min", "Min", sub("Mean", "M", sub("OR", "or", gsub(thrshld.nms.mod, "", names(mods.thrshld$binary[[l]][[comb.plots[1,j]]]))))))) # ".mxnt.preds|.mxnt.pred|mxnt.preds|mxnt.pred"
          # n2 <- gsub(paste0(".",sc), "", sub("_", ".", sub("min", "Min", sub("Mean", "M", sub("OR", "or", gsub(thrshld.nms.mod, "", names(mods.thrshld$binary[[l]][[comb.plots[2,j]]])))))))
          n1 <- gsub(paste0(".",sc), "", gsub(thrshld.nms.mod, "", names(mods.thrshld$binary[[l]][[comb.plots[1,j]]]))  )
          n2 <- gsub(paste0(".",sc), "", gsub(thrshld.nms.mod, "", names(mods.thrshld$binary[[l]][[comb.plots[2,j]]]))  )
          if(n1 %in% a) { n1 <- b[which(a %in% n1)] }
          if(n2 %in% a) { n2 <- b[which(a %in% n2)] }

          main.nms <- paste0( make.underscript(n1), " vs. ", make.underscript(n2))
          # r.dif, breaks= c(-1, -.33, .33, 1), col=c("blue", "white", "red")
          graphics::par(mar=c(2,4,2,4.7)) # par(mar=c(2,4,2,5))
          raster::plot(mods.thrshld$binary[[l]][[comb.plots[1,j]]], breaks= c(0, .5, 1), col=c("white", "gray90"),
                       legend=FALSE) # main= main.nms,
          if(j==1){
            # title(paste0("Climatic scenario: ", sub("mxnt.pred.", "", sub("mxnt.preds", "present", sc)), ".  Threshold criteria: ", thr.crt), line = 2, outer = T, cex.main=2)

            sc1 <- gsub("mxnt.pred.fut.|mxnt.pred.past.", "", sc)
            if(sc1 %in% sca) { sc1 <-  scb[which(sca %in% sc1)] } # {paste0(scb[which(sca %in% sc1)], " (", sc, ")")} else {sc1 <- sc1}
            # title(paste0("Species: " , names(mtp.spl)[sp]), line = 3.5, outer = T, cex.main=2)
            graphics::title(paste0("Climatic Scenario: ", sub("mxnt.pred.", "", sub("mxnt.preds", "bioclim", sc1))), line = 1.5, outer = T, cex.main=2)
            graphics::title(paste0("Threshold Criterion: ", thr.CRT), line = -.5, outer = T, cex.main=2)
          }
          if(!is.null(basemap)) raster::plot(basemap, border="gray50", add= T)
          raster::plot(r.dif, breaks= c(-1, -.33, .33, 1), col=c("blue", grDevices::rgb(0,0,0,0), "red"),
                       legend=FALSE, add=T) # main= main.nms,

          graphics::par(mar=c(2,1,2,6)) # par(mar=c(2,3,2,6))
          raster::plot(r.dif,  legend.only=TRUE, legend.width=1.75, legend.shrink=.75, #smallplot=c(.78, .79, .2, .8),  #horiz=T,#  c(.79, .80, .2, .8)
                       xpd = TRUE, zlim=c(0, 1),#legend.args=list(side=4),
                       breaks= c(-1, -.34, .34, 1), col=c("blue", "gray90", "red"),
                       axis.args=list(at=seq(-1, 1), labels=c(make.underscript(n2), "equal", make.underscript(n1) )))
        }
        grDevices::dev.off()
      }
      # }
    }

    if(numCores>1){

      cl<-parallel::makeCluster(numCores)

      parallel::clusterApply(cl, names(mtp.spl[[sp]]), # climatic scenario
                             function(sc, sp, mtp.spl, mods.thrshld, t.NMS, t.nms, thrshld.path,
                                      comb.plots, thrshld.nms.mod, basemap, make.underscript){

                               f.plot(sc, sp, mtp.spl, mods.thrshld, t.NMS, t.nms, thrshld.path,
                                      comb.plots, thrshld.nms.mod, basemap, make.underscript)

                             }, sp, mtp.spl, mods.thrshld, t.NMS, t.nms, thrshld.path,
                             comb.plots, thrshld.nms.mod, basemap, make.underscript)

      parallel::stopCluster(cl)

    }else{
      lapply(names(mtp.spl[[sp]]), # climatic scenario
             function(sc, sp, mtp.spl, mods.thrshld, t.NMS, t.nms, thrshld.path,
                      comb.plots, thrshld.nms.mod, basemap, make.underscript){

               f.plot(sc, sp, mtp.spl, mods.thrshld, t.NMS, t.nms, thrshld.path,
                      comb.plots, thrshld.nms.mod, basemap)

             }, sp, mtp.spl, mods.thrshld, t.NMS, t.nms, thrshld.path,
             comb.plots, thrshld.nms.mod, basemap, make.underscript)

    }

  }
}



