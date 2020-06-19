#### 4.4 plot differences among model selection criteria predictions

#### 4.8.6 plot prediction diff between models
#' Plot differences in suitable areas between models selected using distinct criteria for multiple species
#'
#' Plot differences between predictions of models selected using distinct
#' criteria (e.g. "AvgAIC", "LowAIC", "OR", "AUC") for multiple climatic scenarios and multiple species
#'
# #' @inheritParams f.plot.mxnt.preds
#' @inheritParams thrshld_b
#' @inheritParams plot_mdl_diff
#' @param mtp.l List of stack or brick of thresholded predictions
#' @seealso \code{\link{plot_scn_diff}}
#' @description Will plot (or save as PDF) differences between predictions of models selected
#'  using distinct model selection criteria (e.g. "AvgAIC", "LowAIC", "OR", "AUC") for
#'  multiple climatic scenarios.
#' @details A panel for each combination of climatic scenario and threshold criteria will be created.
#' Within each panel, the differences between all combinations of 2 models (e.g. AvgAIC vs. OR) will
#' be plotted
#' @examples
#' \dontrun{
#' plot_mdl_diff_b(mcmp.l=mxnt.mdls.preds.lst, mtp.l=mods.thrshld.lst, basemap=NewWorld)
#' }
#' @export
plot_mdl_diff_b <- function(mcmp.l, mtp.l, basemap=NULL, save=FALSE, numCores=1,
                         msnm = c("avg.test.AUC10pct", "avg.test.AUC.MTP", "avg.test.or10pct", "avg.test.orMTP"),
                         msr = c("AUC (OR10p)", "AUC (ORlpt)", "OR10p (AUC)", "ORlpt (AUC)"),
                         scnm = c("cc26bi70", "cc45bi70", "cc60bi70", "cc85bi70", "mp26bi70", "mp45bi70", "mp85bi70",
                                  "mr26bi70", "mr45bi70", "mr60bi50", "mr85bi70", "cclgmbi", "ccmidbi", "lig_30s_bio_",
                                  "melgmbi", "memidbi", "mrlgmbi", "mrmidbi", "mxnt.preds"),
                         scr = c("2070-CCSM4-rcp2.6", "2070-CCSM4-rcp4.5", "2070-CCSM4-rcp6.0", "2070-CCSM4-rcp8.5",
                                 "2070-MPI-ESM-LR-rcp2.6", "2070-MPI-ESM-LR-rcp4.5", "2070-MPI-ESM-LR-rcp8.5",
                                 "2070-MIROC-ESM-rcp2.6", "2070-MIROC-ESM-rcp4.5", "2070-MIROC-ESM-rcp6.0",
                                 "2070-MIROC-ESM-rcp8.5", "LGM-CCSM4", "MH-CCSM4", "LIG-CCSM3", "LGM-MPI-ESM-P",
                                 "MH-MPI-ESM-P", "LGM-MIROC-ESM", "MH-MIROC-ESM", "Present")
){


  sp.nm.l <- names(mcmp.l)

  for(sp in 1:length(mtp.l)){ # species


    mcmp <- mcmp.l[[sp]]
    mtp <- mtp.l[[sp]]
    sp.nm <- sp.nm.l[sp]

    ##### function begins here:
    plot_mdl_diff(mcmp, mtp, sp.nm, basemap, save, numCores, msnm, msr, scnm, scr)

  } # fecha # species
} # fecha function




#### 4.8.6 plot prediction diff between models
#' Plot differences in suitable areas between models selected using distinct criteria
#'
#' Plot differences between predictions of models selected using distinct
#' criteria (e.g. "AvgAIC", "LowAIC", "OR", "AUC") for multiple climatic scenarios
#'
# #' @inheritParams f.plot.mxnt.preds
#' @inheritParams thrshld
#' @inheritParams calib_mdl
#' @param mtp Stack or brick of thresholded predictions
#' @param basemap Shapefile to be plotted with. Usually a continent or country shapefile
#' @param save Logical. If TRUE will save plots in pdf.
#' @param msnm Character vector. Short names of model selection criteria to be replaced. Same as given in
#' Model names (ex. AUC10, AUCmtp, OR10, ORmtp)
#' @param msr Character vector. Long names of model selection criteria to replace the short names. Must
#'  be in same order of argument 'msnm'
#' @param scnm Character vector. Short names of climatic scenarios to be replaced. Ex. "cc26bi70", "cc45bi70"
#' @param scr Character vector. Long names of climatic scenarios to replace the short names. Must
#'  be in same order of argument 'scnm'. Ex. "2070-CCSM4-rcp2.6", "2070-CCSM4-rcp4.5"
#' @seealso \code{\link{plot_scn_diff_b}}
#' @description  Will plot (or save as PDF) differences between predictions of models selected
#'  using distinct model selection criteria (for multiple climatic scenarios). PDFs will be saved
#'  within the folder 'Mdls.thrshld/figs'.
#' @examples
#' \dontrun{
#' plot_mdl_diff(mcmp=mxnt.mdls.preds.lst[[1]], mtp.l=mods.thrshld.lst[[1]], basemap=NewWorld)
#' }
#' @export
plot_mdl_diff <- function(mcmp, mtp, sp.nm="species", basemap=NULL, save=FALSE, numCores=1,
                        msnm = c("AUC10", "AUCmtp", "OR10", "ORmtp"),
                        msr = c("AUC (OR10p)", "AUC (ORlpt)", "OR10p (AUC)", "ORlpt (AUC)"),
                        scnm = c("cc26bi70", "cc45bi70", "cc60bi70", "cc85bi70", "mp26bi70", "mp45bi70", "mp85bi70",
                                 "mr26bi70", "mr45bi70", "mr60bi50", "mr85bi70", "cclgmbi", "ccmidbi", "lig_30s_bio_",
                                 "melgmbi", "memidbi", "mrlgmbi", "mrmidbi", "mxnt.preds"),
                        scr = c("2070-CCSM4-rcp2.6", "2070-CCSM4-rcp4.5", "2070-CCSM4-rcp6.0", "2070-CCSM4-rcp8.5",
                                "2070-MPI-ESM-LR-rcp2.6", "2070-MPI-ESM-LR-rcp4.5", "2070-MPI-ESM-LR-rcp8.5",
                                "2070-MIROC-ESM-rcp2.6", "2070-MIROC-ESM-rcp4.5", "2070-MIROC-ESM-rcp6.0",
                                "2070-MIROC-ESM-rcp8.5", "LGM-CCSM4", "MH-CCSM4", "LIG-CCSM3", "LGM-MPI-ESM-P",
                                "MH-MPI-ESM-P", "LGM-MIROC-ESM", "MH-MIROC-ESM", "Present")
                        ){
  {
    path.res <- "3_out.MaxEnt"
    if(dir.exists(path.res)==FALSE) dir.create(path.res)
    path.mdls <- paste(path.res, paste0("Mdls.", sp.nm), sep="/")
    thrshld.nms.mod <- paste(c("Mod\\.", paste0("\\.",tnm)), collapse ="|")
    thrshld.nms <- paste(tnm, collapse = "|")

    pred.args <- mcmp$pred.args
    outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                    ifelse(grep("logistic", pred.args)==1, 'logistic',
                           ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

    thrshld.path <- paste(path.mdls, outpt, "Mdls.thrshld", "figs", sep='/')
    if(dir.exists(thrshld.path)==FALSE) dir.create(thrshld.path)
    cat(c("\n", "Species: " , sp.nm))
  }

  if(numCores>1 & save){

    cl<-parallel::makeCluster(numCores)
    parallel::clusterExport(cl, list("plot_mdl")) # , "tnm", "tr", "auto_layout", "make_underscript"

    parallel::clusterApply(cl, names(mtp), # climatic scenario
                           function(sc, mcmp, mtp, basemap, thrshld.path, thrshld.nms.mod,
                                    save, msnm, msr, scnm, scr){ # , tr, tnm, auto_layout, make_underscript

                             plot_mdl(sc, mcmp, mtp, basemap, thrshld.path, thrshld.nms.mod,
                                       save, msnm, msr, scnm, scr) # , tr, tnm, auto_layout, make_underscript

                           }, mcmp, mtp, basemap, thrshld.path, thrshld.nms.mod,
                           save, msnm, msr, scnm, scr) # , tr, tnm, auto_layout, make_underscript

    parallel::stopCluster(cl)

  } else {

    lapply(names(mtp), # climatic scenario
         function(sc, mcmp, mtp, basemap, thrshld.path, thrshld.nms.mod,
                  save, msnm, msr, scnm, scr){ # , tr, tnm, auto_layout, make_underscript

           plot_mdl(sc, mcmp, mtp, basemap, thrshld.path, thrshld.nms.mod,
                     save, msnm, msr, scnm, scr) # , tr, tnm, auto_layout, make_underscript

         }, mcmp, mtp, basemap, thrshld.path, thrshld.nms.mod,
         save, msnm, msr, scnm, scr) # , tr, tnm, auto_layout, make_underscript
  } # fecha else

  if(save) cat(c("\n", "Figures saved in:", "\n", thrshld.path))
} # fecha function


#' internal function for \code{\link{plot_mdl_diff}}, \code{\link{plot_scn_diff_b}}
#' @inheritParams plot_mdl_diff
#' @inheritParams plot_mdl_diff_b
#' @param thrshld.path path to threshold projections
#' @param thrshld.nms.mod names of threshold models
#' @param sc Index of climatic scenario to be plotted
#' @keywords internal
plot_mdl <- function(sc, mcmp, mtp, basemap, thrshld.path, thrshld.nms.mod,
                      save, msnm, msr, scnm, scr){ # , tr, tnm, auto_layout, make_underscript
  # for(sc in names(mtp)){ # climatic scenario
  mtp.sc <- mtp[[sc]]
  cat(c("\n", "Climatic Scenario: ", sc))
  cat(c("\n", "Threshold: "))

  N.mdls <- raster::nlayers(mtp.sc$binary[[1]]) # get n models of climatic scenario
  if(N.mdls<2){
    m.nms <- names(mtp.sc$binary[[1]])
    stop("No models to compare. \nJust one model selected: ", m.nms,
         "\nModel selection criteria: ", paste(mcmp$mSel, collapse = " "), sep="")
  }
  comb.plots <- utils::combn(N.mdls, 2)

  for(l in names(mtp.sc$binary)){ # threshold criteria
    thr.crt <- l
    thr.CRT <- tr[which(tnm %in% l)] #}
    cat(paste0(" - ", thr.CRT))

    lm <- auto_layout(ncol(comb.plots), F)
    n.t <- nrow(lm)
    n.scn <- ncol(lm)
    # n.t <- length(names(mtp.sc$binary))
    # n.scn <- ncol(comb.plots)
    if(save){
      grDevices::pdf(paste(thrshld.path, paste0("Mod.diff.bin.", sc, ".", thr.crt, ".pdf"), sep='/'),
                     width = n.scn*5+2, height = n.t*5)
    }
    graphics::par(mfrow=c(n.t, n.scn), oma = c(3.5, 0, 5.5, 2))

    for(j in 1:ncol(comb.plots)){ #ncol(comb.plots)
      r.dif <- raster::overlay(mtp.sc$binary[[l]][[comb.plots[1,j]]], mtp.sc$binary[[l]][[comb.plots[2,j]]], fun=function(r1,r2) {r1-r2})
      n1 <- gsub(paste0("\\.",sc), "", gsub(thrshld.nms.mod, "", names(mtp.sc$binary[[l]][[comb.plots[1,j]]]) )  )
      n2 <- gsub(paste0("\\.",sc), "", gsub(thrshld.nms.mod, "", names(mtp.sc$binary[[l]][[comb.plots[2,j]]]) )  )
      if(n1 %in% msnm) { n1 <- msr[which(msnm %in% n1)] }
      if(n2 %in% msnm) { n2 <- msr[which(msnm %in% n2)] }

      main.nms <- paste0( make_underscript(n1), " vs. ", make_underscript(n2))
      graphics::par(mar=c(2,4,2,4.7)) # par(mar=c(2,4,2,5))
      raster::plot(mtp.sc$binary[[l]][[comb.plots[1,j]]], breaks= c(0, .5, 1), col=c("white", "gray90"),
                   legend=FALSE) # main= main.nms,
      if(j==1){
        sc1 <- gsub("mxnt.pred.fut.|mxnt.pred.past.", "", sc)
        if(sc1 %in% scnm) { sc1 <-  scr[which(scnm %in% sc1)] }
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
                   axis.args=list(at=seq(-1, 1), labels=c(make_underscript(n2), "equal", make_underscript(n1) )))
    }
    if(save){
      grDevices::dev.off()
    }
  } # fecha # threshold criteria

} # fecha # climatic scenario



####
#' Plot differences in suitable areas between climatic scenarios for multiple species
#'
#' Plot differences between a selected climatic scenario and all other climatic scenarios for each species.
#' This function will plota and (optionally) save the figures on pdf files in the folder "Mdls.thrshld/figs".
#'
#' @inheritParams thrshld_b
#' @inheritParams plot_mdl_diff
#' @inheritParams plot_scn_diff
#' @inheritParams plot_mdl_diff_b
#' @inheritParams calib_mdl
#' @seealso \code{\link{plot_scn_diff}}
#' @return won't return any object. Will save pdf's with differences among model predictions (for multiple climatic scenarios)
#' @examples
#' \dontrun{
#' plot_scn_diff_b(mcmp.l=mxnt.mdls.preds.lst, mtp.l=mods.thrshld.lst)
#' }
#' @export
plot_scn_diff_b <- function(mcmp.l, mtp.l, mSel = mcmp.l[[1]]$mSel, ref.scn="current",
                        basemap=NULL, save=FALSE, numCores=1,
                        msnm = c("AUC10", "AUCmtp", "OR10", "ORmtp"),
                        msr = c("AUC (OR10p)", "AUC (ORlpt)", "OR10p (AUC)", "ORlpt (AUC)"),
                        scnm = c("cc26bi70", "cc45bi70", "cc60bi70", "cc85bi70", "mp26bi70", "mp45bi70", "mp85bi70",
                                 "mr26bi70", "mr45bi70", "mr60bi50", "mr85bi70", "cclgmbi", "ccmidbi", "lig_30s_bio_",
                                 "melgmbi", "memidbi", "mrlgmbi", "mrmidbi", "mxnt.preds"),
                        scr = c("2070-CCSM4-rcp2.6", "2070-CCSM4-rcp4.5", "2070-CCSM4-rcp6.0", "2070-CCSM4-rcp8.5",
                                "2070-MPI-ESM-LR-rcp2.6", "2070-MPI-ESM-LR-rcp4.5", "2070-MPI-ESM-LR-rcp8.5",
                                "2070-MIROC-ESM-rcp2.6", "2070-MIROC-ESM-rcp4.5", "2070-MIROC-ESM-rcp6.0",
                                "2070-MIROC-ESM-rcp8.5", "LGM-CCSM4", "MH-CCSM4", "LIG-CCSM3", "LGM-MPI-ESM-P",
                                "MH-MPI-ESM-P", "LGM-MIROC-ESM", "MH-MIROC-ESM", "Present")
                        ){
  if(is.null(mSel)){
    stop("Need to specify 'mSel'")
  }

  sp.nm.l <- names(mcmp.l)

  for(sp in 1:length(mtp.l)){ # species


    mcmp <- mcmp.l[[sp]]
    mtp <- mtp.l[[sp]]
    sp.nm <- sp.nm.l[sp]

    ##### function begins here:
    plot_scn_diff(mcmp, mtp, mSel, ref.scn, sp.nm, basemap, save, numCores, msnm, msr, scnm, scr)

  } # fecha # species
}


#' Plot differences in suitable areas between climatic scenarios
#'
#' Plot differences between a selected climatic scenario and all other climatic scenarios for each species.
#' This function will plota and (optionally) save the figures on pdf files in the folder "Mdls.thrshld/figs".
#'
#' @inheritParams thrshld_b
#' @inheritParams plot_mdl_diff
#' @inheritParams calib_mdl
#' @param ref.scn Selected climatic scenario to compare with all others. Usually "current" one.
#' @param mSel Name of selection criteria to be compared: AvgAIC, LowAIC, avg.test.AUC10pct, avg.test.AUC.MTP,
#' avg.test.or10pct, avg.test.orMTP
#' @param save Export to pdf or not?
#' @seealso \code{\link{plot_scn_diff_b}}
#' @return won't return any object. Will save pdf's with differences among model predictions (for multiple climatic scenarios)
#' @examples
#' \dontrun{
#' plot_scn_diff(mcmp.l=mxnt.mdls.preds.lst, mtp.l=mods.thrshld.lst)
#' }
#' @export
plot_scn_diff <- function(mcmp, mtp, mSel, ref.scn="current", sp.nm="species", basemap=NULL, save=FALSE, numCores=1,
                        msnm = c("AUC10", "AUCmtp", "OR10", "ORmtp"),
                        msr = c("AUC (OR10p)", "AUC (ORlpt)", "OR10p (AUC)", "ORlpt (AUC)"),
                        scnm = c("cc26bi70", "cc45bi70", "cc60bi70", "cc85bi70", "mp26bi70", "mp45bi70", "mp85bi70",
                                 "mr26bi70", "mr45bi70", "mr60bi50", "mr85bi70", "cclgmbi", "ccmidbi", "lig_30s_bio_",
                                 "melgmbi", "memidbi", "mrlgmbi", "mrmidbi", "mxnt.preds"),
                        scr = c("2070-CCSM4-rcp2.6", "2070-CCSM4-rcp4.5", "2070-CCSM4-rcp6.0", "2070-CCSM4-rcp8.5",
                                "2070-MPI-ESM-LR-rcp2.6", "2070-MPI-ESM-LR-rcp4.5", "2070-MPI-ESM-LR-rcp8.5",
                                "2070-MIROC-ESM-rcp2.6", "2070-MIROC-ESM-rcp4.5", "2070-MIROC-ESM-rcp6.0",
                                "2070-MIROC-ESM-rcp8.5", "LGM-CCSM4", "MH-CCSM4", "LIG-CCSM3", "LGM-MPI-ESM-P",
                                "MH-MPI-ESM-P", "LGM-MIROC-ESM", "MH-MIROC-ESM", "Present")
){
  # numCores=1

  if(length(mtp)==1){
    stop("Only one climatic scenario. Nothing to compare to.")
  }

  {
    path.res <- "3_out.MaxEnt"
    if(dir.exists(path.res)==FALSE) dir.create(path.res)
    path.mdls <- paste(path.res, paste0("Mdls.", sp.nm), sep="/")
    thrshld.nms.mod <- paste(c("Mod\\.", paste0("\\.",tnm)), collapse ="|")
    thrshld.nms <- paste(tnm, collapse = "|")

    pred.args <- mcmp$pred.args
    outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                    ifelse(grep("logistic", pred.args)==1, 'logistic',
                           ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

    thrshld.path <- paste(path.mdls, outpt, "Mdls.thrshld", "figs", sep='/')
    if(dir.exists(thrshld.path)==FALSE) dir.create(thrshld.path)
    cat(c("\n", "Species: " , sp.nm))
  }


  # for(sp in 1:length(mtp.l)){ # species
  mSel.crt <- mcmp$mSel
  if(all(!mSel %in% mSel.crt)){
    stop("Need to specify 'mSel' correctly. \nOptions available from selected models are: ", paste(mSel.crt, collapse = " "))
  }
  mSel <- mSel[mSel %in% mSel.crt]
  # m <- mSel[4]
  if(any(grepl("^AUC", mSel))){
    if(sum(grepl("AUCmtp|AUC10", names(mtp[[1]]$binary[[1]]))) == 1){
      mSel[grepl("^AUC", mSel)] <- "AUCmtp.AUC10"
    } else {
    mSel <- c(mSel[!grepl("^AUC", mSel)], c("AUCmtp", "AUC10"))
    }
    # mSel <- c(mSel[!grepl("^AUC", mSel)], c("^AUC"))
  }
  if(any(grepl("^OR", mSel))){
    mSel <- c(mSel[!grepl("^OR", mSel)], c("ORmtp", "OR10"))
  }

  if(numCores>1 & save){

    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl, list("plot_scn"))# , "tnm", "tr", "auto_layout", "make_underscript"))
    # , plot_scn, tnm, tr, auto_layout, make_underscript

    parallel::clusterApply(cl, mSel,  # model selection criteria
                           function(m, mtp, ref.scn, basemap, thrshld.path, thrshld.nms.mod,
                                    save, msnm, msr, scnm, scr){ # , tr, tnm, auto_layout, make_underscript

                             plot_scn(m, mtp, ref.scn, basemap, thrshld.path, thrshld.nms.mod,
                                       save, msnm, msr, scnm, scr) # , tr, tnm, auto_layout, make_underscript

                           }, mtp, ref.scn, basemap, thrshld.path, thrshld.nms.mod,
                           save, msnm, msr, scnm, scr) # , tr, tnm, auto_layout, make_underscript

    parallel::stopCluster(cl)

  } else {

    lapply(mSel, # model selection criteria
           function(m, mtp, ref.scn, basemap, thrshld.path, thrshld.nms.mod,
                    save, msnm, msr, scnm, scr){ # , tr, tnm, auto_layout, make_underscript

             plot_scn(m, mtp, ref.scn, basemap, thrshld.path, thrshld.nms.mod,
                       save, msnm, msr, scnm, scr) # , tr, tnm, auto_layout, make_underscript

           }, mtp, ref.scn, basemap, thrshld.path, thrshld.nms.mod,
           save, msnm, msr, scnm, scr) # , tr, tnm, auto_layout, make_underscript
  }

  if(save) cat(c("\n", "Figures saved in:", "\n", thrshld.path))
}

#' internal function for \code{\link{plot_mdl_diff}}, \code{\link{plot_scn_diff_b}}
#' @inheritParams plot_mdl
#' @inheritParams plot_scn_diff
#' @inheritParams plot_mdl_diff_b
#' @param m Index of model selection criteria to be plotted
#' @keywords internal
plot_scn <- function(m, mtp, ref.scn="current", basemap, thrshld.path, thrshld.nms.mod,
                      save, msnm, msr, scnm, scr){ # , tr, tnm, auto_layout, make_underscript

  # for(m in mSel){ # model selection criteria
    # plot_scn(m, sp, mtp.l, mdl.crit=m, tr, tnm, thrshld.path,
    #           thrshld.nms.mod, basemap, make_underscript, save, msnm, msr, scnm, scr) # comb.plots,
    # mSel <- mdl.crit
    # mdl.crit <- grep(m, names(mtp[[1]]$binary[[1]]))
  mdl.crit <- grep(m, names(mtp[[1]]$binary[[1]]))

  comb.plots <- utils::combn(length(mtp), 2)
  cli.scn.pres <- which(names(mtp) == ref.scn)
  sel.col <- apply(comb.plots == cli.scn.pres, 2, sum)==TRUE # rowsum(comb.plots == cli.scn.pres)==TRUE # comb.plots[, ]
  comb.plots <- matrix(comb.plots[, sel.col], nrow = 2)
  comb.plots[, comb.plots[2,] == cli.scn.pres] <- comb.plots[c(2,1), comb.plots[2,] == cli.scn.pres]

    mods.thrshld <- mtp[[comb.plots[1,1]]]

    n.t <- length(mods.thrshld$binary)
    n.scn <- ncol(comb.plots)
    # lm <- auto_layout(ncol(comb.plots), F)
    lm <- auto_layout((n.t*n.scn), F)
    n.r <- nrow(lm)
    n.col <- ncol(lm)
    if(save){
      grDevices::pdf(paste(thrshld.path, paste0("Mod.clim.scn.diff.bin", m, ".pdf"), sep='/'),
                     width = n.col*5+2, height = n.r*5)
    }
    graphics::par(mfrow=c(n.col, n.r), mar=c(2,4,2,5), oma = c(3.5, 0, 3.5, 2)) #

    for(tc in names(mods.thrshld$binary)){ # threshold criteria
      thr.CRT <- tr[which(tnm %in% tc)] #}
      cat(paste0(" - ", thr.CRT, ":", m))

      for(j in 1:n.scn){ # climatic scenario

        mod.sc1 <- mtp[[comb.plots[1,j]]]$binary[[tc]][[mdl.crit]]
        mod.sc2 <- mtp[[comb.plots[2,j]]]$binary[[tc]][[mdl.crit]]

        r.dif <- raster::overlay(mod.sc1, mod.sc2, fun=function(r1, r2) {r1-r2})

        # clim.scn name
        n1 <- names(mtp)[comb.plots[1,j]]
        n2 <- names(mtp)[comb.plots[2,j]]
        main.nms <- paste0("Threshold: ", thr.CRT, "\nMod. Sel. Crit: ", m) #, "\nClim. scen.: ", n1, " vs. ", n2)

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
  # }
}


############ TODO
#' Plot compare climatic scenarios differences between model selection criteria
#'
#' i.e. plot_scn_diff with panels for each of thresholds, and climatic scenario
#'  difference (i.e. current vs. climScn1) for each model being compared side




##



##### internal functions
#### 4.4 plot differences among model selection criteria predictions
#' Format (underscript) selected texts for plotting
#'
#' Format (underscript) selected texts (criteria used to select models) to be used on plotting.
#'
#' @param x list of text to be formatted
#' @seealso \code{\link{plot_mdl_diff}}, \code{\link{plot_scn_diff_b}}
#' @return list of formatted text
# #' @examples
# #' make_underscript(c("AUC (OR10p)", "AUC (ORlpt)", "OR10p (AUC)", "ORlpt (AUC)"))
#' @keywords internal
# #' @export
make_underscript <- function(x) as.expression(lapply(x, function(y) {
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

#' Determine the arrangement of multiple plots in a single panel
#' Given a particular number of plots, \code{auto_layout} will automatically determine the arrangement of each
#' plot using the \code{layout} function or par(mfrow=c(nrow, ncol)). See examples.
#' modified from: https://github.com/cran/fifer/blob/master/R/auto.layout.R
#' @title Automatically select the layout.
#' @param n the number of plots
#' @param layout should the fuction return a preallocated layout object? If \code{FALSE}, it returns a matrix
#' @return either a matrix or a layout object
#' @author Dustin Fife
#' @examples
#' \dontrun{
#' ## plot six plots
#' auto_layout(6)
#' for (i in 1:6){
#' 	plot(rnorm(100), rnorm(100))
#' }
#' ## same as mar(mfrow=c(3,2))
#' par(mfrow=c(3,2))
#' for (i in 1:6){
#' 	plot(rnorm(100), rnorm(100))
#' }
#' ## default for odd number of plots using mfrow looks terrible
#' par(mfrow=c(3,2))
#' for (i in 1:5){
#' 	plot(rnorm(100), rnorm(100))
#' }
#' ## much better with auto_layout
#' auto_layout(5)
#' for (i in 1:5){
#' 	plot(rnorm(100), rnorm(100))
#' }
#' ## see matrices of layouts for multiple plots
# for(i in 4:12){
#   print(auto_layout(i, layout=F))
# }
#' ##
#' for(i in 2:6){
#'   m <- auto_layout(i, layout=F)
#'   par(mfrow=c(nrow(m), ncol(m)))
#'   for (j in 1:i){
#'     plot(rnorm(100), rnorm(100))
#'   }
#'   Sys.sleep(1)
#' }
#' }
# #' @export
#' @keywords internal
auto_layout <- function(n, layout=T){
  ### figure out how many rows
  sq = sqrt(n)
  rws = round(sq)

  #### if it's a perfect square, fill the matrix
  if (sq == rws){
    numbs = sort(rep(1:n, times=1))
    m = matrix(numbs, nrow=sq, byrow=T)
  } else {

    #### repeat twice the numbers that fit nicely
    topNum = trunc(n/rws)*rws
    numbs = sort(rep(1:topNum, times=1))
    if (topNum==n){
      m = matrix(numbs, nrow=rws, byrow=T)
    } else {
      #### get the rest figured out
      # rem = n-topNum  ### remaining numbers
      rest = (topNum+1):n # sort(rep((topNum+1):n, times=1))
      cols = (topNum/rws)
      if(cols<rws){
        top <- cols*(rws+1)
        rest = c(rest, rep(0, times=(top-n))) # rep(0, times=(cols-length(rest))/2),
        m = matrix(c(numbs, rest), ncol = cols , nrow=rws+1, byrow=T)
      } else {
        top <- (cols+1)*(rws)
        rest = c(rest, rep(0, times=(top-n))) # rep(0, times=(cols-length(rest))/2),
        m = matrix(c(numbs, rest), ncol = cols+1, nrow=rws, byrow=T)
      }
      # rest = c(rest, rep(0, times=(top-n))) # rep(0, times=(cols-length(rest))/2),
      # m = matrix(c(numbs, rest), ncol = cols+1 , nrow=rws, byrow=T)
    }
  }

  if (layout){
    layout(m)
  } else {
    m
  }
}



