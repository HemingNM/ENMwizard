# # ##### 5. METRICS

# TO DO - get_tsa;
# In area.occ.spp[[sp]][] <- array(aperm(ar.mods.t.p, c(3, 2, 1))) :
#   number of items to replace is not a multiple of replacement length
# TO DO - get_fpa
# Error in `[<-.data.frame`(`*tmp*`, , ncol(areas), value = c(0.526, 0.461,  :
#   replacement has 6 rows, data has 8

# # #### 4.6 compute area occupied
# #' Compute total suitable area
# #'
# #' General function description. A short paragraph (or more) describing what the function does.
# #' @inheritParams f.plot.mxnt.preds
# #' @param pred.nm name of prediction to be appended to the final name. Usually "pres", "past" or "fut".
# #' @param thrshld.i List of threshold criteria to be applied
# #' @return Stack or brick of thresholded predictions
# #' @examples
# #' areas.occ.lst <- f.area.occ(mtp.l)
# #' @export


# #### 5.1 compute area occupied at multiple scenarios

#' Compute species' total suitable area
#'
#' Compute total suitable area at multiple climatic scenario, threshold and model criteria.
#'
#' @inheritParams plot_mdl_diff_b
#' @inheritParams get_tsa
#' @seealso \code{\link[raster]{area}}, \code{\link{get_tsa}}, \code{\link{get_cont_permimport}},
#' \code{\link{get_fpa}}, \code{\link{get_cont_permimport_b}}, \code{\link{get_fpa_b}}
#' @return List of arrays containing species' total suitable areas for each climatic scenario, threshold and model criteria
#' @examples
#' \dontrun{
#' areas.occ.lst <- get_tsa_b(mtp.l=mods.thrshld.lst)
#' }
#' @export
get_tsa_b <- function(mtp.l, restrict=NULL, area.raster=NULL, digits=0){

  area.occ.spp <- lapply(seq_along(mtp.l), function(i, mtp.l, restrict, area.raster, digits){
    get_tsa(mtp.l[[i]], restrict, area.raster, digits)
  }, mtp.l, restrict, area.raster, digits) # species, areas

  names(area.occ.spp) <- names(mtp.l)

  area.occ.spp.c <- data.table::rbindlist(area.occ.spp, idcol = "sp")
  utils::write.csv(area.occ.spp.c, paste0("3_out.MaxEnt/metric.totalArea.csv")) # reorder ds

  return(area.occ.spp.c)
}


#' Compute species' total suitable area
#'
#' Compute total suitable area at multiple climatic scenario, threshold and model criteria.
#'
#' @inheritParams plot_mdl_diff
#' @param digits integer indicating the number of decimal places. see ?round for details.
#' @param restrict A raster to select a region to compute area.
#' @param area.raster A raster containing the cell areas to be summed across
#' the suitable pixels. This allows summing areas of habitat when the pixel is
#' partially occupied with the habitat of interest.
#' @seealso \code{\link[raster]{area}}, \code{\link{get_tsa_b}}, \code{\link{get_cont_permimport}}, \code{\link{get_fpa}},
#' \code{\link{get_cont_permimport_b}}, \code{\link{get_fpa_b}}
#' @return List of arrays containing species' total suitable areas for each climatic scenario, threshold and model criteria
#' @examples
#' \dontrun{
#' areas.occ.lst <- get_tsa_b(mtp.l=mods.thrshld.lst)
#' }
#' @export
get_tsa <- function(mtp, restrict = NULL, area.raster=NULL, digits){ # species, areas
  thrshld.nms <- paste(paste0(".", tnm), collapse = "|")

  c.nms <- gsub(paste0("Mod\\.|", gsub("\\.", "\\\\.", thrshld.nms)), "", names(mtp[[1]][[2]][[1]]))
  c.nms2 <- vector("character", length(c.nms))
  s.nms <- c("LowAIC", "ORmtp", "OR10", "AUCmtp", "AUC10", "^AvgAIC", "^EBPM", "^WAAUC", "^ESORIC")
  invisible(sapply(seq_along(s.nms), function(i, x, y, z){
    si <- grepl(s.nms[i], c.nms)
    if(sum(si)>0){
      c.nms2[si] <<- gsub("\\^|^\\.", "", paste(c.nms2[si], s.nms[i], sep = "."))
    }
  }, c.nms, s.nms, c.nms2))
  rep.nm <- find_repeated_characters(gsub(paste(unique(c.nms2), collapse = "|"), "", c.nms))
  masks <- gsub(paste(c(rep.nm,paste(unique(c.nms2), collapse = "|")), collapse = "|"), "", c.nms)
  masks[masks==""] <- "all"
  rep.mdl <- length(c.nms2)/length(unique(c.nms2))
  c.nms <- c.nms2

  thrshld.crit <- names(mtp[[1]][[1]])

  ar.mods.t.p <- lapply(seq_along(mtp), function(sc, mtp, restrict, area.raster, digits){ # , areas  # pred.scenario
    mtp.sc <- mtp[[sc]][[2]]

    ar.mods.t <- sapply(seq_along(mtp.sc), function(t, mtp.sc, sc, restrict, area.raster, digits){ # , areas # threshold criteria
      mtp.sc.t <- mtp.sc[[t]]

      ar.mods <- sapply(1:raster::nlayers(mtp.sc.t), function(m, mtp.sc.t, sc, t, restrict, area.raster, digits){ # , areas # model criteria
        ar <- mtp.sc.t[[m]]

        if(grDevices::is.raster(restrict)){
          if(raster::res(ar)!=raster::res(restrict)){
            ar <- raster::resample(ar, restrict)
            ar <- ar*restrict
          }
        }
        if(is.null(area.raster)){
          area.raster <- raster::area(ar, na.rm=TRUE)
        }
        if(isTRUE(all.equal(raster::extent(ar), raster::extent(area.raster)))){
          area.raster <- raster::crop(area.raster, ar)
        }
        ar <- raster::zonal(area.raster, ar, "sum", digits=digits)
        ar <- empty2zero(ar[ar[,1]==1, 2])
        return(ar) }, mtp.sc.t, sc, t, restrict, area.raster, digits) # , areas # model criteria
      return(ar.mods) }, mtp.sc, sc, restrict, area.raster, digits) # , areas# threshold criteria
    return(ar.mods.t) }, mtp, restrict, area.raster, digits) # , areas # pred.scenario

  ar.mods.t.p <- simplify2array(ar.mods.t.p) # transform list into array
  if(length(dim(ar.mods.t.p))==3){
    ar.mods.t.p <- array(aperm(ar.mods.t.p, c(3,2,1))) #,
  } else if(length(dim(ar.mods.t.p))==2){
    dim(ar.mods.t.p) <- c(dim(ar.mods.t.p), 1)
    ar.mods.t.p <- array(aperm(ar.mods.t.p, c(3,2,1))) #,
  } else if(length(dim(ar.mods.t.p))==1){
    dim(ar.mods.t.p) <- c(dim(ar.mods.t.p), 1, 1)
    ar.mods.t.p <- array(aperm(ar.mods.t.p, c(3,2,1))) #,
  }

  # https://stackoverflow.com/questions/40921426/converting-array-to-matrix-in-r
  areas <- data.frame(expand.grid(Clim.scen=names(mtp), # pred.scenario
                                  threshold=names(mtp[[1]][[2]]), # threshold criteria
                                  Model=c.nms), # model criteria
                      Location=rep(unique(masks), each=length(ar.mods.t.p)/rep.mdl),
                      TotSuitArea=ar.mods.t.p)

  return(areas)
}


# #### 4.7 extract model results
# ### 4.7.1 variable contribution and importance

#' Compute variable contribution and permutation importance
#'
#' Compute variable contribution and importance for each model
#'
# #' @param mcmp.l Stack or brick of predictions to apply the threshold
#' @inheritParams thrshld_b
#' @inheritParams get_tsa_b
#' @seealso \code{\link{get_cont_permimport}}, \code{\link{get_tsa}}, \code{\link{get_fpa}},
#' \code{\link{get_tsa_b}}, \code{\link{get_fpa_b}}, \code{\link[dismo]{maxent}}
#' @return List of arrays containing variable contribution and importance for each species
#' @examples
#' \dontrun{
#' get_cont_permimport_b(mcmp.l = mxnt.mdls.preds.lst)
#' }
#' @export
get_cont_permimport_b <- function(mcmp.l){
  path.res <- "3_out.MaxEnt"
  if(dir.exists(path.res)==FALSE) dir.create(path.res)

  # var.contPermImp <- stats::setNames(vector("list", length(mcmp.l)), names(mcmp.l))

  var.contPermImp <- lapply(seq_along(mcmp.l), function(i, mcmp.l){
    get_cont_permimport(mcmp.l[[i]], names(mcmp.l[i]))
  }, mcmp.l) # species, areas

  names(var.contPermImp) <- names(mcmp.l)

  var.cont.sp <- data.table::rbindlist(lapply(var.contPermImp, function(x) x[[1]]), idcol = "sp", fill=T)
  utils::write.csv(var.cont.sp, paste0("3_out.MaxEnt/metric.var.Contribution.csv")) # reorder ds

  var.permImp.sp <- data.table::rbindlist(lapply(var.contPermImp, function(x) x[[2]]), idcol = "sp", fill=T)
  utils::write.csv(var.permImp.sp, paste0("3_out.MaxEnt/metric.var.PermImportance.csv")) # reorder ds

  # var.contPermImp.c <- data.table::rbindlist(var.contPermImp[[1]], idcol = "sp")
  # colnames(area.occ.spp.c)[1:5] <- c("sp", "Clim.scen", "threshold", "Model", "TotSuitArea")

  return(var.contPermImp)
}


#' Compute variable contribution and permutation importance
#'
#' Compute variable contribution and importance for each model
#'
# #' @param mcmp.l Stack or brick of predictions to apply the threshold
#' @inheritParams thrshld
#' @inheritParams get_tsa
#' @seealso \code{\link{get_cont_permimport_b}}, \code{\link{get_tsa}}, \code{\link{get_fpa}},
#' \code{\link{get_tsa_b}}, \code{\link{get_fpa_b}}, \code{\link[dismo]{maxent}}
#' @return List of arrays containing variable contribution and importance for each species
#' @examples
#' \dontrun{
#' get_cont_permimport(mcmp = mxnt.mdls.preds)
#' }
#' @export
get_cont_permimport <- function(mcmp, sp.nm) {
  mxnt.mdls <- mcmp$mxnt.mdls
  sel.mod.nms <- paste0("Mod.", mcmp$selected.mdls$sel.cri)
  mod.nms <- paste0("Mod_", format(mcmp$selected.mdls[, "rm"], digits=2), "_", mcmp$selected.mdls[, "features"]) #
  # mod.nms <- paste0("Mod.", mcmp$selected.mdls$settings)
  pred.nms <- names(mcmp$mxnt.preds[[1]])
  var.nms <- gsub( ".contribution", "", rownames(mxnt.mdls[[1]]@results)[grepl("contribution", rownames(mxnt.mdls[[1]]@results))])
  # w.mdls <- mcmp$selected.mdls$w.AIC
  if(sum(grepl("AvgAIC", pred.nms))>0) {
    wv.aic <- mcmp[["selected.mdls"]][grep("AIC_", mcmp[["selected.mdls"]]$sel.cri),"w.AIC"]
  }
  if(sum(grepl("WAAUC", pred.nms))>0) {
    wv.wa <- mcmp[["selected.mdls"]][grep("WAAUC_", mcmp[["selected.mdls"]]$sel.cri),"avg.test.AUC"]
  }
  if(sum(grepl("EBPM", pred.nms))>0) {
    wv.bp <- rep(1, length(grep("EBPM", mcmp[["selected.mdls"]]$sel.cri)))
  }
  if(sum(grepl("ESORIC", pred.nms))>0) {
    wv.es <- rep(1, length(grep("ESORIC_", mcmp[["selected.mdls"]]$sel.cri)))
  }

  ## variable contributions and importance
  var.cont.df <- matrix(nrow = length(mxnt.mdls), ncol = length(var.nms))
  rownames(var.cont.df) <- mod.nms
  colnames(var.cont.df) <- var.nms
  var.permImp.df <- var.cont.df

  for(i in 1:nrow(var.cont.df)){
    var.cont.df[i,] <- mxnt.mdls[[i]]@results[grepl("contribution", rownames(mxnt.mdls[[i]]@results))]
    var.permImp.df[i,] <- mxnt.mdls[[i]]@results[grepl("permutation.importance", rownames(mxnt.mdls[[i]]@results))]
  }

  f.wm <- function(pattern="AIC_", pred.nms, sel.mod.nms, var.nms, wv, df, dimnames1="Mod.ensemble" ){
    matrix(apply(data.frame(matrix(df[grep(pattern, sel.mod.nms),],
                                   nrow = sum(grepl(pattern, sel.mod.nms)), byrow = FALSE ) ), 2, function(x, wv) {
                                     stats::weighted.mean(x, wv)
                                   }, wv), nrow = 1, dimnames = list(dimnames1, var.nms) )
  }

  var.cont.df <- as.data.frame(rbind(
    if(sum(grepl("AvgAIC", pred.nms))>0){
      f.wm("AIC_", pred.nms, sel.mod.nms, var.nms, wv.aic, var.cont.df, dimnames1="Mod.AvgAIC")
    },
    if(sum(grepl("WAAUC", pred.nms))>0){
      f.wm("WAAUC_", pred.nms, sel.mod.nms, var.nms, wv.wa, var.cont.df, dimnames1="Mod.WAAUC")
    },
    if(sum(grepl("EBPM", pred.nms))>0){
      f.wm("EBPM_", pred.nms, sel.mod.nms, var.nms, wv.bp, var.cont.df, dimnames1="Mod.EBPM")
    },
    if(sum(grepl("ESORIC", pred.nms))>0){
      f.wm("ESORIC_", pred.nms, sel.mod.nms, var.nms, wv.es, var.cont.df, dimnames1="Mod.ESORIC")
    },
    var.cont.df))

  var.permImp.df <- as.data.frame(rbind(
    if(sum(grepl("AvgAIC", pred.nms))>0){
      f.wm("AIC_", pred.nms, sel.mod.nms, var.nms, wv.aic, var.permImp.df, dimnames1="Mod.AvgAIC")
    },
    if(sum(grepl("WAAUC", pred.nms))>0){
      f.wm("WAAUC_", pred.nms, sel.mod.nms, var.nms, wv.wa, var.permImp.df, dimnames1="Mod.WAAUC")
    },
    if(sum(grepl("EBPM", pred.nms))>0){
      f.wm("EBPM_", pred.nms, sel.mod.nms, var.nms, wv.bp, var.permImp.df, dimnames1="Mod.EBPM")
    },
    if(sum(grepl("ESORIC", pred.nms))>0){
      f.wm("ESORIC_", pred.nms, sel.mod.nms, var.nms, wv.es, var.permImp.df, dimnames1="Mod.ESORIC")
    },
    var.permImp.df))

  mnms.i <- is.na(match(rownames(var.cont.df), mod.nms))
  sel.mod.nms <- c(rownames(var.cont.df)[mnms.i], sel.mod.nms)
  var.cont.df <- cbind(sel.crit=sel.mod.nms, var.cont.df)
  # var.cont.df$sel.crit <- as.character(var.cont.df$sel.crit)
  # var.cont.df$sel.crit[!is.na(match(rownames(var.cont.df), mod.nms))] <- sel.mod.nms
  var.permImp.df <- cbind(sel.crit=sel.mod.nms, var.permImp.df)

  # var.contPermImp[[sp]] <- array(c(as.matrix(var.cont.df), as.matrix(var.permImp.df)), c(nrow(var.cont.df), ncol(var.cont.df), 2), dimnames = c(dimnames(var.cont.df), list(c("contribution", "permutation.importance") )))
  # utils::write.csv(var.cont.df, paste0("3_out.MaxEnt/Mdls.", sp.nm, "/metric.var.Contribution.", sp.nm, ".csv"))
  # utils::write.csv(var.permImp.df, paste0("3_out.MaxEnt/Mdls.", sp.nm, "/metric.var.PermImportance", sp.nm, ".csv"))
  # var.contPermImp[[sp]] <- list(contribution=var.cont.df, permutation.importance=var.permImp.df)
  return(list(contribution=var.cont.df, permutation.importance=var.permImp.df))
}


#' Compute "Fractional predicted area" ('n of occupied pixels'/n)
#'
#' Compute "Fractional predicted area" ('n of occupied pixels'/total n) or ('area of occupied pixels'/total area)
#'
#' @inheritParams get_tsa_b
#' @seealso \code{\link{get_fpa}}, \code{\link{get_tsa}}, \code{\link{get_cont_permimport}}
#' @seealso \code{\link{get_tsa_b}}, \code{\link{get_cont_permimport_b}}
#' @return A list of species' FPAs computed for each climatic scenario, threshold and model criteria
#' @examples
#' \dontrun{
#' get_fpa_b(mtp.l=mods.thrshld.lst)
#' }
#' @export
get_fpa_b <- function(mtp.l, digits = 3){
  # df.FPA <- vector("list", length = length(mtp.l))

  df.FPA <- lapply(seq_along(mtp.l), function(i, mtp.l, digits){
    get_fpa(mtp.l[[i]], digits)
    }, mtp.l, digits) # species, areas

  names(df.FPA) <- names(mtp.l)
  df.FPA.c <- data.table::rbindlist(df.FPA, idcol = "sp")
  utils::write.csv(df.FPA.c, paste0("3_out.MaxEnt/metric.FracPredArea.csv")) # reorder ds

  return(df.FPA)
}

#' Compute "Fractional predicted area" ('n of occupied pixels'/n)
#'
#' Compute "Fractional predicted area" ('n of occupied pixels'/total n) or ('area of occupied pixels'/total area)
#'
#' @inheritParams get_tsa
#' @seealso \code{\link{get_fpa_b}}, \code{\link{get_tsa}}, \code{\link{get_cont_permimport}}
#' @return A list of species' FPAs computed for each climatic scenario, threshold and model criteria
#' @examples
#' \dontrun{
#' get_fpa(mtp.l=mods.thrshld.lst)
#' }
#' @export
get_fpa <- function(mtp, digits){ # species, areas
  # print(sp.nm)
  # areas <- array(dim=c(length(mtp), # rows for pred.scenario
  #                      length(mtp[[1]][[2]]), # cols for threshold criteria
  #                      raster::nlayers(mtp[[1]][[2]][[1]])), # sheet (3rd dim) for model criteria
  #                dimnames = list(names(mtp), # pred.scenario
  #                                names(mtp[[1]][[2]]), # threshold criteria
  #                                gsub(paste(c(".mxnt.pred.", ".current.", "Mod.", "fcv1", "fcv5",
  #                                             "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto",
  #                                             "eetd", paste0(".", names(mtp), ".") ), collapse = "|"), "", names(mtp[[1]][[2]][[1]]))
  #                )) # model criteria
  # areas <- data.table::melt(areas)
  # colnames(areas)[1:4] <- c("Clim.scen", "threshold", "Model", "FPA")
  #
  # # areas <- areas
  # # mtp.l.sp <- mtp
  # areas <- expand.grid(Clim.scen=names(mtp),
  #                      threshold=names(mtp[[1]][[2]]),
  #                      Model=gsub(paste(c(".mxnt.pred.", ".current.", "Mod.", "fcv1", "fcv5",
  #                                         "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto",
  #                                         "eetd", paste0(".", names(mtp), ".") ), collapse = "|"), "", names(mtp[[1]][[2]][[1]])),
  #                      FPA=NA)

  fpa.mods.t.p <- lapply(seq_along(mtp), function(sc, mtp, digits){  # pred.scenario
    mtp.sc <- mtp[[sc]][[2]]

    fpa.mods.t <- sapply(seq_along(mtp.sc), function(t, mtp.sc, sc, digits){ # threshold criteria
      mtp.sc.t <- mtp.sc[[t]]

      fpa.mods <- sapply(1:raster::nlayers(mtp.sc.t), function(m, mtp.sc.t, sc, t, digits){ # model criteria
        ar <- mtp.sc.t[[m]]

        FPA <- (sum(raster::area(ar, na.rm=TRUE)[raster::getValues(ar)==1], na.rm=TRUE)/
                  sum(raster::area(ar, na.rm=TRUE)[!is.na(raster::getValues(ar))], na.rm=TRUE) )

        return(FPA) }, mtp.sc.t, sc, t, digits) # model criteria
      return(fpa.mods) }, mtp.sc, sc, digits) # threshold criteria
    return(fpa.mods.t) }, mtp, digits) # pred.scenario

  fpa.mods.t.p <- simplify2array(fpa.mods.t.p)
  if(length(dim(fpa.mods.t.p))==3){
    fpa.mods.t.p <- round(array(aperm(fpa.mods.t.p, c(3,2,1))), digits = digits) #,
  } else if(length(dim(fpa.mods.t.p))==2){
    dim(fpa.mods.t.p) <- c(dim(fpa.mods.t.p), 1)
    fpa.mods.t.p <- round(array(aperm(fpa.mods.t.p, c(3,2,1))), digits = digits) #,
  } else if(length(dim(fpa.mods.t.p))==1){
    dim(fpa.mods.t.p) <- c(dim(fpa.mods.t.p), 1, 1)
    fpa.mods.t.p <- round(array(aperm(fpa.mods.t.p, c(3,2,1))), digits = digits) #,
  } #else if(is.null(dim(fpa.mods.t.p))){
  #   fpa.mods.t.p <- fpa.mods.t.p
  # }
  areas <- data.frame(expand.grid(Clim.scen=names(mtp),
              threshold=names(mtp[[1]][[2]]),
              Model=gsub(paste(c(".mxnt.pred.", ".current.", "Mod.", "fcv1", "fcv5",
                                 "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto",
                                 "eetd", paste0(".", names(mtp), ".") ), collapse = "|"), "", names(mtp[[1]][[2]][[1]]))),
              FPA=fpa.mods.t.p)

  # utils::write.csv(areas, paste0("3_out.MaxEnt/Mdls.", sp.nm, "/metric.FracPredArea.", sp.nm, ".csv")) # reorder ds
  return(areas)
}




#' Compute "Omission Rate"
#'
#' Compute "Omission Rate" of species occurence points for a climatic scenario (usually "current")
#'
#' @inheritParams get_tsa_b
#' @inheritParams ENMevaluate_b
#' @param clim.scn.nm name to locate climatic scenario from which Omission Rate will
#' be extracted. Usually the scenario used to calibrate maxent models
#' @seealso \code{\link{get_tsa}}, \code{\link{get_cont_permimport}}, \code{\link{get_fpa}}
#' @return A list of species' ORs computed for the selected (current) climatic scenario and
#' each threshold and model criteria
#' @examples
#' \dontrun{
#' get_OR(mtp.l=mods.thrshld.lst, occ.l=occ.locs)
#' }
# #'@export
get_OR <- function(mtp.l, occ.l, clim.scn.nm = "current", digits = 3){ # , save=TRUE
  if(is.null(clim.scn.nm)){
    stop("Need to specify 'clim.scn.nm'")
  }
  df.OmR <- vector("list")
  for(sp in names(mtp.l)){ # species
    occ.spdf <- occ.l[[sp]]
    if(!class(occ.spdf) %in% c("SpatialPoints", "SpatialPointsDataFrame")){
      lon.col <- colnames(occ.spdf)[grep("^lon$|^long$|^longitude$", colnames(occ.spdf), ignore.case = T, fixed = F)][1]
      lat.col <- colnames(occ.spdf)[grep("^lat$|^latitude$", colnames(occ.spdf), ignore.case = T)][1]
      sp::coordinates(occ.spdf) <- c(lon.col, lat.col)
    }
    N.pts <- length(occ.spdf)
    ci <- grep(clim.scn.nm, names(mtp.l[[sp]]))
    if(length(ci)<1){
      stop("No climatic scenario named as: ", clim.scn.nm)
    }
    trlds <- names(mtp.l[[sp]][[ci]]$binary)
    thrshld.nms <- paste0(".", trlds, collapse = "|") # c("fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd")
    mdls <- gsub(paste(c(thrshld.nms, "Mod.", ".current"), collapse = "|"), "", names(mtp.l[[sp]][[ci]]$binary[[1]]))
    nr <- length(mdls)
    nc <- length(trlds)
    df.OmR[[sp]] <- data.frame(matrix(nrow=nr, ncol=nc), Model=mdls)
    colnames(df.OmR[[sp]])[1:length(trlds)] <- trlds

    for(t in names(mtp.l[[sp]][[ci]]$binary)){ # threshold criteria
      for(m in 1:raster::nlayers(mtp.l[[sp]][[ci]]$binary[[t]])){ # model criteria
        df.OmR[[sp]][m, t] <- round((1-(sum(raster::extract(mtp.l[[sp]][[ci]]$binary[[t]][[m]], occ.spdf), na.rm = T)/N.pts) ), digits)
      } # model criteria
    } # threshold criteria

    df.OmR[[sp]] <- data.table::melt(data.table::data.table(df.OmR[[sp]]), id.vars="Model", variable.name="threshold", value.name="OmR") # reshape2::melt(df.OmR[[sp]], id="Model") #
    # colnames(df.OmR[[sp]])[1:3] <- c("Model", "threshold", "OmR")
    # utils::write.csv(as.data.frame(df.OmR[[sp]]), paste0("3_out.MaxEnt/Mdls.", sp, "/metric.OmRate", sp, ".csv")) # reorder ds
  }
  df.OmR.c <- data.table::rbindlist(df.OmR, idcol = "sp")
  utils::write.csv(df.OmR.c, paste0("3_out.MaxEnt/metric.OmRate.csv")) # reorder ds
  return(OmR = df.OmR)
}
