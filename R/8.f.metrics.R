# # ##### 5. METRICS

# TO DO - f.area.occ.mscn;
# In area.occ.spp[[sp]][] <- array(aperm(ar.mods.t.p, c(3, 2, 1))) :
#   number of items to replace is not a multiple of replacement length
# TO DO - f.FPA
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
# f.area.occ <- function(mtp.l){
#   area.occ.spp <- vector("list", length = length(mtp.l))
#   names(area.occ.spp) <- names(mtp.l)
#   areas <- matrix(ncol=length(mtp.l[[1]][[2]]), nrow = nlayers(mtp.l[[1]][[2]][[1]]))
#   rownames(areas) <- gsub(".mtp", "", names(mtp.l[[1]][[2]][[1]]))
#   colnames(areas) <- names(mtp.l[[1]][[2]])
#   for(sp in 1:length(mtp.l)){
#     print(paste(names(mtp.l)[sp]))
#     area.occ.spp[[sp]] <- areas
#     for(t in 1:length(mtp.l[[sp]][[2]])){
#       # print(paste(names(mtp.l[[sp]][[2]])[t]))
#       for(m in 1:nlayers(mtp.l[[sp]][[2]][[t]])){
#         print(paste(names(mtp.l[[sp]][[2]][[t]])[m]))
#         ar <- sum(area(mtp.l[[sp]][[2]][[t]][[m]], na.rm=TRUE)[getValues(mtp.l[[sp]][[2]][[t]][[m]])==1], na.rm=TRUE)
#         print(ar)
#         area.occ.spp[[sp]][m,t] <- ar
#       }
#     }
#   }
#
#   return(area.occ.spp)
# }
#
# areas.occ.lst <- f.area.occ(mtp.l)


# #### 5.1 compute area occupied for multiple scenarios

#' Compute species' total suitable area
#'
#' Compute total suitable area for multiple climatic scenario, threshold and model criteria.
#'
#' @inheritParams plot.mdl.diff
#' @param digits integer indicating the number of decimal places. see ?round for details.
#' @param restrict a raster to select a region to compute area.
#' @seealso \code{\link[raster]{area}}, \code{\link{f.var.ci}}, \code{\link{f.OR}}, \code{\link{f.FPA}}, \code{\link{f.raster.overlap.mscn}}
#' @return List of arrays containing species' total suitable areas for each climatic scenario, threshold and model criteria
#' @examples
#' areas.occ.lst <- f.area.occ.mscn(mtp.l=mods.thrshld.lst)
#' @export
f.area.occ.mscn <- function(mtp.l, restrict=NULL, digits=0){
  area.occ.spp <- vector("list", length = length(mtp.l))
  names(area.occ.spp) <- names(mtp.l)
  thrshld.nms <- c("fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd")
  thrshld.nms <- paste(paste0(".",thrshld.nms), collapse = "|")

  areas.occ.df <- vector("list")

  area.occ.spp <- lapply(names(mtp.l), function(sp, mtp.l, restrict, digits){ # species, areas

    c.nms <- names(mtp.l[[sp]][[1]][[2]][[1]])
    m.nms <- c("LowAICc", "ORmin", "OR10", "AUCmin", "AUC10", "AvgAICc") # , "test"
    invisible(sapply(seq_along(m.nms), function(i, x, y){
      if(sum(grepl(m.nms[i], c.nms))>0){
        c.nms[grepl(m.nms[i], c.nms)] <<- m.nms[i]
      }
    }, c.nms, m.nms))

    areas <- array(dim=c(length(mtp.l[[sp]]), # rows for pred.scenario
                         length(mtp.l[[sp]][[1]][[2]]), # cols for threshold criteria
                         raster::nlayers(mtp.l[[sp]][[1]][[2]][[1]])), # sheet (3rd dim) for model criteria
                   dimnames = list(names(mtp.l[[sp]]), # pred.scenario
                                   names(mtp.l[[sp]][[1]][[2]]), # threshold criteria
                                   c.nms )) # model criteria

    thrshld.crit <- names(mtp.l[[sp]][[1]][[1]])

    print(sp)
    area.occ.spp[[sp]] <- areas

    mtp.l.sp <- mtp.l[[sp]]

    ar.mods.t.p <- lapply(seq_along(mtp.l.sp), function(sc, mtp.l.sp, sp, restrict, digits, area.occ.spp){  # pred.scenario
      mtp.l.sp.sc <- mtp.l.sp[[sc]][[2]]

      ar.mods.t <- sapply(seq_along(mtp.l.sp.sc), function(t, mtp.l.sp.sc, sp,sc, restrict, digits, area.occ.spp){ # threshold criteria
        mtp.l.sp.sc.t <- mtp.l.sp.sc[[t]]

        ar.mods <- sapply(1:raster::nlayers(mtp.l.sp.sc.t), function(m, mtp.l.sp.sc.t, sp,sc,t, restrict, digits, area.occ.spp){ # model criteria
          ar <- mtp.l.sp.sc.t[[m]]

          if(grDevices::is.raster(restrict)){
            if(raster::res(ar)!=raster::res(restrict)){
              ar <- raster::resample(ar, restrict)
              ar <- ar*restrict
            }
          }
          ar <- sum(raster::area(ar, na.rm=TRUE)[raster::getValues(ar)==1], na.rm=TRUE)
          ar <- round(ar, digits = digits)

          area.occ.spp[[sp]][sc,t,m] <<- ar
          return(ar) }, mtp.l.sp.sc.t, sp,sc,t, restrict, digits, area.occ.spp) # model criteria
        return(ar.mods) }, mtp.l.sp.sc, sp,sc, restrict, digits, area.occ.spp) # threshold criteria
      return(ar.mods.t) }, mtp.l.sp, sp, restrict, digits, area.occ.spp) # pred.scenario

    ar.mods.t.p <- simplify2array(ar.mods.t.p)
    if(length(dim(ar.mods.t.p))==3){
      area.occ.spp[[sp]][] <- array(aperm(ar.mods.t.p, c(3,2,1))) #,
    } else if(length(dim(ar.mods.t.p))==2){
      dim(ar.mods.t.p) <- c(dim(ar.mods.t.p), 1)
      area.occ.spp[[sp]][] <- array(aperm(ar.mods.t.p, c(3,2,1))) #,
    } else if(length(dim(ar.mods.t.p))==1){
      dim(ar.mods.t.p) <- c(dim(ar.mods.t.p), 1, 1)
      area.occ.spp[[sp]][] <- array(aperm(ar.mods.t.p, c(3,2,1))) #,
    } else if(is.null(dim(ar.mods.t.p))){
      area.occ.spp[[sp]][] <- ar.mods.t.p
    }

    areas.occ.df[[sp]] <- as.data.frame(area.occ.spp[[sp]]) #
    colnames(areas.occ.df[[sp]]) <- paste(thrshld.crit, rep(c.nms, each=length(thrshld.crit)), sep = ".")
    utils::write.csv(areas.occ.df[[sp]], paste0("3_out.MaxEnt/Mdls.", sp, "/totalArea", sp, ".csv"))
    return(area.occ.spp[[sp]]) }, mtp.l, restrict, digits) # species, areas

  names(area.occ.spp) <- names(mtp.l)
  area.occ.spp.c <- data.table::rbindlist(lapply(area.occ.spp, function(x) data.table::melt(x)), idcol = "sp")
  colnames(area.occ.spp.c)[1:5] <- c("sp", "Clim.scen", "threshold", "Model", "TotSuitArea")
  utils::write.csv(area.occ.spp.c, paste0("3_out.MaxEnt/totalArea.csv")) # reorder ds

  return(area.occ.spp)
}



# #### 4.7 extract model results
# ### 4.7.1 variable contribution and importance

#' Compute variable contribution and importance
#'
#' Compute variable contribution and importance for each model
#'
# #' @param mcmp.l Stack or brick of predictions to apply the threshold
#' @inheritParams f.thr.batch
#' @seealso \code{\link[dismo]{maxent}}, \code{\link{f.area.occ.mscn}}, \code{\link{f.OR}}, \code{\link{f.FPA}}, \code{\link{f.raster.overlap.mscn}}
#' @return List of arrays containing variable contribution and importance for each species
#' @examples
#' f.var.ci(mcmp.l = mxnt.mdls.preds.lst)
#' @export
f.var.ci <- function(mcmp.l){
  path.res <- "3_out.MaxEnt"
  if(dir.exists(path.res)==FALSE) dir.create(path.res)

  var.contPermImp <- stats::setNames(vector("list", length(mcmp.l)), names(mcmp.l))
  for(sp in names(mcmp.l)){
    mxnt.mdls <- mcmp.l[[sp]]$mxnt.mdls
    mod.nms <- paste0("Mod.", mcmp.l[[sp]]$selected.mdls$sel.cri)
    pred.nms <- names(mcmp.l[[sp]]$mxnt.preds[[1]])
    var.nms <- gsub( ".contribution", "", rownames(mxnt.mdls[[1]]@results)[grepl("contribution", rownames(mxnt.mdls[[1]]@results))])
    w.mdls <- mcmp.l[[sp]]$selected.mdls$w.AIC

    ## variable contributions and importance
    var.cont.df <- matrix(nrow = length(mxnt.mdls), ncol = length(var.nms))
    rownames(var.cont.df) <- mod.nms
    colnames(var.cont.df) <- var.nms
    var.permImp.df <- var.cont.df

    for(i in 1:nrow(var.cont.df)){
      var.cont.df[i,] <- mxnt.mdls[[i]]@results[grepl("contribution", rownames(mxnt.mdls[[i]]@results))]
      var.permImp.df[i,] <- mxnt.mdls[[i]]@results[grepl("permutation.importance", rownames(mxnt.mdls[[i]]@results))]
    }

    var.cont.df <- rbind( if(sum(grepl("AvgAIC", pred.nms))>0){
      matrix(apply(data.frame(var.cont.df[grep("AICc_",mod.nms),]), 2,
                           function(x) sum(x*w.mdls[grep("AICc_", mod.nms)])), nrow = 1, dimnames = list("Mod.Avg.AICc", var.nms) )
    },
    var.cont.df)

    var.permImp.df <- rbind(if(sum(grepl("AvgAIC", pred.nms))>0){
      matrix(apply(data.frame(var.permImp.df[grep("AICc_", mod.nms),]), 2,
                           function(x) sum(x*w.mdls[grep("AICc_", mod.nms)])), nrow = 1, dimnames = list("Mod.Avg.AICc", var.nms) )
    },
    var.permImp.df)
    var.contPermImp[[sp]] <- array(c(as.matrix(var.cont.df), as.matrix(var.permImp.df)), c(nrow(var.cont.df), ncol(var.cont.df), 2), dimnames = c(dimnames(var.cont.df), list(c("contribution", "permutation.importance") )))
    utils::write.csv(var.cont.df, paste0("3_out.MaxEnt/Mdls.", sp, "/var.Contribution.", sp, ".csv"))
    utils::write.csv(var.permImp.df, paste0("3_out.MaxEnt/Mdls.", sp, "/var.PermImportance", sp, ".csv"))
  }
  return(var.contPermImp)
}



#' Compute "Omission Rate"
#'
#' Compute "Omission Rate" of species occurence points for a climatic scenario (usually "current")
#'
#' @inheritParams f.area.occ.mscn
#' @param occ.l list of species occurrence data.
#' @param clim.scn.nm name to locate climatic scenario from which Omission Rate will
#' be extracted. Usually the scenario used to calibrate maxent models
#' @seealso \code{\link{f.area.occ.mscn}}, \code{\link{f.var.ci}}, \code{\link{f.FPA}}, \code{\link{f.raster.overlap.mscn}}
#' @return A list of species' ORs computed for the selected (current) climatic scenario and
#' each threshold and model criteria
##' @examples
##' f.OR(mtp.l=mods.thrshld.lst, occ.l=occ.locs, "current")
#' @export
f.OR <- function(mtp.l, occ.l, clim.scn.nm = NULL, digits = 3){ # , save=TRUE
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

    df.OmR[[sp]] <- data.table::melt(df.OmR[[sp]], id.vars="Model") # reshape2::melt(df.OmR[[sp]], id="Model") #
    colnames(df.OmR[[sp]])[1:3] <- c("Model", "threshold", "OmR")
    utils::write.csv(df.OmR[[sp]], paste0("3_out.MaxEnt/Mdls.", sp, "/OmRate", sp, ".csv")) # reorder ds
  }
  df.OmR.c <- data.table::rbindlist(df.OmR, idcol = "sp")
  utils::write.csv(df.OmR.c, paste0("3_out.MaxEnt/OmRate.csv")) # reorder ds
  return(OmR = df.OmR)
}



#' Compute "Fractional predicted area" ('n of occupied pixels'/n) for multiple scenarios
#'
#' Compute "Fractional predicted area" ('n of occupied pixels'/total n) or ('area of occupied pixels'/total area)
#'
#' @inheritParams f.OR
#' @seealso \code{\link{f.area.occ.mscn}}, \code{\link{f.var.ci}}, \code{\link{f.OR}}, \code{\link{f.raster.overlap.mscn}}
#' @return A list of species' FPAs computed for each climatic scenario, threshold and model criteria
#' @examples
#' f.FPA(mtp.l=mods.thrshld.lst)
#' @export
f.FPA <- function(mtp.l, digits = 3){
  df.FPA <- vector("list", length = length(mtp.l))
  names(df.FPA) <- names(mtp.l)

  df.FPA <- lapply(names(mtp.l), function(sp, mtp.l, digits){ # species, areas
    print(sp)
    areas <- array(dim=c(length(mtp.l[[sp]]), # rows for pred.scenario
                         length(mtp.l[[sp]][[1]][[2]]), # cols for threshold criteria
                         raster::nlayers(mtp.l[[sp]][[1]][[2]][[1]])), # sheet (3rd dim) for model criteria
                   dimnames = list(names(mtp.l[[sp]]), # pred.scenario
                                   names(mtp.l[[sp]][[1]][[2]]), # threshold criteria
                                   gsub(paste(c(".mxnt.pred.", ".current.", "Mod.", "fcv1", "fcv5",
                                                "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto",
                                                "eetd", paste0(".", names(mtp.l[[sp]]), ".") ), collapse = "|"), "", names(mtp.l[[sp]][[1]][[2]][[1]]))
                   )) # model criteria
    areas <- data.table::melt(areas) # reshape2::melt(areas)
    colnames(areas)[1:4] <- c("Clim.scen", "threshold", "Model", "FPA")

    df.FPA[[sp]] <- areas
    mtp.l.sp <- mtp.l[[sp]]

    fpa.mods.t.p <- lapply(seq_along(mtp.l.sp), function(sc, mtp.l.sp, sp, digits, df.FPA){  # pred.scenario
      mtp.l.sp.sc <- mtp.l.sp[[sc]][[2]]

      fpa.mods.t <- sapply(seq_along(mtp.l.sp.sc), function(t, mtp.l.sp.sc, sp,sc, digits, df.FPA){ # threshold criteria
        mtp.l.sp.sc.t <- mtp.l.sp.sc[[t]]

        fpa.mods <- sapply(1:raster::nlayers(mtp.l.sp.sc.t), function(m, mtp.l.sp.sc.t, sp,sc,t, digits, df.FPA){ # model criteria
          ar <- mtp.l.sp.sc.t[[m]]

          FPA <- (sum(raster::area(ar, na.rm=TRUE)[raster::getValues(ar)==1], na.rm=TRUE)/
                    sum(raster::area(ar, na.rm=TRUE)[!is.na(raster::getValues(ar))], na.rm=TRUE) )

          return(FPA) }, mtp.l.sp.sc.t, sp,sc,t, digits, df.FPA) # model criteria
        return(fpa.mods) }, mtp.l.sp.sc, sp,sc, digits, df.FPA) # threshold criteria
      return(fpa.mods.t) }, mtp.l.sp, sp, digits, df.FPA) # pred.scenario

    fpa.mods.t.p <- simplify2array(fpa.mods.t.p)
    if(length(dim(fpa.mods.t.p))==3){
      df.FPA[[sp]][,ncol(areas)] <- round(array(aperm(fpa.mods.t.p, c(3,2,1))), digits = digits) #,
    } else if(length(dim(fpa.mods.t.p))==2){
      dim(fpa.mods.t.p) <- c(dim(fpa.mods.t.p), 1)
      df.FPA[[sp]][,ncol(areas)] <- round(array(aperm(fpa.mods.t.p, c(3,2,1))), digits = digits) #,
    } else if(length(dim(fpa.mods.t.p))==1){
      dim(fpa.mods.t.p) <- c(dim(fpa.mods.t.p), 1, 1)
      df.FPA[[sp]][,ncol(areas)] <- round(array(aperm(fpa.mods.t.p, c(3,2,1))), digits = digits) #,
    } else if(is.null(dim(fpa.mods.t.p))){
      df.FPA[[sp]][,ncol(areas)] <- fpa.mods.t.p
    }

    utils::write.csv(df.FPA[[sp]], paste0("3_out.MaxEnt/Mdls.", sp, "/FracPredArea.", sp, ".csv")) # reorder ds
    return(df.FPA[[sp]]) }, mtp.l, digits) # species, areas

  names(df.FPA) <- names(mtp.l)
  df.FPA.c <- data.table::rbindlist(df.FPA, idcol = "sp")
  utils::write.csv(df.FPA.c, paste0("3_out.MaxEnt/FracPredArea.csv")) # reorder ds
  return(df.FPA)
}


