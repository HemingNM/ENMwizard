# # ##### 5. METRICS
# # #### 4.6 compute area occupied
# #' Compute total suitable area
# #'
# #' General function description. A short paragraph (or more) describing what the function does.
# #' @inheritParams f.plot.mxnt.preds
# #' @param pred.nm name of prediction to be appended to the final name. Usually "pres", "past" or "fut".
# #' @param thrshld.i List of threshold criteria to be applied
# #' @return Stack or brick of thresholded predictions
# #' @examples
# #' areas.occ.lst <- f.area.occ(mtp.spl)
# #' @export
# f.area.occ <- function(mtp.spl){
#   area.occ.spp <- vector("list", length = length(mtp.spl))
#   names(area.occ.spp) <- names(mtp.spl)
#   areas <- matrix(ncol=length(mtp.spl[[1]][[2]]), nrow = nlayers(mtp.spl[[1]][[2]][[1]]))
#   rownames(areas) <- gsub(".mtp", "", names(mtp.spl[[1]][[2]][[1]]))
#   colnames(areas) <- names(mtp.spl[[1]][[2]])
#   for(sp in 1:length(mtp.spl)){
#     print(paste(names(mtp.spl)[sp]))
#     area.occ.spp[[sp]] <- areas
#     for(t in 1:length(mtp.spl[[sp]][[2]])){
#       # print(paste(names(mtp.spl[[sp]][[2]])[t]))
#       for(m in 1:nlayers(mtp.spl[[sp]][[2]][[t]])){
#         print(paste(names(mtp.spl[[sp]][[2]][[t]])[m]))
#         ar <- sum(area(mtp.spl[[sp]][[2]][[t]][[m]], na.rm=TRUE)[getValues(mtp.spl[[sp]][[2]][[t]][[m]])==1], na.rm=TRUE)
#         print(ar)
#         area.occ.spp[[sp]][m,t] <- ar
#       }
#     }
#   }
#
#   return(area.occ.spp)
# }
#
# areas.occ.lst <- f.area.occ(mtp.spl)


# #### 5.1 compute area occupied for multiple scenarios

#' Compute species' total suitable area
#'
#' Compute total suitable area for multiple climatic scenario, threshold and model criteria
#' @inheritParams f.plot.mxnt.preds
#' @param digits integer indicating the number of decimal places. see ?round for details.
#' @param restrict a raster to select a region to compute area.
#' @return List of arrays containing species' total suitable areas for each climatic scenario, threshold and model criteria
#' @examples
#' areas.occ.lst <- f.area.occ.mscn(mods.thrshld.lst)
#' @export
f.area.occ.mscn <- function(mtp.spl, restrict=NULL, digits=0){
  area.occ.spp <- vector("list", length = length(mtp.spl))
  names(area.occ.spp) <- names(mtp.spl)
  thrshld.nms <- c("fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd")
  thrshld.nms <- paste(paste0(".",thrshld.nms), collapse = "|")
  # c.nms <- gsub(".mxnt.preds","",gsub(thrshld.nms,"",names(mtp.spl[[1]][[1]][[2]][[1]]) ))
  c.nms <- names(mtp.spl[[1]][[1]][[2]][[1]])
  m.nms <- c("LowAICc", "Mean.ORmin", "Mean.OR10", "Mean.AUCmin", "Mean.AUC10", "test", "AvgAICc")
  invisible(sapply(seq_along(m.nms), function(i, x, y){
    if(sum(grepl(m.nms[i], c.nms))>0){
      c.nms[grepl(m.nms[i], c.nms)] <<- m.nms[i]
    }
  }, c.nms, m.nms))

  areas <- array(dim=c(length(mtp.spl[[1]]), # rows for pred.scenario
                       length(mtp.spl[[1]][[1]][[2]]), # cols for threshold criteria
                       raster::nlayers(mtp.spl[[1]][[1]][[2]][[1]])), # sheet (3rd dim) for model criteria
                 dimnames = list(names(mtp.spl[[1]]), # pred.scenario
                                 names(mtp.spl[[1]][[1]][[2]]), # threshold criteria
                                 c.nms )) # model criteria
  #
  thrshld.crit <- names(mtp.spl[[1]][[1]][[1]])
  areas.occ.df <- vector("list")

  area.occ.spp <- lapply(names(mtp.spl), function(sp, mtp.spl, areas, restrict, digits){ # species
    # for(sp in names(mtp.spl)){ # species
    print(sp)
    area.occ.spp[[sp]] <- areas

    mtp.spl.sp <- mtp.spl[[sp]]

    # area.occ.spp[[sp]] <-
    ar.mods.t.p <- lapply(seq_along(mtp.spl.sp), function(sc, mtp.spl.sp, sp, restrict, digits, area.occ.spp){  # pred.scenario
      # for(sc in seq_along(mtp.spl.sp)){  # pred.scenario

      mtp.spl.sp.sc <- mtp.spl.sp[[sc]][[2]]

      ar.mods.t <- sapply(seq_along(mtp.spl.sp.sc), function(t, mtp.spl.sp.sc, sp,sc, restrict, digits, area.occ.spp){ # threshold criteria
        # for(t in seq_along(mtp.spl.sp.sc)){ # threshold criteria

        mtp.spl.sp.sc.t <- mtp.spl.sp.sc[[t]]

        ar.mods <- sapply(1:raster::nlayers(mtp.spl.sp.sc.t), function(m, mtp.spl.sp.sc.t, sp,sc,t, restrict, digits, area.occ.spp){ # model criteria
          # for(m in 1:raster::nlayers(mtp.spl.sp.sc.t)){ # model criteria

          ar <- mtp.spl.sp.sc.t[[m]]

          if(grDevices::is.raster(restrict)){
            if(raster::res(ar)!=raster::res(restrict)){
              ar <- raster::resample(ar, restrict)
              ar <- ar*restrict
            }
          }
          ar <- sum(raster::area(ar, na.rm=TRUE)[raster::getValues(ar)==1], na.rm=TRUE)
          ar <- round(ar, digits = digits)
          print(paste(names(mtp.spl.sp.sc.t)[m], " area is ", ar))
          area.occ.spp[[sp]][sc,t,m] <<- ar
          # return(ar)
          # } # model criteria
          return(ar) }, mtp.spl.sp.sc.t, sp,sc,t, restrict, digits, area.occ.spp) # model criteria

        # } # threshold criteria
        return(ar.mods) }, mtp.spl.sp.sc, sp,sc, restrict, digits, area.occ.spp) # threshold criteria

      # return(area.occ.spp[[sp]])
      # } # pred.scenario
      return(ar.mods.t) }, mtp.spl.sp, sp, restrict, digits, area.occ.spp) # pred.scenario

    area.occ.spp[[sp]][] <- array(aperm(simplify2array(ar.mods.t.p), c(3,2,1))) #,
    # dim = dim(areas),
    # dimnames = list(names(mtp.spl[[1]]), # pred.scenario
    #                 names(mtp.spl[[1]][[1]][[2]]), # threshold criteria
    #                 c.nms )) # model criteria


    #### TO DO export areas
    # c.nms <- names(mtp.spl[[sp]][[1]][[2]][[1]])
    # if(length(c.nms)>1){
    #   c.nms <- c("Total", gsub(c.nms[1], "", c.nms[2:length(c.nms)]))
    # } else {c.nms <- "Total"}
    # sp.nm <- names(mtp.spl)[sp]
    areas.occ.df[[sp]] <- as.data.frame(area.occ.spp[[sp]]) #
    colnames(areas.occ.df[[sp]]) <- paste(thrshld.crit, rep(c.nms, each=length(thrshld.crit)), sep = ".")
    xlsx::write.xlsx(areas.occ.df[[sp]], paste0("4_ENMeval.results/Mdls.", sp, "/areaT.", sp, ".xlsx"))

    # } # species
    return(area.occ.spp[[sp]]) }, mtp.spl, areas, restrict, digits) # species

  names(area.occ.spp) <- names(mtp.spl)
  return(area.occ.spp)
}



# #### 4.7 extract model results
# ### 4.7.1 variable contribution and importance
#' Compute variable contribution and importance
#'
#' General function description. A short paragraph (or more) describing what the function does.
# #' @param mmp.spl Stack or brick of predictions to apply the threshold
#' @inheritParams f.thr.batch
#' @return List of arrays containing variable contribution and importance for each species
#' @examples
#' f.var.ci(mxnt.mdls.preds.lst)
#' @export
f.var.ci <- function(mmp.spl){
  path.res <- "4_ENMeval.results"
  if(dir.exists(path.res)==F) dir.create(path.res)
  path.sp.m <- paste0("Mdls.", names(mmp.spl))
  path.mdls <- paste(path.res, path.sp.m, sep="/")

  var.contPermImp <- stats::setNames(vector("list", length(mmp.spl)), names(mmp.spl))
  for(sp in seq_along(mmp.spl)){
    mxnt.mdls <- mmp.spl[[sp]]$mxnt.mdls
    mod.nms <- mmp.spl[[sp]]$selected.mdls$sel.cri
    var.nms <- gsub( ".contribution", "", rownames(mxnt.mdls[[1]]@results)[grepl("contribution", rownames(mxnt.mdls[[1]]@results))])
    w.mdls <- mmp.spl[[sp]]$selected.mdls$w.AIC

    ## variable contributions and importance
    var.cont.df <- matrix(nrow = length(mxnt.mdls), ncol = length(var.nms))
    rownames(var.cont.df) <- mod.nms
    colnames(var.cont.df) <- var.nms
    var.permImp.df <- var.cont.df

    for(i in 1:nrow(var.cont.df)){
      var.cont.df[i,] <- mxnt.mdls[[i]]@results[grepl("contribution", rownames(mxnt.mdls[[i]]@results))]
      var.permImp.df[i,] <- mxnt.mdls[[i]]@results[grepl("permutation.importance", rownames(mxnt.mdls[[i]]@results))]
    }

    var.cont.df <- rbind(Mod.Avg.AICc = apply(var.cont.df[grep("Mod.AICc",mod.nms),], 2,
                                              function(x) sum(x*w.mdls[grep("Mod.AICc",mod.nms)])),
                         var.cont.df[grep("Mod.Mean.ORmin|Mod.Mean.OR10|Mod.Mean.AUCmin|Mod.Mean.AUC10",mod.nms),],
                         var.cont.df[grep("Mod.AICc",mod.nms),])

    var.permImp.df <- rbind(Mod.Avg.AICc = apply(var.permImp.df[grep("Mod.AICc",mod.nms),], 2,
                                                 function(x) sum(x*w.mdls[grep("Mod.AICc",mod.nms)])),
                            var.permImp.df[grep("Mod.Mean.ORmin|Mod.Mean.OR10|Mod.Mean.AUCmin|Mod.Mean.AUC10",mod.nms),],
                            var.permImp.df[grep("Mod.AICc",mod.nms),])
    # save into list subitem
    # var.cont.df <- var.cont.df[c(1,2,(nrow(var.cont.df)-3):nrow(var.cont.df)),]
    # var.permImp.df <- var.permImp.df[c(1,2,(nrow(var.permImp.df)-3):nrow(var.permImp.df)),]
    var.contPermImp[[sp]] <- array(c(var.cont.df,var.permImp.df), c(nrow(var.cont.df), ncol(var.cont.df), 2), dimnames = c(dimnames(var.cont.df), list(c("contribution", "permutation.importance") )))
    xlsx::write.xlsx(var.contPermImp[[sp]][,,1], paste0(path.mdls[sp],"/var.contPermImp.", names((mmp.spl)[sp]), ".xlsx"), sheetName="contribution")
    xlsx::write.xlsx(var.contPermImp[[sp]][,,2], paste0(path.mdls[sp],"/var.contPermImp.", names((mmp.spl)[sp]), ".xlsx"), append=T, sheetName="permutation.importance")
  }
  # var.cont.df
  # var.permImp.df

  return(var.contPermImp)
}









# #### 5.2 calcular "Fractional predicted area" (n de pixels ocupados/n)
#' Compute "Fractional predicted area" ('n of occupied pixels'/n)
#'
#' Compute "Fractional predicted area" ('n of occupied pixels'/total n) or  ('area of occupied pixels'/total area)
#' @inheritParams f.area.occ.mscn
#' @param occ.l list of species occurrence data.
#' @param current.pred.nm name to locate climatic scenario (usually "current") used to calibrate maxent models
#' @return A list of species' ORs computed for the selected (current) climatic scenario and
#' each threshold and model criteria
#' @examples
#' f.OR(mods.thrshld.lst, occ.locs, "current")
#' @export
f.OR <- function(mtp.spl, occ.l, current.pred.nm = "current", digits = 3){ # , save=T
  # library(data.table)
  df.OmR <- vector("list")
  # df.FPA <- df.OmR
  for(sp in names(mtp.spl)){ # species
    extr.OBS <- as.data.frame(occ.l[[sp]])
    N.pts <- nrow(extr.OBS)
    sp::coordinates(extr.OBS) <- ~LONG+LAT
    ci <- grep(current.pred.nm, names(mtp.spl[[sp]]))
    trlds <- names(mtp.spl[[sp]][[ci]]$binary)
    thrshld.nms <- paste0(".", trlds, collapse = "|") # c("fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd")
    mdls <- gsub(paste(c(thrshld.nms, "Mod."), collapse = "|"), "", names(mtp.spl[[1]][[ci]]$binary[[1]]))
    nr <- length(mdls)
    nc <- length(trlds)
    df.OmR[[sp]] <- as.data.frame(matrix(nrow=nr, ncol=nc))
    rownames(df.OmR[[sp]]) <- mdls
    colnames(df.OmR[[sp]]) <- trlds
    # df.FPA[[sp]] <- df.OmR[[sp]]

    for(t in names(mtp.spl[[sp]][[ci]]$binary)){ # threshold criteria
      for(m in 1:raster::nlayers(mtp.spl[[sp]][[ci]]$binary[[t]])){ # model criteria
        df.OmR[[sp]][m, t] <- (1-(sum(raster::extract(mtp.spl[[sp]][[ci]]$binary[[t]][[m]], extr.OBS), na.rm = T)/N.pts) )
      } # model criteria
    } # threshold criteria

    xlsx::write.xlsx(round(df.OmR[[sp]],digits), paste0("4_ENMeval.results/Mdls.", sp, "/OmRate", sp, ".xlsx")) # reorder ds
    # write.xlsx(round(df.FPA[[sp]],digits), paste0("4_ENMeval.results/Mdls.", sp, "/FracPredArea", sp, ".xlsx")) # reorder ds
  }
  # if(save){
    df.OmR.c <- data.table::rbindlist(lapply(lapply(df.OmR, round, digits=digits), data.table::setDT, keep.rownames = TRUE), idcol = TRUE)
    # df.FPA.c <- data.table::rbindlist(lapply(lapply(df.FPA, round, digits=digits), data.table::setDT, keep.rownames = TRUE), idcol = TRUE)
    colnames(df.OmR.c)[1:2] <- c("sp", "Model")
    # colnames(df.FPA.c)[1:2] <- c("sp", "Model")

    xlsx::write.xlsx(df.OmR.c, paste0("4_ENMeval.results/OmRate.xlsx")) # reorder ds
    # xlsx::write.xlsx(df.FPA.c, paste0("4_ENMeval.results/FracPredArea.xlsx")) # reorder ds
  # }
  return(OmR = df.OmR)
}



#' Compute "Fractional predicted area" ('n of occupied pixels'/n) for multiple scenarios
#'
#' General function description. A short paragraph (or more) describing what the function does.
#' @inheritParams f.OR
#' @return A list of species' FPAs computed for each climatic scenario, threshold and model criteria
#' @examples
#' f.FPA.mscn(mtp.spl)
#' @export
f.FPA <- function(mtp.spl, digits = 3){
  df.FPA <- vector("list", length = length(mtp.spl))
  names(df.FPA) <- names(mtp.spl)

  areas <- array(dim=c(length(mtp.spl[[1]]), # rows for pred.scenario
                       length(mtp.spl[[1]][[1]][[2]]), # cols for threshold criteria
                       raster::nlayers(mtp.spl[[1]][[1]][[2]][[1]])), # sheet (3rd dim) for model criteria
                 dimnames = list(names(mtp.spl[[1]]), # pred.scenario
                                 names(mtp.spl[[1]][[1]][[2]]), # threshold criteria
                                 gsub(paste(c(".mxnt.pred.", "fcv1", "fcv5", "fcv10", "mtp", "x10ptp", "etss", "mtss", "bto", "eetd"), collapse = "|"), "", names(mtp.spl[[1]][[1]][[2]][[1]]))
                 )) # model criteria
  #


  df.FPA <- lapply(names(mtp.spl), function(sp, mtp.spl, areas, digits){ # species
    # for(sp in names(mtp.spl)){ # species
    print(sp)
    df.FPA[[sp]] <- areas

    mtp.spl.sp <- mtp.spl[[sp]]

    # df.FPA[[sp]] <-
    fpa.mods.t.p <- lapply(seq_along(mtp.spl.sp), function(sc, mtp.spl.sp, sp, digits, df.FPA){  # pred.scenario
      # for(sc in seq_along(mtp.spl.sp)){  # pred.scenario

      mtp.spl.sp.sc <- mtp.spl.sp[[sc]][[2]]

      fpa.mods.t <- sapply(seq_along(mtp.spl.sp.sc), function(t, mtp.spl.sp.sc, sp,sc, digits, df.FPA){ # threshold criteria
        # for(t in seq_along(mtp.spl.sp.sc)){ # threshold criteria

        mtp.spl.sp.sc.t <- mtp.spl.sp.sc[[t]]

        fpa.mods <- sapply(1:raster::nlayers(mtp.spl.sp.sc.t), function(m, mtp.spl.sp.sc.t, sp,sc,t, digits, df.FPA){ # model criteria
          # for(m in 1:raster::nlayers(mtp.spl.sp.sc.t)){ # model criteria
          ar <- mtp.spl.sp.sc.t[[m]]

          FPA <- (sum(raster::area(ar, na.rm=TRUE)[raster::getValues(ar)==1], na.rm=TRUE)/
                    sum(raster::area(ar, na.rm=TRUE)[!is.na(raster::getValues(ar))], na.rm=TRUE) )

          # FPA <- (length(ar[ar==1]) /
          #           length(ar[!is.na(ar)]) )

          print(paste(names(mtp.spl[[sp]][[sc]][[2]][[t]])[m], " FPA is ", FPA))

          # return(FPA)
          # } # model criteria
          return(FPA) }, mtp.spl.sp.sc.t, sp,sc,t, digits, df.FPA) # model criteria

        # } # threshold criteria
        return(fpa.mods) }, mtp.spl.sp.sc, sp,sc, digits, df.FPA) # threshold criteria

      # return(df.FPA[[sp]])
      # } # pred.scenario
      return(fpa.mods.t) }, mtp.spl.sp, sp, digits, df.FPA) # pred.scenario

    df.FPA[[sp]][] <- round(array(aperm(simplify2array(fpa.mods.t.p), c(3,2,1))), digits = digits) #,
    xlsx::write.xlsx(df.FPA[[sp]], paste0("4_ENMeval.results/Mdls.", sp, "/FracPredArea.", sp, ".xlsx")) # reorder ds

    # areas.occ.df[[sp]] <- as.data.frame(df.FPA[[sp]]) #
    # colnames(areas.occ.df[[sp]]) <- paste(thrshld.crit, rep(c.nms, each=length(thrshld.crit)), sep = ".")
    # xlsx::write.xlsx(areas.occ.df[[sp]], paste0("4_ENMeval.results/Mdls.", sp, "/areaT.", sp, ".xlsx"))

    # } # species
    return(df.FPA[[sp]]) }, mtp.spl, areas, digits) # species

  names(df.FPA) <- names(mtp.spl)
  return(df.FPA)
}


