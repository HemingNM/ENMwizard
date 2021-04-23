#### 4.9 function to remove predicted areas outside protected areas (all and integral)

#' Remove areas from projections using polygons or binary rasters
#'
#' This function will remove (mask) areas from suitability projections using polygons or
#' binary rasters. This is usefull to filter predicted areas overlapping (e.g.) protected
#' areas or species habitats (forests, open areas, etc.).
#'
#' @param sp.mask spatial object or list of spatial objects to use for masking suitability projections
#' @param append Logical. Should append or not the resulting masked rasters?
#' Append (append=T): apply each mask on the raster and save each as separate layer.
#' Cascade (append=F): apply masks sequentially on the same raster and save as a single layer.
#' @param mask.nm Character vector. Names to identify the resulting masked rasters.
#' Should have the same length of 'sp.mask'.
#' @param pred.args arguments used for MaxEnt prediction. Can be found in mcmp$pred.args.
#' See resulting object from \code{\link{calib_mdl}} or \code{\link{calib_mdl_b}} functions
#' @inheritParams plot_mdl_diff_b
#' @examples
#' \dontrun{
#' mask_thr_projs_mscn_b(mtp.l=mods.thrshld.lst, sp.mask=NewWorld)
#' }
#' @export
mask_thr_projs_mscn_b <- function(mtp.l, pred.args='cloglog', sp.mask, append=F, mask.nm="msk"){
  mods.maskd <- mtp.l
  sp.mask <- as.list(sp.mask)
  r.msk <- raster::raster(raster::extent(mtp.l[[1]][[1]][[2]][[1]]),
                  crs = raster::crs(mtp.l[[1]][[1]][[2]][[1]]),
                  resolution = raster::res(mtp.l[[1]][[1]][[2]][[1]]))
  r.msk[] <- 1
  r.msk.sh <- raster::stack()
  for(sh in seq_along(sp.mask)){
    if(any(grepl("SpatialPolygons", class(sp.mask[[sh]])))){
      r.msk.sh <- raster::addLayer(r.msk.sh, raster::mask(r.msk, sp.mask[[sh]]))
    } else if(any(grepl("Raster", class(sp.mask[[sh]])))){
      r <- r.msk*sp.mask[[sh]] ## mask raster must be binary (0, 1)
      r[r==0] <- NA
      r.msk.sh <- raster::addLayer(r.msk.sh, r)
    }

  }

  # pred.args <- filename(mtp.l[[1]][[1]][[2]][[1]])
  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

  for(sp in names(mtp.l)){ # species
    path.nm <- paste0("3_out.MaxEnt/Mdls.", sp, "/", outpt, "/Mdls.mask")
    if(dir.exists(path.nm)==F) dir.create(path.nm)
    for(sc in seq_along(mtp.l[[sp]])){ # Climatic SCN
      for(t in seq_along(mtp.l[[sp]][[sc]][[2]])){ # threshold
        if(append) {
          mods.maskd[[sp]][[sc]][[1]][[t]] <- mtp.l[[sp]][[sc]][[1]][[t]]
          mods.maskd[[sp]][[sc]][[2]][[t]] <- mtp.l[[sp]][[sc]][[2]][[t]]
          for(sh in seq_along(sp.mask)){
            fl.nm.b <- paste0(path.nm, "/", gsub("Mod_|ESOR.|AvgAIC.|LowAICc.|Mean.ORmin.|Mean.OR10.|Mean.AUC_min.|Mean.AUC_10.",
                                                 "", names(mtp.l[[sp]][[sc]][[2]][[t]])[1]), mask.nm[sh], ".b", ".grd")
            fl.nm <- paste0(path.nm, "/", gsub("Mod_|ESOR.|AvgAIC.|LowAICc.|Mean.ORmin.|Mean.OR10.|Mean.AUC_min.|Mean.AUC_10.",
                                               "", names(mtp.l[[sp]][[sc]][[2]][[t]])[1]), mask.nm[sh], ".grd")

            mskd.r <- mtp.l[[sp]][[sc]][[1]][[t]] * r.msk.sh[[sh]]
            names(mskd.r) <- paste0(names(mtp.l[[sp]][[sc]][[1]][[t]]), mask.nm[sh])
            mods.maskd[[sp]][[sc]][[1]][[t]] <- raster::addLayer(mods.maskd[[sp]][[sc]][[1]][[t]], mskd.r)

            mskd.r <- mtp.l[[sp]][[sc]][[2]][[t]] * r.msk.sh[[sh]]
            names(mskd.r) <- paste0(names(mtp.l[[sp]][[sc]][[2]][[t]]), mask.nm[sh])
            mods.maskd[[sp]][[sc]][[2]][[t]] <- raster::addLayer(mods.maskd[[sp]][[sc]][[2]][[t]], mskd.r)
          }
          mods.maskd[[sp]][[sc]][[1]][[t]] <- raster::brick(mods.maskd[[sp]][[sc]][[1]][[t]], filename=fl.nm, format="raster", overwrite=T)
          mods.maskd[[sp]][[sc]][[2]][[t]] <- raster::brick(mods.maskd[[sp]][[sc]][[2]][[t]], filename=fl.nm.b, format="raster", overwrite=T)
        } else { # cascade
          for(sh in seq_along(sp.mask)){
            fl.nm.b <- paste0(path.nm, "/", gsub("Mod_|ESOR.|AvgAIC.|LowAICc.|Mean.ORmin.|Mean.OR10.|Mean.AUC_min.|Mean.AUC_10.",
                                                 "", names(mtp.l[[sp]][[sc]][[2]][[t]])[1]), mask.nm[sh], ".b", ".grd")
            fl.nm <- paste0(path.nm, "/", gsub("Mod_|ESOR.|AvgAIC.|LowAICc.|Mean.ORmin.|Mean.OR10.|Mean.AUC_min.|Mean.AUC_10.",
                                               "", names(mtp.l[[sp]][[sc]][[2]][[t]])[1]), mask.nm[sh], ".grd")
            mods.maskd[[sp]][[sc]][[1]][[t]] <- mods.maskd[[sp]][[sc]][[1]][[t]]*r.msk.sh[[sh]]
            mods.maskd[[sp]][[sc]][[2]][[t]] <- mods.maskd[[sp]][[sc]][[2]][[t]]*r.msk.sh[[sh]]
          }
          names(mods.maskd[[sp]][[sc]][[1]][[t]]) <- paste0(names(mtp.l[[sp]][[sc]][[1]][[t]]), mask.nm[sh])
          mods.maskd[[sp]][[sc]][[1]][[t]] <- raster::brick(raster::addLayer(mtp.l[[sp]][[sc]][[1]][[t]], mods.maskd[[sp]][[sc]][[1]][[t]]), filename=fl.nm, format="raster", overwrite=T)

          names(mods.maskd[[sp]][[sc]][[2]][[t]]) <- paste0(names(mtp.l[[sp]][[sc]][[2]][[t]]), mask.nm[sh])
          mods.maskd[[sp]][[sc]][[2]][[t]] <- raster::brick(raster::addLayer(mtp.l[[sp]][[sc]][[2]][[t]], mods.maskd[[sp]][[sc]][[2]][[t]]), filename=fl.nm.b, format="raster", overwrite=T)

        } # else append

      }  # threshold
    } # Climatic SCN
  } # species
  return(mods.maskd)
}
