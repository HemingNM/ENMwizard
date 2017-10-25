# working on this. Replacing functions from 2.f.proj.prep.R

# 4. Model Fitting (Calibration) and Projection
# 4.1 save rasters onto which the model will be projected in an elements called "area_proj"
# 4.1.1 select area for projection based on the extent of occ points
# "pred.a" names changed to "pred.a"

#' Select area for projection based on the extent of occ points
#'
#' This function will create SpatialPolygons that will be used to crop/mask raster/brick objects to be used on model projections. It has several options
#' The user can crop a squared area with an extent larger than the extent of occ_poly. By default, the "extent increase"
#' is the maximum of latitudinal and longitudinal extent "max(abs(ext_proj[1]-ext_proj[2]), abs(ext_proj[3]-ext_proj[4]))".
#' The result is added to each side of the occ_poly extent. This may be changed by setting "mult", which will be multiplied
#' by the "extent increase". Latitudinal and longitudinal increase may also vary independently by setting "same=F".
#'
#' The user can also set a buffer around occ_poly to cut the raster/brick object. Buffer value is, by default,
#' calculated in the same way as "extent increase", using "max(abs(ext_proj[1]-ext_proj[2]), abs(ext_proj[3]-ext_proj[4]))".
#' But the exact value can be defined through "deg.incr" and "mult". This method takes longer to run.
#'
#' @param occ_poly list of species occurencies SpatialPolygons
#' @param sp.nm name (of species) to give to saved object
# @param o.path Output path
#' @param deg.incr used to manually set the increase in prediction area, relative to occ_poly
#' @param mult how much increase or decrease buffer
#' @param buffer should the area be cut using a buffer around occ_poly?
#' @param same should latitudinal and longitudinal increase vary independently?
#' @return  SpatialPolygons (enlarged occ_poly)
# #' @examples
#'
#' @export
pred.a.poly <- function(occ_poly, sp.nm="sp", deg.incr=NULL, mult=1, buffer=F, same=T){ #, o.path = "occ_poly"
  path.proj <- "1_sppData/area_proj_poly"
  if(dir.exists("1_sppData")==F) dir.create("1_sppData")
  # if(dir.exists("1_sppData/area_proj")==F) dir.create("1_sppData/area_proj")
  if(dir.exists(path.proj)==F) dir.create(path.proj)
  # if(dir.exists(paste(path.proj, sp.nm, sep="/"))==F) dir.create(paste(path.proj, sp.nm, sep="/"))

  # pensar melhor como projetar
  #if(is.na(projection(occ_poly))){print} # projection(occ_poly) <- crs.set
  #if(is.na(projection(env_uncut))){} # projection(env_uncut) <- crs.set

  ext_proj <- raster::extent(occ_poly)
  cat(paste("Old extent"), "\n")
  print(ext_proj)
  cat("\n")
  if(is.null(deg.incr)){
    deg.incr <- max(c(abs(ext_proj[1]-ext_proj[2]), abs(ext_proj[3]-ext_proj[4])))*mult
    #deg.incr <- mean(c(abs(ext_proj[1]-ext_proj[2]), abs(ext_proj[3]-ext_proj[4])))*mult
  }
  cat(c("Buffer width:", deg.incr, "\n"))

  if(buffer == F){
    if(same == T){# increase the area to project by max() range of lat and lon
      for(i in 1:length(as.vector(ext_proj))){
        # cat(paste("Change in",slotNames(ext_proj)[i], ":", deg.incr, "\n"))
        ifelse(i %% 2 == 0, ext_proj[i] <- ext_proj[i] + deg.incr, ext_proj[i] <- ext_proj[i] - deg.incr)
      }
    } else {
      deg.incr.x <- (abs(ext_proj[1]-ext_proj[2]))
      deg.incr.y <- (abs(ext_proj[3]-ext_proj[4]))
      # increase the area to project by lat and lon range separately
      for(i in 1:length(as.vector(ext_proj))){
        ifelse(grepl("y", methods::slotNames(ext_proj)[i]), deg.incr <- deg.incr.y*mult, deg.incr <- deg.incr.x*mult)
        cat(paste("Change in", methods::slotNames(ext_proj)[i], ":", deg.incr, "\n"))
        ifelse(i %% 2 == 0, ext_proj[i] <- ext_proj[i] + deg.incr, ext_proj[i] <- ext_proj[i] - deg.incr)
      }
    }
    cat("\n", paste("New extent"), "\n")
    print(ext_proj)
    area_p <- methods::as(ext_proj, "SpatialPolygons")
    raster::projection(area_p) <- raster::projection(occ_poly)
    # area_p <- raster::shapefile(paste(o.path, paste0(sp.nm, "_sel.area", ".shp"), sep = "/" ))
  } else {
    area_p <- rgeos::gBuffer(occ_poly, width = deg.incr*mult, quadsegs=100)
    cat("\n", paste("New extent"), "\n")
    print(raster::extent(area_p))
    raster::projection(area_p) <- raster::projection(occ_poly)
  }
  raster::shapefile(area_p, filename = paste(path.proj, paste0(sp.nm, "_a.proj", ".shp"), sep = "/" ), overwrite=TRUE)
  # area_p <- raster::shapefile(paste(path.proj, paste0(sp.nm, "_a.proj", ".shp"), sep = "/" ))
  return(area_p)
}


#' Select area for projection based on the extent of occ points
#'
#' This function is a wrapper for "pred.a". See ?pred.a
#' It works with a named list of occ_polys to delimitate the projection area for each of the species
#' @param occ_polys list of occ_poly SpatialPolygons.  see ?pred.a.poly for details
#' @inheritParams pred.a.poly
#' @return  named list of SpatialPolygons (enlarged occ_poly)
# TODO - examples
#' @examples
#'area_projection <- pred.a.polyb(occ_polys, mult=.55, buffer=F)
#'plot(area_projection[[1]][[1]])
#' @export
pred.a.polyb <- function(occ_polys, deg.incr=NULL, mult=1, buffer=F, same=T){ # , o.path = "occ_poly"
  area_pl <- vector("list", length(occ_polys))
  names(area_pl) <- names(occ_polys)
  # cat("\n", deg.incr, "\n")
  if(!is.null(deg.incr) & length(deg.incr) == 1) { deg.incr <- rep(deg.incr, length(occ_polys))}
  if(!is.null(deg.incr) & (length(deg.incr) != length(occ_polys))) {stop("length of deg.incr must be of same of occ_polys")}
  # cat("\n", deg.incr, "\n")
  for(i in 1:length(occ_polys)){
    if(is.null(deg.incr)){
      deg.incr.i <- NULL
    }else{
      deg.incr.i <- deg.incr[i]
    }
    cat(c("\n", "Creating projection area for", names(occ_polys)[i],"\n",  "Species",  i, "of", length(occ_polys),"\n"))
    area_pl[[i]] <- pred.a.poly(occ_polys[[i]], deg.incr=deg.incr.i, mult=mult, sp.nm = names(occ_polys)[i], buffer = buffer) #, o.path = o.path
  }
  return(area_pl)
}


#' Cut area for projection based on SpatialPolygons
#'
#' This function will use SpatialPolygons to crop/mask raster/brick objects to be used on model projections.
#' @param pred_poly list of SpatialPolygons (usually of based on species occ points)
#' @param env_uncut raster/brick of environmental variables to be cut
#' @param prj.nm climatic scenario name, usually "fut" or "past" is used to indicate a general timing of scenario
#' @param sp.nm name (of species) to give to saved object
#' @return  raster or brick cropped based on SpatialPolygons
# TODO - examples
# #' @examples
#'
#' @export
pred.a <- function(pred_poly, env_uncut, prj.nm="", sp.nm="sp"){
  path.proj <- "2_envData/area_proj"
  area_p <- pred_poly # TODO - change variable name below
  if(dir.exists(path.proj)==F) dir.create(path.proj)
  if(dir.exists(paste(path.proj, sp.nm, sep="/"))==F) dir.create(paste(path.proj, sp.nm, sep="/"))

  ext_proj <- raster::extent(area_p)
  if(all.equal(area_p, methods::as(ext_proj, "SpatialPolygons"))){
    area_p <- raster::crop(env_uncut, ext_proj,
                           file= paste(path.proj, sp.nm, paste0("areaProj.", sp.nm, prj.nm,".grd"), sep="/"),
                           format="raster", overwrite=T)
  } else {
    env_crp <- raster::crop(env_uncut, ext_proj)
    area_p <- raster::mask(env_crp, area_p,
                           file= paste(path.proj, sp.nm, paste0("areaProj.", sp.nm, prj.nm,".grd"), sep="/"),
                           format="raster", overwrite=T)
  }
  return(area_p)
}



### TODO -  descriptions
#' Cut area for projection based on a list of SpatialPolygons
#'
#' This function is a wrapper for "pred.a". See ?pred.a
#' It works with a named list of pred_poly to delimitate the projection area for each of the species
#' @inheritParams pred.a
#' @param pred_polys list of SpatialPolygons (usually of based on species occ points)
#' @return  named list of cropped raster or brick
#' @examples
#' area_projection <- pred.ab(pred_polys, env_uncut, mult=.55, buffer=F)
#' plot(area_projection[[1]][[1]])
#' @export
pred.ab <- function(pred_polys, env_uncut, prj.nm=""){ # pred_poly, env_uncut, prj.nm="", sp.nm="sp"
  if(prj.nm != ""){ prj.nm <- paste0(".", prj.nm)}
  area_pl <- vector("list", length(pred_polys))
  names(area_pl) <- names(pred_polys)
  for(i in 1:length(pred_polys)){
    cat(c("\n", "Creating projection area for", names(pred_polys)[i],"\n",  "Species",  i, "of", length(pred_polys),"\n"))
    area_pl[[i]] <- pred.a(pred_polys[[i]], env_uncut, prj.nm = prj.nm, sp.nm = names(pred_polys)[i])
  }
  return(area_pl)
}


#### 4.8.3 function to cut multiple environmental layers based on pred_polys
### TODO - update description
#' Cut multiple projection areas (climatic scenarios) for multiple species (list of SpatialPolygons)
#'
#' This function is a wrapper for "pred.a". See ?pred.a. This function delimitates the projection area for
#'  each of the species contained in the pred_polys named list and crops
#'  multiple rasters/bricks (i.e. representing distinct climatic scenaries) based on the same criteria for each species
#' @inheritParams pred.ab
#' @param env_uncut.l list of raster/brick of environmental variables to be cut
#' @param cores specify number of cores if aim run in parallel
# #' @param ext_proj
#' @return  named list of cropped list of raster/brick of environmental variables
# TODO - examples
# #' @examples
#'
#' @export
pred.ab.mscn <- function(pred_polys, env_uncut.l, prj.nm="", cores=1){ # , ext_proj=NULL
  # if(prj.nm != ""){ prj.nm <- paste0(".", prj.nm)}
  path.proj <- "2_envData/area_proj"
  if(dir.exists(path.proj)==F) dir.create(path.proj)
  area_pl <- vector("list", length(pred_polys))
  names(area_pl) <- names(pred_polys)
  for(i in seq_along(pred_polys)){
    area_p.spi <- vector("list", length(env_uncut.l))
    # names(area_p.spi) <- names(env_uncut.l)
    # for(j in seq_along(area_p.spi)){
    if(cores>1){
      area_pl[[i]] <- unlist(parallel::mclapply(seq_along(area_p.spi), mc.cores = getOption("mc.cores", as.integer(cores)), function(j){
        prj.nm.j <- gsub("[..]",".",paste("", prj.nm, names(env_uncut.l)[j], sep="."))
        area_p.spi[[j]] <- pred.a(pred_polys[[i]], env_uncut.l[[j]], prj.nm = prj.nm.j, sp.nm = names(pred_polys)[i]) # ext_proj,
        area_p.spi[[j]] <- stats::setNames(area_p.spi[j], paste0(prj.nm, ".", names(env_uncut.l)[j]) )
      } ) )
    } else {
      area_pl[[i]] <- unlist(lapply(seq_along(area_p.spi), function(j){
        prj.nm.j <- gsub("[..]",".",paste("", prj.nm, names(env_uncut.l)[j], sep="."))
        area_p.spi[[j]] <- pred.a(pred_polys[[i]], env_uncut.l[[j]], prj.nm = prj.nm.j, sp.nm = names(pred_polys)[i]) # ext_proj,
        area_p.spi[[j]] <- stats::setNames(area_p.spi[j], paste0(prj.nm, ".", names(env_uncut.l)[j]) )
      } ) )
    }

    # area_pl[[i]] <- area_p.spi
  }
  return(area_pl)
}

### TODO - update description
### TODO - check this function became very similar to pred.a
#' Cut a projection area based on a SpatialPolygon (e.g. Ecoregion)
#'
#' This function will use a single SpatialPolygon to crop/mask raster/brick objects to be used on model projections.
#' @param area_p SpatialPolygon to be used as reference to crop/mask environmental variables of all species
#' @inheritParams pred.a
#' @inheritParams pred.ab
#' @param mask should use area_p to "mask" or "crop" env_uncut? See ?raster::mask and ?raster::crop for details
#' @return environmental layers (raster/brick) cutted
# TODO - examples
# #' @examples
#'
#' @export
pred.a.rst <- function(area_p, env_uncut, mask=F, prj.nm="", sp.nm="sp"){ # , crs.set
  path.proj <- "2_envData/area_proj"
  if(dir.exists(path.proj)==F) dir.create(path.proj)
  if(dir.exists(paste(path.proj, sp.nm, sep="/"))==F) dir.create(paste(path.proj, sp.nm, sep="/"))

  ext_proj <- raster::extent(area_p)
  if(mask==F){
    area_p <- raster::crop(env_uncut, ext_proj,
                   file= paste(path.proj, sp.nm, paste0("areaProj.", sp.nm, prj.nm,".grd"), sep="/"),
                   format="raster", overwrite=T)
  } else {
    env_crp <- raster::crop(env_uncut, ext_proj)
    area_p <- raster::mask(env_crp, area_p,
                   file= paste(path.proj, sp.nm, paste0("areaProj.", sp.nm, prj.nm,".grd"), sep="/"),
                   format="raster", overwrite=T)
  }
  raster::projection(area_p) <- raster::projection(env_uncut)
  return(area_p)
}

### TODO - update description
#' Cut projection areas of multiple species based on a unique SpatialPolygon (e.g. Ecoregion)
#'
#' This function will use a single SpatialPolygon to crop/mask raster/brick objects to be used on model projections.
#' @param area_p SpatialPolygon to be used as reference to crop/mask environmental variables of all species
#' @inheritParams pred.a
#' @inheritParams pred.ab
#' @inheritParams pred.a.polyb
#' @param mask Should mask raster? (i.e. only use area inside polygon. See ?raster::mask for details) or use all spatial extent of area_p
#' @return list of environmental layers (raster/brick) cutted
# TODO - examples
# #' @examples
#'
#' @export
pred.ab.rst <- function(area_p, env_uncut, occ_polys, mask=F, prj.nm="", sp.nm = "a.proj4mult.spp"){
  if(prj.nm != ""){ prj.nm <- paste0(".", prj.nm)}
  area_pl <- vector("list", length(occ_polys))
  names(area_pl) <- names(occ_polys)

  if(mask==F){
    area_p <- methods::as(raster::extent(area_p), "SpatialPolygons")
  }

  # for(i in 1:length(occ_polys)){
  #   cat(c("\n","Creating projection area for", names(occ_polys)[i],"\n",  "Species",  i, "of", length(occ_polys),"\n"))
  #   area_pl[[i]] <- pred.a.rst(area_p, env_uncut, occ_polys = occ_polys[[i]], mask=mask, prj.nm = prj.nm, sp.nm = names(occ_polys)[i])
  # }

  cat(c("\n","Creating projection area","\n"))
  # area_pl[[1]] <- pred.a.rst(area_p, env_uncut, mask=mask, prj.nm = prj.nm, sp.nm = sp.nm)
  area_pl[[1]] <- pred.a(area_p, env_uncut, prj.nm = prj.nm, sp.nm = sp.nm)
  for(i in 2:length(occ_polys)){
    # cat(c("\n","Creating projection area for", names(occ_polys)[i],"\n",  "Species",  i, "of", length(occ_polys),"\n"))
    area_pl[[i]] <- area_pl[[1]]
  }
  return(area_pl)
}

#### 4.8.3 function to cut multiple environmental layers based on pred_polys
### TODO - update description
#' Cut multiple projection areas of multiple species based on a unique SpatialPolygon (e.g. Ecoregion)
#'
#' This function will use a single SpatialPolygon to crop/mask multiple raster/brick objects
#' to be used on model projections.
#' @inheritParams pred.ab.mscn
#' @inheritParams pred.a.polyb
#' @inheritParams pred.ab.rst
#' @return list of list with multiple environmental layers (raster/brick) cutted
# TODO - examples
# #' @examples
#'
#' @export
pred.ab.rst.mscn <- function(area_p, env_uncut.l, occ_polys, mask=F, prj.nm="", sp.nm = "a.proj4mult.spp", cores=1){
  # if(prj.nm != ""){ prj.nm <- paste0(".", prj.nm)}
  path.proj <- "2_envData/area_proj"
  if(dir.exists(path.proj)==F) dir.create(path.proj)
  area_pl <- vector("list", length(occ_polys))
  names(area_pl) <- names(occ_polys)

  if(mask==F){
    area_p <- methods::as(raster::extent(area_p), "SpatialPolygons")
  }

  # for(i in 1:length(occ_polys)){
  #   cat(c("\n","Creating projection area for", names(occ_polys)[i],"\n",  "Species",  i, "of", length(occ_polys),"\n"))
  #   area_p.spi <- vector("list", length(env_uncut.l))
  #
  #
  #   area_pl[[i]] <- unlist(parallel::mclapply(seq_along(area_p.spi), mc.cores = getOption("mc.cores", as.integer(cores)), function(j){
  #     prj.nm.j <- paste("", prj.nm, names(env_uncut.l)[j], sep=".")
  #     area_p.spi[[j]] <- pred.a.rst(area_p, env_uncut.l[[j]], occ_polys[[i]], mask=mask, prj.nm = prj.nm.j, sp.nm = names(occ_polys)[i])
  #     area_p.spi[[j]] <- stats::setNames(area_p.spi[j], paste0(prj.nm, ".", names(env_uncut.l)[j]) )
  #   } ) )
  #
  # }
  cat(c("\n","Creating projection area","\n"))
  area_p.spi <- vector("list", length(env_uncut.l))

  area_pl[[1]] <- unlist(parallel::mclapply(seq_along(area_p.spi), mc.cores = getOption("mc.cores", as.integer(cores)), function(j){
    prj.nm.j <- paste("", prj.nm, names(env_uncut.l)[j], sep=".")
    # area_p.spi[[j]] <- pred.a.rst(area_p, env_uncut.l[[j]], occ_polys[[1]], mask=mask, prj.nm = prj.nm.j, sp.nm = sp.nm)
    area_p.spi[[j]] <- pred.a(area_p, env_uncut.l[[j]], prj.nm = prj.nm.j, sp.nm = sp.nm)
    area_p.spi[[j]] <- stats::setNames(area_p.spi[j], paste0(prj.nm, ".", names(env_uncut.l)[j]) )
  } ) )

  for(i in 2:length(occ_polys)){
    # cat(c("\n","Creating projection area for", names(occ_polys)[i],"\n",  "Species",  i, "of", length(occ_polys),"\n"))
    area_pl[[i]] <- area_pl[[1]]
  }
  return(area_pl)
}

