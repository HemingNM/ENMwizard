# working on this. Replacing functions from 2.f.proj.prep.R

# 4. Model Fitting (Calibration) and Projection
# 4.1 save rasters onto which the model will be projected in an elements called "area.proj"
# 4.1.1 select area for projection based on the extent of occ points

# "f.a.proj" names changed to "pred.a"

#' Select area for projection based on the extent of occ points
#'
#' This function will create SpatialPolygons that will be used to crop/mask raster/brick objects to be used on model projections. It has several options
#' The user can crop a squared area with an extent larger than the extent of occ.poly. By default, the "extent increase"
#' is the maximum of latitudinal and longitudinal extent "max(abs(ext.proj[1]-ext.proj[2]), abs(ext.proj[3]-ext.proj[4]))".
#' The result is added to each side of the occ.poly extent. This may be changed by setting "mult", which will be multiplied
#' by the "extent increase". Latitudinal and longitudinal increase may also vary independently by setting "same=F".
#'
#' The user can also set a buffer around occ.poly to cut the raster/brick object. Buffer value is, by default,
#' calculated in the same way as "extent increase", using "max(abs(ext.proj[1]-ext.proj[2]), abs(ext.proj[3]-ext.proj[4]))".
#' But the exact value can be defined through "deg.incr" and "mult". This method takes longer to run.
#'
#' @param occ.poly list of species SpatialPolygons
#' @param sp.nm name (of species) to give to saved object
# @param o.path Output path
#' @param deg.incr used to manually set the increase in prediction area, relative to occ.poly
#' @param mult how much increase or decrease buffer
#' @param buffer should the area be cut using a buffer around occ.poly?
#' @param same should latitudinal and longitudinal increase vary independently?
#' @seealso \code{\link{pred.a.poly.batch}}
#' @return  SpatialPolygons (enlarged occ.poly)
# #' @examples
#'
#' @export
pred.a.poly <- function(occ.poly, sp.nm="sp", deg.incr=NULL, mult=1, buffer=F, same=T){ #, o.path = "occ.poly"
  path.proj <- "1_sppData/area.proj.poly"
  if(dir.exists("1_sppData")==F) dir.create("1_sppData")
  # if(dir.exists("1_sppData/area.proj")==F) dir.create("1_sppData/area.proj")
  if(dir.exists(path.proj)==F) dir.create(path.proj)
  # if(dir.exists(paste(path.proj, sp.nm, sep="/"))==F) dir.create(paste(path.proj, sp.nm, sep="/"))

  # pensar melhor como projetar
  #if(is.na(projection(occ.poly))){print} # projection(occ.poly) <- crs.set
  #if(is.na(projection(env.uncut))){} # projection(env.uncut) <- crs.set

  ext.proj <- raster::extent(occ.poly)
  cat(paste("Old extent"), "\n")
  print(ext.proj)
  cat("\n")
  if(is.null(deg.incr)){
    deg.incr <- max(c(abs(ext.proj[1]-ext.proj[2]), abs(ext.proj[3]-ext.proj[4])))*mult
    #deg.incr <- mean(c(abs(ext.proj[1]-ext.proj[2]), abs(ext.proj[3]-ext.proj[4])))*mult
  }
  cat(c("Buffer width:", deg.incr, "\n"))

  if(buffer == F){
    if(same == T){# increase the area to project by max() range of lat and lon
      for(i in 1:length(as.vector(ext.proj))){
        # cat(paste("Change in",slotNames(ext.proj)[i], ":", deg.incr, "\n"))
        ifelse(i %% 2 == 0, ext.proj[i] <- ext.proj[i] + deg.incr, ext.proj[i] <- ext.proj[i] - deg.incr)
      }
    } else {
      deg.incr.x <- (abs(ext.proj[1]-ext.proj[2]))
      deg.incr.y <- (abs(ext.proj[3]-ext.proj[4]))
      # increase the area to project by lat and lon range separately
      for(i in 1:length(as.vector(ext.proj))){
        ifelse(grepl("y", methods::slotNames(ext.proj)[i]), deg.incr <- deg.incr.y*mult, deg.incr <- deg.incr.x*mult)
        cat(paste("Change in", methods::slotNames(ext.proj)[i], ":", deg.incr, "\n"))
        ifelse(i %% 2 == 0, ext.proj[i] <- ext.proj[i] + deg.incr, ext.proj[i] <- ext.proj[i] - deg.incr)
      }
    }
    cat("\n", paste("New extent"), "\n")
    print(ext.proj)
    area.p <- methods::as(ext.proj, "SpatialPolygons")
    raster::projection(area.p) <- raster::projection(occ.poly)
    # area.p <- raster::shapefile(paste(o.path, paste0(sp.nm, ".sel.area", ".shp"), sep = "/" ))
  } else {
    area.p <- rgeos::gBuffer(occ.poly, width = deg.incr*mult, quadsegs=100)
    cat("\n", paste("New extent"), "\n")
    print(raster::extent(area.p))
    raster::projection(area.p) <- raster::projection(occ.poly)
  }
  raster::shapefile(area.p, filename = paste(path.proj, paste0(sp.nm, ".a.proj", ".shp"), sep = "/" ), overwrite=TRUE)
  # area.p <- raster::shapefile(paste(path.proj, paste0(sp.nm, ".a.proj", ".shp"), sep = "/" ))
  return(area.p)
}


#' Select area for projection based on the extent of occ points
#'
#' This function is a wrapper for "pred.a". See ?pred.a
#' It works with a named list of occ.polys to delimitate the projection area for each of the species.
#'
#' @param occ.polys list of occ.poly SpatialPolygons.  see ?pred.a.poly for details.
#' @inheritParams pred.a.poly
#' @seealso \code{\link{pred.a.poly}}
#' @return  named list of SpatialPolygons (enlarged occ.poly)
# TODO - examples
#' @examples
#'area.projection <- pred.a.poly.batch(occ.polys, mult=.55, buffer=F)
#'plot(area.projection[[1]][[1]])
#' @export
pred.a.poly.batch <- function(occ.polys, deg.incr=NULL, mult=1, buffer=F, same=T){ # , o.path = "occ.poly"
  area.pl <- vector("list", length(occ.polys))
  names(area.pl) <- names(occ.polys)
  # cat("\n", deg.incr, "\n")
  if(!is.null(deg.incr) & length(deg.incr) == 1) { deg.incr <- rep(deg.incr, length(occ.polys))}
  if(!is.null(deg.incr) & (length(deg.incr) != length(occ.polys))) {stop("length of deg.incr must be of same of occ.polys")}
  # cat("\n", deg.incr, "\n")
  for(i in 1:length(occ.polys)){
    if(is.null(deg.incr)){
      deg.incr.i <- NULL
    }else{
      deg.incr.i <- deg.incr[i]
    }
    cat(c("\n", "Creating projection area for", names(occ.polys)[i],"\n",  "Species",  i, "of", length(occ.polys),"\n"))
    area.pl[[i]] <- pred.a.poly(occ.polys[[i]], deg.incr=deg.incr.i, mult=mult, sp.nm = names(occ.polys)[i], buffer = buffer) #, o.path = o.path
  }
  return(area.pl)
}


#' Cut area for projection based on SpatialPolygons
#'
#' This function will use SpatialPolygons to crop/mask raster/brick objects to be used on model projections.
#'
#' @param pred.poly list of SpatialPolygons (usually of based on species occ points)
#' @param env.uncut raster/brick of environmental variables to be cut
#' @param prj.nm climatic scenario name, usually "fut" or "past" is used to indicate a general timing of scenario
#' @param sp.nm name (of species) to give to saved object
#' @seealso \code{\link{pred.a.poly}}, \code{\link{pred.a.poly.batch}},
#' \code{\link{pred.a.batch}}, \code{\link{pred.a.batch.mscn}},
#' \code{\link{pred.a.rst}}, \code{\link{pred.a.batch.rst}}, \code{\link{pred.a.batch.rst.mscn}}
#' @return  raster or brick cropped based on SpatialPolygons
# TODO - examples
# #' @examples
#'
#' @export
pred.a <- function(pred.poly, env.uncut, prj.nm="", sp.nm="sp"){
  path.proj <- "2_envData/area.proj"
  if(dir.exists(path.proj)==F) dir.create(path.proj)
  if(dir.exists(paste(path.proj, sp.nm, sep="/"))==F) dir.create(paste(path.proj, sp.nm, sep="/"))

  # ext.proj <- raster::extent(pred.poly)
  # add polygon
  # x.min<-bbox(pred.poly)[1,1]
  # x.max<-bbox(pred.poly)[1,2]
  # y.min<-bbox(pred.poly)[2,1]
  # y.max<-bbox(pred.poly)[2,2]
  # x = c(x.min, x.min, x.max, x.max, x.min)
  # y = c(y.min, y.max, y.max, y.min, y.min)
  # # m <- matrix(c(x,y), ncol = 2, byrow = F)
  # # pm <- sapply(slot(pred.poly, "polygons"), function(x) lapply(slot(x,"Polygons"), function(y) slot(y, "coords")))[[1]]
  # p <- Polygon(matrix(c(x,y), ncol = 2, byrow = F))
  p <- sp::Polygon(rbind(c(sp::bbox(pred.poly)[1,1], sp::bbox(pred.poly)[2,1]),
                     c(sp::bbox(pred.poly)[1,1], sp::bbox(pred.poly)[2,2]),
                     c(sp::bbox(pred.poly)[1,2], sp::bbox(pred.poly)[2,2]),
                     c(sp::bbox(pred.poly)[1,2], sp::bbox(pred.poly)[2,1]),
                     c(sp::bbox(pred.poly)[1,1], sp::bbox(pred.poly)[2,1])))
  p <- sp::SpatialPolygons(list(sp::Polygons(list(p), ID = 1)))
  # # p <- methods::as(ext.proj, "SpatialPolygons")

  # if(all.equal(pm, m)){
  if(identical(pred.poly, p)){
  # if(mask==F){
    area.p <- raster::crop(env.uncut, p,
                           file= paste(path.proj, sp.nm, paste0("areaProj.", sp.nm, prj.nm,".grd"), sep="/"),
                           format="raster", overwrite=T)
  } else {
    env.crp <- raster::crop(env.uncut, p)
    area.p <- raster::mask(env.crp, pred.poly,
                           file= paste(path.proj, sp.nm, paste0("areaProj.", sp.nm, prj.nm,".grd"), sep="/"),
                           format="raster", overwrite=T)
  }
  return(area.p)
}



### TODO -  descriptions
#' Cut area for projection based on a list of SpatialPolygons
#'
#' This function is a wrapper for "pred.a". See ?pred.a
#' It works with a named list of pred.poly to delimitate the projection area for each of the species.
#'
#' @inheritParams pred.a
#' @param pred.polys list of SpatialPolygons (usually of based on species occ points)
#' #' @seealso \code{\link{pred.a.poly}}, \code{\link{pred.a.poly.batch}},
#' \code{\link{pred.a}}, \code{\link{pred.a.batch.mscn}},
#' \code{\link{pred.a.rst}}, \code{\link{pred.a.batch.rst}}, \code{\link{pred.a.batch.rst.mscn}}
#' @return  named list of cropped raster or brick
#' @examples
#' area.projection <- pred.a.batch(pred.polys, env.uncut)
#' plot(area.projection[[1]][[1]])
#' @export
pred.a.batch <- function(pred.polys, env.uncut, prj.nm=""){ # pred.poly, env.uncut, prj.nm="", sp.nm="sp"
  if(prj.nm != ""){ prj.nm <- paste0(".", prj.nm)}
  area.pl <- vector("list", length(pred.polys))
  names(area.pl) <- names(pred.polys)
  for(i in 1:length(pred.polys)){
    cat(c("\n", "Creating projection area for", names(pred.polys)[i],"\n",  "Species",  i, "of", length(pred.polys),"\n"))
    area.pl[[i]] <- pred.a(pred.polys[[i]], env.uncut, prj.nm = prj.nm, sp.nm = names(pred.polys)[i])
  }
  return(area.pl)
}


#### 4.8.3 function to cut multiple environmental layers based on pred.polys
### TODO - update description
#' Cut multiple projection areas (climatic scenarios) for multiple species (list of SpatialPolygons)
#'
#' This function is a wrapper for "pred.a". See ?pred.a. This function delimitates the projection area for
#' each of the species contained in the pred.polys named list and crops
#' multiple rasters/bricks (i.e. representing distinct climatic scenaries) based on the same criteria for each species.
#'
#' @inheritParams pred.a.batch
#' @param env.uncut.l list of raster/brick of environmental variables to be cut
#' @param numCores specify number of cores if aim run in parallel
#' @seealso \code{\link{pred.a.poly}}, \code{\link{pred.a.poly.batch}},
#' \code{\link{pred.a}}, \code{\link{pred.a.batch}},
#' \code{\link{pred.a.rst}}, \code{\link{pred.a.batch.rst}}, \code{\link{pred.a.batch.rst.mscn}}
# #' @param ext.proj
#' @return  named list of cropped list of raster/brick of environmental variables
# TODO - examples
# #' @examples
#'
#' @export
pred.a.batch.mscn <- function(pred.polys, env.uncut.l, prj.nm="", numCores=1){ # , ext.proj=NULL
  # if(prj.nm != ""){ prj.nm <- paste0(".", prj.nm)}
  path.proj <- "2_envData/area.proj"
  if(dir.exists(path.proj)==F) dir.create(path.proj)
  area.pl <- vector("list", length(pred.polys))
  names(area.pl) <- names(pred.polys)
  for(i in base::seq_along(pred.polys)){
    area.p.spi <- vector("list", length(env.uncut.l))
    # names(area.p.spi) <- names(env.uncut.l)
    # for(j in seq_along(area.p.spi)){
    if(numCores>1){
      area.pl[[i]] <- unlist(parallel::mclapply(base::seq_along(area.p.spi), mc.cores = getOption("mc.cores", as.integer(numCores)), function(j){
        prj.nm.j <- gsub("^\\.|\\.\\.\\.|\\.\\.",".",paste(prj.nm, names(env.uncut.l)[j], sep="."))
        area.p.spi[[j]] <- pred.a(pred.polys[[i]], env.uncut.l[[j]], prj.nm = prj.nm.j, sp.nm = names(pred.polys)[i]) # ext.proj,
        area.p.spi[[j]] <- stats::setNames(area.p.spi[j], gsub("^\\.", "", prj.nm.j)) # paste0(prj.nm, ".", names(env.uncut.l)[j]) )
      } ) )
    } else {
      area.pl[[i]] <- unlist(lapply(base::seq_along(area.p.spi), function(j){
        prj.nm.j <- gsub("^\\.|\\.\\.\\.|\\.\\.",".",paste(prj.nm, names(env.uncut.l)[j], sep="."))
        area.p.spi[[j]] <- pred.a(pred.polys[[i]], env.uncut.l[[j]], prj.nm = prj.nm.j, sp.nm = names(pred.polys)[i]) # ext.proj,
        area.p.spi[[j]] <- stats::setNames(area.p.spi[j], gsub("^\\.", "", prj.nm.j)) #paste0(prj.nm, ".", names(env.uncut.l)[j]) )
      } ) )
    }

    # area.pl[[i]] <- area.p.spi
  }
  return(area.pl)
}

### TODO - update description
### TODO - check this function became very similar to pred.a
#' Cut a projection area based on a SpatialPolygon (e.g. Ecoregion)
#'
#' This function will use a single SpatialPolygon to crop/mask raster/brick objects to be used on model projections.
#'
#' @param area.p SpatialPolygon to be used as reference to crop/mask environmental variables of all species
#' @inheritParams pred.a
#' @inheritParams pred.a.batch
#' @param mask should use area.p to "mask" or "crop" env.uncut? See ?raster::mask and ?raster::crop for details
#' @seealso \code{\link{pred.a.poly}}, \code{\link{pred.a.poly.batch}},
#' \code{\link{pred.a}}, \code{\link{pred.a.batch}}, \code{\link{pred.a.batch.mscn}},
#' \code{\link{pred.a.batch.rst}}, \code{\link{pred.a.batch.rst.mscn}}
#' @return environmental layers (raster/brick) cutted
# TODO - examples
# #' @examples
#'
#' @export
pred.a.rst <- function(area.p, env.uncut, mask=F, prj.nm="", sp.nm="sp"){ # , crs.set
  path.proj <- "2_envData/area.proj"
  if(dir.exists(path.proj)==F) dir.create(path.proj)
  if(dir.exists(paste(path.proj, sp.nm, sep="/"))==F) dir.create(paste(path.proj, sp.nm, sep="/"))

  ext.proj <- raster::extent(area.p)
  if(mask==F){
    area.p <- raster::crop(env.uncut, ext.proj,
                           file= paste(path.proj, sp.nm, paste0("areaProj.", sp.nm, prj.nm,".grd"), sep="/"),
                           format="raster", overwrite=T)
  } else {
    env.crp <- raster::crop(env.uncut, ext.proj)
    area.p <- raster::mask(env.crp, area.p,
                           file= paste(path.proj, sp.nm, paste0("areaProj.", sp.nm, prj.nm,".grd"), sep="/"),
                           format="raster", overwrite=T)
  }
  raster::projection(area.p) <- raster::projection(env.uncut)
  return(area.p)
}

### TODO - update description
#' Cut projection areas of multiple species based on a unique SpatialPolygon (e.g. Ecoregion)
#'
#' This function will use a single SpatialPolygon to crop/mask raster/brick objects to be used on model projections.
#'
#' @param area.p SpatialPolygon to be used as reference to crop/mask environmental variables of all species
#' @inheritParams pred.a
#' @inheritParams pred.a.batch
#' @inheritParams pred.a.poly.batch
#' @param mask Should mask raster? (i.e. only use area inside polygon. See ?raster::mask for details) or use all spatial extent of area.p
#' #' @seealso \code{\link{pred.a.poly}}, \code{\link{pred.a.poly.batch}},
#' \code{\link{pred.a}}, \code{\link{pred.a.batch}}, \code{\link{pred.a.batch.mscn}},
#' \code{\link{pred.a.rst}}, \code{\link{pred.a.batch.rst.mscn}}
#' @return list of environmental layers (raster/brick) cutted
# TODO - examples
# #' @examples
#'
#' @export
pred.a.batch.rst <- function(area.p, env.uncut, occ.polys, mask=F, prj.nm="", sp.nm = "a.proj4mult.spp"){
  if(prj.nm != ""){ prj.nm <- paste0(".", prj.nm)}
  area.pl <- vector("list", length(occ.polys))
  names(area.pl) <- names(occ.polys)

  if(mask==F){
    area.p <- methods::as(raster::extent(area.p), "SpatialPolygons")
  }

  # for(i in 1:length(occ.polys)){
  #   cat(c("\n","Creating projection area for", names(occ.polys)[i],"\n",  "Species",  i, "of", length(occ.polys),"\n"))
  #   area.pl[[i]] <- pred.a.rst(area.p, env.uncut, occ.polys = occ.polys[[i]], mask=mask, prj.nm = prj.nm, sp.nm = names(occ.polys)[i])
  # }

  cat(c("\n","Creating projection area","\n"))
  # area.pl[[1]] <- pred.a.rst(area.p, env.uncut, mask=mask, prj.nm = prj.nm, sp.nm = sp.nm)
  area.pl[[1]] <- pred.a(area.p, env.uncut, prj.nm = prj.nm, sp.nm = sp.nm)
  for(i in 2:length(occ.polys)){
    # cat(c("\n","Creating projection area for", names(occ.polys)[i],"\n",  "Species",  i, "of", length(occ.polys),"\n"))
    area.pl[[i]] <- area.pl[[1]]
  }
  return(area.pl)
}

#### 4.8.3 function to cut multiple environmental layers based on pred.polys
### TODO - update description
#' Cut multiple projection areas of multiple species based on a unique SpatialPolygon (e.g. Ecoregion)
#'
#' This function will use a single SpatialPolygon to crop/mask multiple raster/brick objects
#' to be used on model projections.
#'
#' @inheritParams pred.a.batch.mscn
#' @inheritParams pred.a.poly.batch
#' @inheritParams pred.a.batch.rst
#' @inheritParams mxnt.c.batch
#' @seealso \code{\link{pred.a.poly}}, \code{\link{pred.a.poly.batch}},
#' \code{\link{pred.a}}, \code{\link{pred.a.batch}}, \code{\link{pred.a.batch.mscn}},
#' \code{\link{pred.a.rst}}, \code{\link{pred.a.batch.rst}}
#' @return list of list with multiple environmental layers (raster/brick) cutted
# TODO - examples
# #' @examples
#'
#' @export
pred.a.batch.rst.mscn <- function(area.p, env.uncut.l, occ.polys, mask=F, prj.nm="", sp.nm = "a.proj4mult.spp", numCores=1){
  # if(prj.nm != ""){ prj.nm <- paste0(".", prj.nm)}
  path.proj <- "2_envData/area.proj"
  if(dir.exists(path.proj)==F) dir.create(path.proj)
  area.pl <- vector("list", length(occ.polys))
  names(area.pl) <- names(occ.polys)

  if(mask==F){
    area.p <- methods::as(raster::extent(area.p), "SpatialPolygons")
  }

  # for(i in 1:length(occ.polys)){
  #   cat(c("\n","Creating projection area for", names(occ.polys)[i],"\n",  "Species",  i, "of", length(occ.polys),"\n"))
  #   area.p.spi <- vector("list", length(env.uncut.l))
  #
  #
  #   area.pl[[i]] <- unlist(parallel::mclapply(seq_along(area.p.spi), mc.cores = getOption("mc.cores", as.integer(cores)), function(j){
  #     prj.nm.j <- paste("", prj.nm, names(env.uncut.l)[j], sep=".")
  #     area.p.spi[[j]] <- pred.a.rst(area.p, env.uncut.l[[j]], occ.polys[[i]], mask=mask, prj.nm = prj.nm.j, sp.nm = names(occ.polys)[i])
  #     area.p.spi[[j]] <- stats::setNames(area.p.spi[j], paste0(prj.nm, ".", names(env.uncut.l)[j]) )
  #   } ) )
  #
  # }
  cat(c("\n","Creating projection area","\n"))
  area.p.spi <- vector("list", length(env.uncut.l))

  area.pl[[1]] <- unlist(parallel::mclapply(base::seq_along(area.p.spi), mc.cores = getOption("mc.cores", as.integer(numCores)), function(j){
    prj.nm.j <- paste("", prj.nm, names(env.uncut.l)[j], sep=".")
    # area.p.spi[[j]] <- pred.a.rst(area.p, env.uncut.l[[j]], occ.polys[[1]], mask=mask, prj.nm = prj.nm.j, sp.nm = sp.nm)
    area.p.spi[[j]] <- pred.a(area.p, env.uncut.l[[j]], prj.nm = prj.nm.j, sp.nm = sp.nm)
    # area.p.spi[[j]] <- stats::setNames(area.p.spi[j], paste0(prj.nm, ".", names(env.uncut.l)[j]) )
    area.p.spi[[j]] <- stats::setNames(area.p.spi[j], gsub("^\\.|\\.\\.\\.|\\.\\.","", paste0(prj.nm, ".", names(env.uncut.l)[j])) )
  } ) )

  for(i in 2:length(occ.polys)){
    # cat(c("\n","Creating projection area for", names(occ.polys)[i],"\n",  "Species",  i, "of", length(occ.polys),"\n"))
    area.pl[[i]] <- area.pl[[1]]
  }
  return(area.pl)
}

