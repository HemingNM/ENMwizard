# working on this. Replacing functions from 2.f.proj.prep.R

# 4. Model Fitting (Calibration) and Projection
# 4.1 save rasters onto which the model will be projected in an elements called "area.proj"
# 4.1.1 select area for projection based on the extent of occ points

# "f.a.proj" names changed to "cut_projarea"

#' Select area for projection based on the extent of occ points
#'
#' This function will create SpatialPolygons that will be used to crop/mask raster/brick objects to be used on model projections. It has several options
#' The user can crop a squared area with an extent larger than the extent of occ.poly. By default, the "extent increase"
#' is the maximum of latitudinal and longitudinal extent "max(abs(ext.proj[1]-ext.proj[2]), abs(ext.proj[3]-ext.proj[4]))".
#' The result is added to each side of the occ.poly extent. This may be changed by setting "mult", which will be multiplied
#' by the "extent increase". Latitudinal and longitudinal increase may also vary independently by setting "same=FALSE".
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
#'
#' @seealso \code{\link{set_projarea_b}}
#' @return  SpatialPolygons (enlarged occ.poly)
# #' @examples
#'
#' @export
set_projarea <- function(occ.poly, sp.nm="sp", deg.incr=NULL, mult=1, buffer=FALSE, same=TRUE){ #, o.path = "occ.poly"
  path.proj <- "1_sppData/area.proj.poly"
  if(dir.exists("1_sppData")==FALSE) dir.create("1_sppData")
  if(dir.exists(path.proj)==FALSE) dir.create(path.proj)

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
  } else {
    area.p <- raster::buffer(occ.poly, width = deg.incr*mult)
    cat("\n", paste("New extent"), "\n")
    print(raster::extent(area.p))
    raster::projection(area.p) <- raster::projection(occ.poly)
  }
  raster::shapefile(area.p, filename = paste(path.proj, paste0(sp.nm, ".a.proj", ".shp"), sep = "/" ), overwrite=TRUE)
  return(area.p)
}


#' Select area for projection based on the extent of occ points for multiple species
#'
#' This function is a wrapper for "cut_projarea". See ?cut_projarea
#' It works with a named list of occ.polys to delimit the projection area for each of the species.
#'
#' @param occ.polys list of occ.poly SpatialPolygons.  see ?set_projarea for details.
#' @inheritParams set_projarea
#'
#' @seealso \code{\link{set_projarea}}
#' @return  named list of SpatialPolygons (enlarged occ.poly)
# TODO - examples
#' @examples
#' \dontrun{
#' area.projection <- set_projarea_b(occ.polys, mult=.55, buffer=FALSE)
#' plot(area.projection[[1]][[1]])
#' }
#' @export
set_projarea_b <- function(occ.polys, deg.incr=NULL, mult=1, buffer=FALSE, same=TRUE){ # , o.path = "occ.poly"
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
    area.pl[[i]] <- set_projarea(occ.polys[[i]], deg.incr=deg.incr.i, mult=mult, sp.nm = names(occ.polys)[i], buffer = buffer) #, o.path = o.path
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
#' @seealso \code{\link{set_projarea}}, \code{\link{set_projarea_b}},
#' \code{\link{cut_projarea_b}}, \code{\link{cut_projarea_mscn_b}},
#' \code{\link{cut_projarea_rst}}, \code{\link{cut_projarea_rst_b}}, \code{\link{cut_projarea_rst_mscn_b}}
#' @return  raster or brick cropped based on SpatialPolygons
# TODO - examples
# #' @examples
#'
#' @keywords internal
cut_projarea <- function(pred.poly, env.uncut, prj.nm="", sp.nm="species"){
  path.proj <- "2_envData/area.proj"
  if(dir.exists(path.proj)==FALSE) dir.create(path.proj)
  if(dir.exists(paste(path.proj, sp.nm, sep="/"))==FALSE) dir.create(paste(path.proj, sp.nm, sep="/"))

  p <- sp::Polygon(rbind(c(sp::bbox(pred.poly)[1,1], sp::bbox(pred.poly)[2,1]),
                     c(sp::bbox(pred.poly)[1,1], sp::bbox(pred.poly)[2,2]),
                     c(sp::bbox(pred.poly)[1,2], sp::bbox(pred.poly)[2,2]),
                     c(sp::bbox(pred.poly)[1,2], sp::bbox(pred.poly)[2,1]),
                     c(sp::bbox(pred.poly)[1,1], sp::bbox(pred.poly)[2,1])))
  p <- sp::SpatialPolygons(list(sp::Polygons(list(p), ID = 1)))

  if(identical(pred.poly, p)){
    area.p <- raster::crop(env.uncut, p,
                           file= paste(path.proj, sp.nm, gsub("^\\.|\\.\\.\\.|\\.\\.", ".", paste0("areaProj.", sp.nm, ".", prj.nm,".grd")), sep="/"),
                           format="raster", overwrite=TRUE)
  } else {
    area.p <- raster::crop(env.uncut, p)
    area.p <- raster::mask(area.p, pred.poly,
                           file= paste(path.proj, sp.nm, gsub("^\\.|\\.\\.\\.|\\.\\.", ".", paste0("areaProj.", sp.nm, ".", prj.nm,".grd")), sep="/"),
                           format="raster", overwrite=TRUE)
  }
  return(area.p)
}



### TODO -  descriptions
#' Cut area for projection based on a list of SpatialPolygons
#'
#' This function is a wrapper for "cut_projarea". See ?cut_projarea
#' It works with a named list of pred.poly to delimit the projection area for each of the species.
#'
#' @inheritParams cut_projarea
#' @param pred.polys list of SpatialPolygons (usually of based on species occ points)
#'
#' #' @seealso \code{\link{set_projarea}}, \code{\link{set_projarea_b}},
#' \code{\link{cut_projarea}}, \code{\link{cut_projarea_mscn_b}},
#' \code{\link{cut_projarea_rst}}, \code{\link{cut_projarea_rst_b}}, \code{\link{cut_projarea_rst_mscn_b}}
#' @return  named list of cropped raster or brick
#' @examples
#' \dontrun{
#' area.projection <- cut_projarea_b(pred.polys, env.uncut)
#' plot(area.projection[[1]][[1]])
#' }
#' @export
cut_projarea_b <- function(pred.polys, env.uncut, prj.nm=""){ # pred.poly, env.uncut, prj.nm="", sp.nm="sp"
  if(prj.nm != ""){ prj.nm <- paste0(".", prj.nm)}
  area.pl <- vector("list", length(pred.polys))
  names(area.pl) <- names(pred.polys)
  for(i in 1:length(pred.polys)){
    cat(c("\n", "Creating projection area for", names(pred.polys)[i],"\n",  "Species",  i, "of", length(pred.polys),"\n"))
    area.pl[[i]] <- cut_projarea(pred.polys[[i]], env.uncut, prj.nm = prj.nm, sp.nm = names(pred.polys)[i])
  }
  return(area.pl)
}


#### 4.8.3 function to cut multiple environmental layers based on pred.polys
### TODO - update description
#' Cut multiple projection areas (climatic scenarios) for multiple species (list of SpatialPolygons)
#'
#' This function is a wrapper for "cut_projarea". See ?cut_projarea. This function delimits the projection area for
#' each of the species contained in the pred.polys named list and crops
#' multiple rasters/bricks (i.e. representing distinct climatic scenaries) based on the same criteria for each species.
#'
#' @inheritParams cut_projarea_b
#' @param env.uncut.l list of raster/brick of environmental variables to be cut
#' @param numCores specify number of cores if aim run in parallel
#'
#' @seealso \code{\link{set_projarea}}, \code{\link{set_projarea_b}},
#' \code{\link{cut_projarea}}, \code{\link{cut_projarea_b}},
#' \code{\link{cut_projarea_rst}}, \code{\link{cut_projarea_rst_b}}, \code{\link{cut_projarea_rst_mscn_b}}
# #' @param ext.proj
#' @return  named list of cropped list of raster/brick of environmental variables
# TODO - examples
# #' @examples
#'
#' @export
cut_projarea_mscn_b <- function(pred.polys, env.uncut.l, numCores=1){ # , ext.proj=NULL, prj.nm=""
  # if(prj.nm != ""){ prj.nm <- paste0(".", prj.nm)}
  path.proj <- "2_envData/area.proj"
  if(dir.exists(path.proj)==FALSE) dir.create(path.proj)
  area.pl <- vector("list", length(pred.polys))
  names(area.pl) <- names(pred.polys)
  for(i in base::seq_along(pred.polys)){
    area.p.spi <- vector("list", length(env.uncut.l))
    if(numCores>1){
      check_install_pkg("parallel")

      area.pl[[i]] <- unlist(parallel::mclapply(base::seq_along(area.p.spi), mc.cores = getOption("mc.cores", as.integer(numCores)),
                                                function(j, pred.polys, env.uncut.l, i){
        prj.nm.j <- names(env.uncut.l)[j]
        area.p.spi[[j]] <- cut_projarea(pred.polys[[i]], env.uncut.l[[j]], prj.nm = prj.nm.j, sp.nm = names(pred.polys)[i]) # ext.proj,
        area.p.spi[[j]] <- stats::setNames(area.p.spi[j], gsub("^\\.", "", prj.nm.j)) # paste0(prj.nm, ".", names(env.uncut.l)[j]) )
      }, pred.polys, env.uncut.l, i ) )
    } else {
      area.pl[[i]] <- unlist(lapply(base::seq_along(area.p.spi), function(j, pred.polys, env.uncut.l, i){
        prj.nm.j <- names(env.uncut.l)[j]
        area.p.spi[[j]] <- cut_projarea(pred.polys[[i]], env.uncut.l[[j]], prj.nm = prj.nm.j, sp.nm = names(pred.polys)[i]) # ext.proj,
        area.p.spi[[j]] <- stats::setNames(area.p.spi[j], gsub("^\\.", "", prj.nm.j)) #paste0(prj.nm, ".", names(env.uncut.l)[j]) )
      }, pred.polys, env.uncut.l, i ) )
    }
  }
  return(area.pl)
}

### TODO - update description
### TODO - check this function became very similar to cut_projarea
#' Cut a projection area based on a SpatialPolygon (e.g. Ecoregion)
#'
#' This function will use a single SpatialPolygon to crop/mask raster/brick objects to be used on model projections.
#'
#' @param area.p SpatialPolygon to be used as reference to crop/mask environmental variables of all species
#' @inheritParams cut_projarea
#' @inheritParams cut_projarea_b
#' @param mask should use area.p to "mask" or "crop" env.uncut? See ?raster::mask and ?raster::crop for details
#' @param sp.mask an object of class: SpatialPolygons or SpatialPolygonsDataFrame, used as additional mask to be
#'  applied on the predictor variables (env.uncut)
#'
#' @seealso \code{\link{set_projarea}}, \code{\link{set_projarea_b}},
#' \code{\link{cut_projarea}}, \code{\link{cut_projarea_b}}, \code{\link{cut_projarea_mscn_b}},
#' \code{\link{cut_projarea_rst_b}}, \code{\link{cut_projarea_rst_mscn_b}}
#' @return environmental layers (raster/brick) cutted
# TODO - examples
# #' @examples
#'
#' @export
cut_projarea_rst <- function(area.p, env.uncut, mask=FALSE, sp.mask=NULL, prj.nm="", sp.nm="species"){ # , crs.set
  path.proj <- "2_envData/area.proj"
  if(dir.exists(path.proj)==FALSE) dir.create(path.proj)
  if(dir.exists(paste(path.proj, sp.nm, sep="/"))==FALSE) dir.create(paste(path.proj, sp.nm, sep="/"))

  ext.proj <- raster::extent(area.p)
  if(mask==FALSE){
    area.p <- raster::crop(env.uncut, ext.proj,
                           file = paste(path.proj, sp.nm, gsub("^\\.|\\.\\.\\.|\\.\\.", ".", paste0("areaProj.", sp.nm, prj.nm,".grd")), sep="/"),
                           format = "raster", overwrite=TRUE)
  } else {
    env.crp <- raster::crop(env.uncut, ext.proj)
    if(!is.null(sp.mask)){
      if(grepl("SpatialPolygons", class(sp.mask))){
        env.crp <- raster::mask(env.crp, sp.mask)
      } else {
        stop("sp.mask need to be an object of class: SpatialPolygons or SpatialPolygonsDataFrame")
      }
    }
    area.p <- raster::mask(env.crp, area.p,
                           file = paste(path.proj, sp.nm, gsub("^\\.|\\.\\.\\.|\\.\\.", ".", paste0("areaProj.", sp.nm, prj.nm,".grd")), sep="/"),
                           format = "raster", overwrite=TRUE)
  }
  raster::projection(area.p) <- raster::projection(env.uncut)
  return(area.p)
}

### TODO - update description
#' Cut projection areas of multiple species based on a single SpatialPolygon object (e.g. Ecoregion)
#'
#' This function will use a single SpatialPolygon to crop/mask raster/brick objects to be used on model projections.
#'
#' @param area.p SpatialPolygon to be used as reference to crop/mask environmental variables of all species
#' @inheritParams cut_projarea
#' @inheritParams cut_projarea_b
#' @inheritParams set_projarea_b
#' @param mask Should mask raster? (i.e. only use area inside polygon. See ?raster::mask for details) or
#' use all spatial extent of area.p
#' @param sp.mask.l a list of objects of class: SpatialPolygons or SpatialPolygonsDataFrame, used as additional mask to be
#'  applied on the predictor variables (env.uncut)
#'
#' #' @seealso \code{\link{set_projarea}}, \code{\link{set_projarea_b}},
#' \code{\link{cut_projarea}}, \code{\link{cut_projarea_b}}, \code{\link{cut_projarea_mscn_b}},
#' \code{\link{cut_projarea_rst}}, \code{\link{cut_projarea_rst_mscn_b}}
#' @return list of environmental layers (raster/brick) cutted
# TODO - examples
# #' @examples
#'
#' @export
cut_projarea_rst_b <- function(area.p, env.uncut, occ.polys, mask=FALSE, sp.mask.l=NULL, prj.nm="", sp.nm = "mult.spp"){
  if(prj.nm != ""){ prj.nm <- paste0(".", prj.nm)}
  area.pl <- vector("list", length(occ.polys))
  names(area.pl) <- names(occ.polys)

  if(mask==FALSE){
    area.p <- methods::as(raster::extent(area.p), "SpatialPolygons")
  }

  if(is.null(sp.mask.l)){
    cat(c("\n","Creating projection area","\n"))
    area.pl[[1]] <- cut_projarea(area.p, env.uncut, prj.nm = prj.nm, sp.nm = sp.nm)
    for(i in 2:length(occ.polys)){
      area.pl[[i]] <- area.pl[[1]]
    }
  } else {
    if(length(sp.mask.l) != length(occ.polys)){
      stop("length of 'sp.mask.l' must be the same of 'occ.polys'")
    }
    for(i in 1:length(occ.polys)){
      area.pl[[i]] <- cut_projarea_rst(area.p, env.uncut, mask=mask,
                                       sp.mask=sp.mask.l[[i]],
                                       prj.nm=prj.nm, sp.nm=sp.nm)
    }
  }
  return(area.pl)
}

#### 4.8.3 function to cut multiple environmental layers based on pred.polys
### TODO - update description
#' Cut multiple projection areas of multiple species based on a single SpatialPolygon object (e.g. Ecoregion)
#'
#' This function will use a single SpatialPolygon to crop/mask multiple raster/brick objects
#' to be used on model projections.
#'
#' @inheritParams cut_projarea_mscn_b
#' @inheritParams set_projarea_b
#' @inheritParams cut_projarea_rst_b
#' @inheritParams calib_mdl_b
#'
#' @seealso \code{\link{set_projarea}}, \code{\link{set_projarea_b}},
#' \code{\link{cut_projarea}}, \code{\link{cut_projarea_b}}, \code{\link{cut_projarea_mscn_b}},
#' \code{\link{cut_projarea_rst}}, \code{\link{cut_projarea_rst_b}}
#' @return list of list with multiple environmental layers (raster/brick) cutted
# TODO - examples
# #' @examples
#'
#' @export
cut_projarea_rst_mscn_b <- function(area.p, env.uncut.l, occ.polys, mask=FALSE,
                                    sp.mask.l=NULL, prj.nm="", sp.nm = "mult.spp", numCores=1){
  path.proj <- "2_envData/area.proj"
  if(dir.exists(path.proj)==FALSE) dir.create(path.proj)
  area.pl <- vector("list", length(occ.polys))
  names(area.pl) <- names(occ.polys)

  if(mask==FALSE){
    area.p <- methods::as(raster::extent(area.p), "SpatialPolygons")
  }

  cat(c("\n","Creating projection area","\n"))
  area.p.spi <- vector("list", length(env.uncut.l))
  check_install_pkg("parallel")

  if(is.null(sp.mask.l)){
    p.area <- unlist(parallel::mclapply(base::seq_along(env.uncut.l), mc.cores = getOption("mc.cores", as.integer(numCores)), function(j){
      prj.nm.j <- paste("", names(env.uncut.l)[j], sep=".")
      area.p.spi[[j]] <- cut_projarea(area.p, env.uncut.l[[j]], prj.nm = prj.nm.j, sp.nm = sp.nm)
      area.p.spi[[j]] <- stats::setNames(area.p.spi[j], gsub("^\\.|\\.\\.\\.|\\.\\.","", paste0(names(env.uncut.l)[j])) )
    } ) )

    for(i in 1:length(occ.polys)){
      area.pl[[i]] <- p.area
    }
  } else {
    for(i in 1:length(occ.polys)){
      area.pl[[i]] <- unlist(parallel::mclapply(base::seq_along(env.uncut.l), mc.cores = getOption("mc.cores", as.integer(numCores)), function(j){
        prj.nm.j <- paste("", names(env.uncut.l)[j], sep=".")

        area.p.spi[[j]] <-# cut_projarea(area.p, env.uncut.l[[j]], prj.nm = prj.nm.j, sp.nm = names(occ.polys)[i])
          cut_projarea_rst(area.p, env.uncut.l[[j]],
                           mask=mask, sp.mask=sp.mask.l[[i]],
                           prj.nm=prj.nm.j, sp.nm=names(occ.polys)[i])
        area.p.spi[[j]] <- stats::setNames(area.p.spi[j], gsub("^\\.|\\.\\.\\.|\\.\\.","", paste0(names(env.uncut.l)[j])) )
      } ) )
    }
  }
  return(area.pl)
}

