# NOT working on this. Functions are being replaced by those in 2.f.proj.prepb.R

# 4. Model Fitting (Calibration) and Projection
# 4.1 save rasters onto which the model will be projected in an elements called "area_projection"
# 4.1.1 select area for projection based on the extent of occ points
# "pred.a" names changed to "pred.a"

#' Select area for projection based on the extent of occ points
#'
#' This function will crop/mask raster/brick objects to be used on model projections. It has several options
#' The user can crop a squared area with an extent larger than the extent of occ_poly. By default, the "extent increase"
#' is the maximum of latitudinal and longitudinal extent "max(abs(ext_proj[1]-ext_proj[2]), abs(ext_proj[3]-ext_proj[4]))".
#' The result is added to each side of the occ_poly extent. This may be changed by setting "mult", which will be multiplied
#' by the "extent increase". Latitudinal and longitudinal increase may also vary independently by setting "same=F".
#'
#' The user can also set a buffer around occ_poly to cut the raster/brick object. Buffer value is, by default,
#' calculated in the same way as "extent increase", using "max(abs(ext_proj[1]-ext_proj[2]), abs(ext_proj[3]-ext_proj[4]))".
#' But the exact value can be defined through "deg.incr" and "mult". This method takes longer to run.
#'
#' @param occ_poly SpatialPolygon (usually of based on species occ points)
#' @param env_uncut raster of brick of environmental variables
#' @param sp.nm name (of species) to give to saved object
#' @param crs.set
#' @param mult how much increase or decrease buffer
#' @param buffer should the area be cut using a buffer around occ_poly?
#' @param deg.incr used to manually set the increase in prediction area, relative to occ_poly
#' @param same should latitudinal and longitudinal increase vary independently?
#' @return  raster or brick cropped based on occ points
#' @examples
#'
#' @export
pred.a <- function(occ_poly, env_uncut, prj.nm="", sp.nm="sp", crs.set, mult=1, buffer=F, deg.incr=NULL, same=T){
  path.proj <- "2_envData/area_projection"
  if(dir.exists(path.proj)==F) dir.create(path.proj)
  if(dir.exists(paste(path.proj, sp.nm, sep="/"))==F) dir.create(paste(path.proj, sp.nm, sep="/"))
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
      # increase the area to project by lat and lon range separately
      for(i in 1:length(as.vector(ext_proj))){
        ifelse(grepl("y", methods::slotNames(ext_proj)[i]), deg.incr <- deg.incr.y*mult, deg.incr <- deg.incr.x*mult)
        cat(paste("Change in", methods::slotNames(ext_proj)[i], ":", deg.incr, "\n"))
        ifelse(i %% 2 == 0, ext_proj[i] <- ext_proj[i] + deg.incr, ext_proj[i] <- ext_proj[i] - deg.incr)
      }
    }
    cat("\n", paste("New extent"), "\n")
    print(ext_proj)
    area_p <- raster::crop(env_uncut, ext_proj,
                   file= paste(path.proj, sp.nm, paste0("areaProj.", sp.nm, prj.nm,".grd"), sep="/"),
                   format="raster", overwrite=T)
  } else {
    area_p <- rgeos::gBuffer(occ_poly, width = deg.incr*mult, quadsegs=100)
    ext_proj <- raster::extent(area_p)
    cat("\n", paste("New extent"), "\n")
    print(ext_proj)
    env_crp <- raster::crop(env_uncut, ext_proj)#,
    area_p <- raster::mask(env_crp, area_p,
                   file= paste(path.proj, sp.nm, paste0("areaProj.", sp.nm, prj.nm,".grd"), sep="/"),
                   format="raster", overwrite=T)
  }
  raster::projection(area_p) <- raster::projection(env_uncut)
  return(area_p)
}




#' Select area for projection based on the extent of occ points
#'
#' This function is a wrapper for "pred.a". See ?pred.a
#' It works with a named list of occ_polys to delimitate the projection area for each of the species
#' @inheritParams pred.a
#' @param occ_polys list of SpatialPolygons (usually of based on species occ points)
#' @param prj.nm
#' @return  named list of cropped raster or brick
#' @examples
#'area_projection <- pred.ab(occ_polys, env_uncut, mult=.55, buffer=F)
#'plot(area_projection[[1]][[1]])
#' @export
pred.ab <- function(occ_polys, env_uncut, prj.nm="", deg.incr=NULL, mult=1, buffer=F, same=T, ...){
  if(prj.nm != ""){ prj.nm <- paste0(".", prj.nm)}
  area_pl <- vector("list", length(occ_polys))
  names(area_pl) <- names(occ_polys)
  # cat("\n", deg.incr, "\n")
  if(!is.null(deg.incr) & (length(deg.incr) != length(occ_polys))) { deg.incr <- rep(deg.incr, length(occ_polys))}
  if(length(deg.incr) != length(occ_polys)) {stop("length of deg.incr must be of same of occ_polys")}
  # cat("\n", deg.incr, "\n")
  for(i in 1:length(occ_polys)){
    cat(c("\n", "Creating projection area for", names(occ_polys)[i],"\n",  "Species",  i, "of", length(occ_polys),"\n"))
    area_pl[[i]] <- pred.a(occ_polys[[i]], env_uncut, prj.nm = prj.nm, deg.incr=deg.incr[i], mult=mult, sp.nm = names(occ_polys)[i], buffer = buffer)
  }
  return(area_pl)
}

#### 4.8.3 function to cut multiple environmental layers based on occ_polys
### TODO - option to run without "parallel"
#' Select area for projection based on the extent of occ points
#'
#' This function is a wrapper for "pred.a". See ?pred.a. This function delimitates the projection area for
#'  each of the species contained in the occ_polys named list and crops
#'  multiple rasters/bricks (i.e. representing distinct climatic scenaries) based on the same criteria for each species
#' @inheritParams pred.ab
#' @param prj.nm climatic scenario name, usually "fut" or "past" is used to indicate a general timing of scenario
#' @param env_crop.l
#' @param ext_proj
#' @param cores
#' @return  named list of cropped raster or brick
pred.ab.mscn <- function(occ_polys, env_crop.l, ext_proj=NULL, prj.nm="", mult=1, buffer=F, same=T, cores=1, ...){
  # if(prj.nm != ""){ prj.nm <- paste0(".", prj.nm)}
  path.proj <- "2_envData/area_projection"
  if(dir.exists(path.proj)==F) dir.create(path.proj)
  area_pl <- vector("list", length(occ_polys))
  names(area_pl) <- names(occ_polys)
  for(i in seq_along(occ_polys)){
    area_p.spi <- vector("list", length(env_crop.l))
    # names(area_p.spi) <- names(env_crop.l)
    # for(j in seq_along(area_p.spi)){
    if(cores>1){
      area_pl[[i]] <- unlist(parallel::mclapply(seq_along(area_p.spi), mc.cores = getOption("mc.cores", as.integer(cores)), function(j){
        prj.nm.j <- gsub("[..]",".",paste("", prj.nm, names(env_crop.l)[j], sep="."))
        area_p.spi[[j]] <- pred.a(occ_polys[[i]], env_crop.l[[j]], prj.nm = prj.nm.j, mult = mult, sp.nm = names(occ_polys)[i], buffer = buffer) # ext_proj,
        area_p.spi[[j]] <- setNames(area_p.spi[j], paste0(prj.nm, ".", names(env_crop.l)[j]) )
      } ) )
    } else {
      area_pl[[i]] <- unlist(lapply(seq_along(area_p.spi), function(j){
        prj.nm.j <- gsub("[..]",".",paste("", prj.nm, names(env_crop.l)[j], sep="."))
        area_p.spi[[j]] <- pred.a(occ_polys[[i]], env_crop.l[[j]], prj.nm = prj.nm.j, mult = mult, sp.nm = names(occ_polys)[i], buffer = buffer) #ext_proj,
        area_p.spi[[j]] <- setNames(area_p.spi[j], paste0(prj.nm, ".", names(env_crop.l)[j]) )
      } ) )
    }

    # area_pl[[i]] <- area_p.spi
  }
  return(area_pl)
}


pred.a.rst <- function(area_p, env_uncut, occ_polys, mask=F, prj.nm="", sp.nm="sp", crs.set){
  library(raster)
  library(rgeos)
  path.proj <- "2_envData/area_projection"
  if(dir.exists(path.proj)==F) dir.create(path.proj)
  if(dir.exists(paste(path.proj, sp.nm, sep="/"))==F) dir.create(paste(path.proj, sp.nm, sep="/"))
  # pensar melhor como projetar
  #if(is.na(projection(occ_polys))){print} # projection(occ_polys) <- crs.set
  #if(is.na(projection(env_uncut))){} # projection(env_uncut) <- crs.set
  cat(paste("Old extent"), "\n")
  print(extent(occ_polys))
  cat("\n")
  ext_proj <- extent(area_p)
  if(mask==F){
    cat("\n", paste("New extent"), "\n")
    print(ext_proj)
    area_p <- crop(env_uncut, ext_proj,
                   file= paste(path.proj, sp.nm, paste0("areaProj.", sp.nm, prj.nm,".grd"), sep="/"),
                   format="raster", overwrite=T)
  } else {
    # area_p <- gBuffer(area_p, width = deg.incr*mult, quadsegs=100)
    # ext_proj <- extent(area_p)
    cat("\n", paste("New extent"), "\n")
    print(ext_proj)
    env_crp <- crop(env_uncut, ext_proj)#,
    area_p <- mask(env_crp, area_p,
                   file= paste(path.proj, sp.nm, paste0("areaProj.", sp.nm, prj.nm,".grd"), sep="/"),
                   format="raster", overwrite=T)
    # }
  }
  projection(area_p) <- projection(env_uncut)
  return(area_p)
}

pred.ab.rst <- function(area_p, env_uncut, occ_polys, mask=F, prj.nm="", ...){
  if(prj.nm != ""){ prj.nm <- paste0(".", prj.nm)}
  area_pl <- vector("list", length(occ_polys))
  names(area_pl) <- names(occ_polys)

  for(i in 1:length(occ_polys)){
    cat(c("\n","Creating projection area for", names(occ_polys)[i],"\n",  "Species",  i, "of", length(occ_polys),"\n"))
    area_pl[[i]] <- pred.a.rst(area_p, env_uncut, occ_polys = occ_polys[[i]], mask=mask, prj.nm = prj.nm, sp.nm = names(occ_polys)[i])
  }
  return(area_pl)
}

#### 4.8.3 function to cut multiple environmental layers based on occ_polyss

pred.ab.rst.mscn <- function(area_p, env_crop.l, occ_polys, mask=F, prj.nm="", cores=2, ...){
  # if(prj.nm != ""){ prj.nm <- paste0(".", prj.nm)}
  path.proj <- "2_envData/area_projection"
  if(dir.exists(path.proj)==F) dir.create(path.proj)
  area_pl <- vector("list", length(occ_polys))
  names(area_pl) <- names(occ_polys)
  for(i in 1:length(occ_polys)){
    cat(c("\n","Creating projection area for", names(occ_polys)[i],"\n",  "Species",  i, "of", length(occ_polys),"\n"))
    area_p.spi <- vector("list", length(env_crop.l))


    area_pl[[i]] <- unlist(mclapply(seq_along(area_p.spi), mc.cores = getOption("mc.cores", as.integer(cores)), function(j){
      prj.nm.j <- paste("", prj.nm, names(env_crop.l)[j], sep=".")
      # area_pl[[i]] <- pred.a.rst(area_p, env_uncut, occ_polys[[i]], mask=mask, prj.nm = prj.nm, sp.nm = names(occ_polys)[i])
      area_p.spi[[j]] <- pred.a.rst(area_p, env_crop.l[[j]], occ_polys[[i]], mask=mask, prj.nm = prj.nm.j, sp.nm = names(occ_polys)[i])
      area_p.spi[[j]] <- setNames(area_p.spi[j], paste0(prj.nm, ".", names(env_crop.l)[j]) )
    } ) )

  }
  return(area_pl)
}




