#' Create polygon based on species occurence coordinates
#'
#' @param occ.spdf A spatial data.frame of coordinates, usually species occurence coordinates
#' @param o.path Output path
#' @param lr.nm Polygon output name
#' @param convex Concave or convex polygon (T or F)
#' @param alpha see ?alphahull::ashape
#' @param crs.set crs value used ???
#' @return spatial polygon built using coordinates
#' @examples
#' occ_poly <- poly.c(occ.spdf, o.path="occ_poly", lr.nm="occ_poly")
#' plot(occ_poly)
#' @export
poly.c <- function(occ.spdf, o.path = NULL, lr.nm="sp_nm", convex=T, alpha=10, crs.set = NULL ){
  if(convex==F){ # convex hulls to crop rasters
    # http://r.789695.n4.nabble.com/Concave-hull-td863710.html#a4688606
    # https://rpubs.com/geospacedman/alphasimple

    ch <- alphahull::ashape(unique(sp::coordinates(occ.spdf)), alpha=alpha)
    chg <- igraph::graph.edgelist(cbind(as.character(ch$edges[, "ind1"]),
                               as.character(ch$edges[, "ind2"])), directed = FALSE)
    if (!igraph::is.connected(chg)) {
      stop("Graph not connected")
    } else if (any(igraph::degree(chg) != 2)) {
      stop("Graph not circular")
    } else if (igraph::clusters(chg)$no > 1) {
      stop("Graph composed of more than one circle")
    }
    # In the next step, we delete one edge of the circle to create a chain,
    # and find the path from one end of the chain to the other.
    cutg <- chg - igraph::E(chg)[1]
    # find chain end points
    ends <- names(which(igraph::degree(cutg) == 1))
    path <- igraph::get.shortest.paths(cutg, ends[1], ends[2])[[1]]
    # this is an index into the points
    pathX <- as.numeric(igraph::V(chg)[path[[1]]]$name)
    # join the ends
    pathX <- c(pathX, pathX[1])
    # now show the alpha shape plot with our poly on top
    # plot(ch, lwd = 10, col = "gray")
    # get the points from the ashape object
    # lines(ch$x[pathX, ], lwd = 2)
    coords <- ch$x[pathX, ]
  } else {
    ch <- grDevices::chull(sp::coordinates(occ.spdf))
    coords <- sp::coordinates(occ.spdf)[c(ch, ch[1]),]
  }
  occ_poly <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(coords)), ID=1)))
  occ_poly <- sp::SpatialPolygonsDataFrame(occ_poly, data=data.frame(ID=1))
  # raster::crs(occ_poly) <- crs.set
  # if(!is.null(crs.set)){raster::projection(occ_poly) <- crs.set}
  if(!is.null(o.path)){
    raster::shapefile(occ_poly, filename = paste(o.path, paste0(lr.nm,".shp"), sep = "/" ), overwrite=TRUE)
  }
  return(occ_poly)
}


#' Create occ polygon for several species
#'
#' @param spp.occ.list A named list of spatial data.frame of species occurence points
#' @inheritParams poly.c
#' @return A named list of spatial polygons built using coordinates
#' @examples
#' occ_polys <- poly.c.batch(spp.occ.list, o.path="occ_poly", convex=T, alpha=10)
#' @export
poly.c.batch <- function(spp.occ.list, o.path=NULL, crs.set = NULL, convex=T, alpha=10, plot=T){
  occ.pgns <- vector("list", length(spp.occ.list)) # , names=
  lr.nm <- paste(names(spp.occ.list), o.path, sep = ".")
  if(!is.null(o.path)) {if(dir.exists(o.path)==F) dir.create(o.path)}
  for(i in 1:length(spp.occ.list)){
    occ.spdf <- data.frame(spp.occ.list[[i]])
    sp::coordinates(occ.spdf) <- ~LONG+LAT
    # raster::crs(occ.spdf) <- crs.set
    # if(!is.null(crs.set)){raster::projection(occ.spdf) <- crs.set}
    occ.pgns[[i]] <- poly.c(occ.spdf, o.path=o.path, lr.nm=lr.nm[i], convex=convex, alpha=alpha, crs.set=crs.set)
    if(plot){
      sp::plot(occ.pgns[[i]], main=names(spp.occ.list)[i])
      sp::plot(occ.spdf, col="red", add=T)
    }
  }
  names(occ.pgns) <- names(spp.occ.list)
  return(occ.pgns)
}


#' bind list of SpatialPolygons into a single SpatialPolygon
#'
#' @param files ??
#' @param sp.nm name used to save file
#' @inheritParams poly.c
#' @return shapefile with binded polygons
#' @export
bind.shp <- function(files, sp.nm="sp", o.path = "occ_poly", crs.set = NULL ){
  # http://r-sig-geo.2731867.n2.nabble.com/merging-several-shapefiles-into-one-td6401613.html
  # Get polygons and change IDs
  uid<-1
  poly.l <- vector("list", length(files))
  for (i in 1:length(files)) {
    temp.data <- files[[i]]
    n <- length(methods::slot(temp.data, "polygons"))
    temp.data <- sp::spChFIDs(temp.data, as.character(uid:(uid+n-1)))
    uid <- uid + n
    poly.l[[i]] <- temp.data
    # poly.data <- spRbind(poly.data,temp.data)
  }

  # mapunit polygoan: combin remaining  polygons with first polygoan
  poly.data <- do.call(raster::bind, poly.l)
  # names(poly.data)
  # raster::crs(poly.data) <- crs.set
  # if(!is.null(crs.set)){raster::projection(poly.data) <- crs.set}
  sp.nm <- paste0(sp.nm, ".occ_poly")
  if(!is.null(o.path)){
    raster::shapefile(poly.data, filename = paste(o.path, paste0(sp.nm,".shp"), sep = "/" ), overwrite=TRUE)
    # writeOGR(poly.data, dsn=o.path, layer=paste0(sp.nm), overwrite_layer=T, driver="ESRI Shapefile")
    return(rgdal::readOGR(paste(o.path, paste0(sp.nm, ".shp"), sep="/")) )
  } else {
    return(poly.data)
  }
}

#' split a species occ polygon (whenever distribution seems disjoint) into K polygons and save in a single .shp
#'
#' @param spp.occ species occurence coordinates
#' @param k number of polygons to create based on coordinates
#' @inheritParams poly.c
#' @inheritParams bind.shp
#' @return spatial polygons built using coordinates
#' @examples
#' occ_polys$Bvarieg <- poly.splt(spp.occ = Bvarieg.occ, k=5)
#' @export
poly.splt <- function(spp.occ, k=2, convex=T, alpha=10, sp.nm = "sp1", o.path = "occ_poly", crs.set = NULL){
  hc <- stats::hclust(stats::dist(cbind(spp.occ$LONG, spp.occ$LAT)))
  # plot(hc)
  clust <- stats::cutree(hc, k)

  # create one polygon for each set of points
  spp.k.list <- lapply(1:k, function(i){spp.occ[clust==i,]})
  occ_polys.lst <- poly.c.batch(spp.k.list, convex=convex, alpha=alpha, plot=F)
  occ_polys.sp <- bind.shp(occ_polys.lst, sp.nm = sp.nm, o.path = o.path, crs.set = crs.set)
  # raster::crs(occ_polys.sp) <- crs.set
  # if(!is.null(crs.set)){raster::projection(occ_polys.sp) <- crs.set}
  sp::plot(occ_polys.sp)
  spp.occ <- as.data.frame(spp.occ)
  sp::coordinates(spp.occ) <- ~LONG+LAT
  sp::plot(spp.occ, col="red", add=T)
  return(occ_polys.sp)
}





#' Create buffer based on species polygons
#'
#' @param occ_polys list of SpatialPolygons, usually obj returned from poly.c.batch()
#' @param bffr.width Buffer width. See 'width' of ?rgeos::gBuffer
#' @param mult How much expand bffr.width
#' @param plot Boolean, to draw plots or not
#' @param quadsegs see ?rgeos::gBuffer
#' @inheritParams poly.c
#' @return A named list of SpatialPolygons
#' @examples
#' Bvarieg.occ <- read.table(paste(system.file(package="dismo"), "/ex/bradypus.csv", sep=""), header=TRUE, sep=",")
#' colnames(Bvarieg.occ) <- c("SPEC", "LONG", "LAT")
#' spp.occ.list <- list(Bvarieg = Bvarieg.occ)
#' occ_polys <- poly.c.batch(spp.occ.list, o.path="occ_poly")
#' occ_b <- bffr.b(occ_polys, bffr.width=1.5)
#' @export
bffr.b <- function(occ_polys, bffr.width = NULL, mult = .2, quadsegs = 100, o.path = "occ_poly", crs.set = NULL, plot = T){
  # https://gis.stackexchange.com/questions/194848/creating-outside-only-buffer-around-polygon-using-r
  occ_b <- vector("list", length(occ_polys))
  names(occ_b) <- names(occ_polys)
  if(dir.exists(o.path)==F) dir.create(o.path)
  bf.path <- paste(o.path,"bffr", sep = "/" )
  if(dir.exists(bf.path)==F) dir.create(bf.path)
  if(length(mult)==1){ mult <- rep(mult, length(occ_polys))}
  TF.b.w <- is.null(bffr.width)
  for(i in 1:length(occ_polys)){
    if(!is.null(crs.set) & is.null(raster::projection(occ_polys[[i]]))){raster::projection(occ_polys[[i]]) <- crs.set}
    if(TF.b.w){
      ext_proj <- raster::extent(occ_polys[[i]])
      bffr.width <- mean(c(abs(ext_proj[1]-ext_proj[2]), abs(ext_proj[3]-ext_proj[4])))*mult[i]
    }
    cat(c("Buffer width for", names(occ_b)[i], "is", bffr.width, "\n"))
    occ_b[[i]] <- rgeos::gBuffer(occ_polys[[i]], width=bffr.width, quadsegs=quadsegs)
    # raster::crs(occ_b[[i]]) <- crs.set
    # if(!is.null(crs.set)){raster::projection(occ_b[[i]]) <- crs.set}
    raster::shapefile(occ_b[[i]], filename = paste(o.path, "bffr", paste0(names(occ_b)[i], "_bffr", ".shp"), sep = "/" ), overwrite=TRUE)
    occ_b[[i]] <- raster::shapefile(paste(o.path, "bffr", paste0(names(occ_b)[i], "_bffr", ".shp"), sep = "/" ))

    if(plot == T){
      plot(occ_b[[i]], col="green", main=names(occ_b)[i])
      plot(occ_polys[[i]], border="red", add=T)
    }
  }
  return(occ_b)
}




#' Crop environmental variables for each species
#'
#' @param occ_b polygon, usually a buffer
#' @param env_uncut raster brick or stack to be cropped
#' @return list [for each species] of cropped environmental variables. Details in ?raster::crop
#' @examples
#' env_uncut <- brick("path/to/env")
#' occ_b_env <- env.cut(occ_b, env_uncut)
#' for(i in 1:length(occ_b_env)){
#'    plot(occ_b_env[[i]][[1]])
#'    plot(occ_b[[i]], add=T)
#' }
#' @export
env.cut <- function(occ_b, env_uncut){
  path.env.out <- "2_envData"
  if(dir.exists(path.env.out)==F) dir.create(path.env.out)
  ## Clipping rasters for each species
  occ_b_env <- vector("list", length(occ_b))
  names(occ_b_env) <- names(occ_b)

  for(i in 1:length(occ_b)){
    cat(c("Cutting environmental variables of species", i, "of", length(occ_b), "\n"))
    # occ_b[[i]] <- sp::spTransform(occ_b[[i]], raster::crs(env_uncut))
    occ_b_env[[i]] <- raster::crop(env_uncut, raster::extent(occ_b[[i]]))
    occ_b_env[[i]] <- raster::mask(occ_b_env[[i]], occ_b[[i]])
    # raster::crs(occ_b_env[[i]]) <- raster::crs(env_uncut)
    # if(dir.exists(paste("2_envData", names(spp.occ.list)[i], sep = "/") )==F) dir.create(paste("2_envData", names(spp.occ.list)[i], sep = "/"))
    occ_b_env[[i]] <- raster::writeRaster(occ_b_env[[i]],
                                  filename = paste(path.env.out, paste("envData.", names(occ_b_env)[i], ".grd", sep=''), sep='/'),
                                  format = "raster", overwrite = T)
    # plot(occ_b_env[[i]])
  }

  return(occ_b_env)
}



#' Filter each species' occurence dataset
#'
#' @param loc.data.lst Named list containing species occ data
#' @param lat.col See ?thin of spThin package
#' @param long.col See ?thin of spThin package
#' @param spec.col See ?thin of spThin package
#' @param thin.par See ?thin of spThin package
#' @param reps See ?thin of spThin package
#' @param locs.thinned.list.return See ?thin of spThin package
#' @param write.files See ?thin of spThin package
#' @param max.files See ?thin of spThin package
#' @param write.log.file See ?thin of spThin package
#' @return Named list containing thinned datasets for each species. See ?thin of spThin package. Also, by default it saves log file and the first thinned dataset in the folder "occ_thinned_full".
#' @examples
#' thinned_dataset_batch <- thin.batch(loc.data.lst = spp.occ.list)
#' plotThin(thinned_dataset_batch[[1]])
#' length(thinned_dataset_batch[[1]])
#' @export
thin.batch <- function(loc.data.lst, lat.col = "LAT", long.col = "LONG", spec.col = "SPEC",
                       thin.par = 10, reps = 10, locs.thinned.list.return = TRUE,
                       write.files = TRUE, max.files = 1, write.log.file = TRUE) {

  out.dir <- "occ_thinned_full"
  if(dir.exists(out.dir)==F) dir.create(out.dir)
  spp <- names(loc.data.lst)

  t.loc <- function(i, loc.data.lst,  spp, ...){
    spThin::thin(loc.data.lst[[i]],
                 lat.col = lat.col, long.col = long.col,
                 spec.col = spec.col,
                 thin.par = thin.par, reps = reps, # reps = 1000 thin.par 'é a distancia min (km) para considerar pontos distintos
                 locs.thinned.list.return = locs.thinned.list.return,
                 write.files = write.files,
                 max.files = max.files,
                 out.dir = out.dir,
                 out.base = paste0(spp[i], ".occ_thinned"),
                 log.file = paste0(out.dir, "/", spp[i], ".occ_thinned_full_log_file.txt"),
                 write.log.file = write.log.file)
  }

  thinned_dataset_full <- vector(mode = "list", length = length(loc.data.lst))
  thinned_dataset_full <- lapply(1:length(loc.data.lst), t.loc, loc.data.lst=loc.data.lst, spp=spp )

  names(thinned_dataset_full) <- spp

  return(thinned_dataset_full)
}




#' Load occurrence data filtered using "thin.batch"
#' @param occ.list.thin named list returned from "thin.batch"
#' @param from.disk boolean. Read from disk or from one of thinned datasets stored in 'occ.list.thin' obj
#' @param wtd : which thinned dataset?
#' @return named list of thinned occurence data for each species in the list
#' @examples
#' occ_locs <- loadTocc(thinned_dataset_batch)
#' @export
loadTocc <- function(occ.list.thin, from.disk=F, wtd=1){

  if (from.disk){ # retrieve from disk
    out.dir <- "occ_thinned_full"

    occ_l <- vector("list", length(occ.list.thin))
    names(occ_l) <- names(occ.list.thin)
    for(i in 1:length(occ.list.thin)){
      occ_l[[i]] <- utils::read.csv(paste0(out.dir, "/", names(occ.list.thin)[i], ".occ_thinned", "_thin1.csv"),
                             header=TRUE, sep=',', stringsAsFactors=F)[2:3]
    }
  } else { # retrieve from thinned obj
    for(i in 1:length(occ.list.thin)){
      if(wtd > length(occ.list.thin[[i]])) stop(paste("There are only", length(occ.list.thin[[i]]), "thinned datasets. 'wtd' was", wtd))

      occ_l[[i]] <- occ.list.thin[[i]][[wtd]]
      colnames(occ_l[[i]]) <- c("LONG", "LAT")
    }
  }

  return(occ_l)
}
