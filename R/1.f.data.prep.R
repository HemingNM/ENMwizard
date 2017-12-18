#' Create polygon based on species occurence coordinates
#'
#' @param occ.spdf A spatial data.frame of coordinates, usually species occurence coordinates
# @param o.path Output path
#' @param lr.nm Polygon output name
#' @param convex Concave or convex polygon (T or F)
# #' @param alpha see ?alphahull::ashape
#' @param crs.set set the coordinate reference system (CRS) of the polygons
#' @inheritParams alphahull::ashape
#' @return spatial polygon of occurencies built using coordinates
#' @examples
#' occ.poly <- poly.c(occ.spdf, lr.nm="occ.poly")
#' plot(occ.poly)
#' @export
poly.c <- function(occ.spdf, lr.nm="sp.nm", convex=T, alpha=10, crs.set = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"){ # , o.path = NULL
  o.path <- "1_sppData/occ.poly"
  if(dir.exists("1_sppData")==F) dir.create("1_sppData")
  if(dir.exists(o.path)==F) dir.create(o.path)


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
  occ.poly <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(coords)), ID=1)))
  occ.poly <- sp::SpatialPolygonsDataFrame(occ.poly, data=data.frame(ID=1))
  raster::crs(occ.poly) <- raster::crs(crs.set)
  # if(!is.null(crs.set)){raster::projection(occ.poly) <- crs.set}
  # if(!is.null(o.path)){
    raster::shapefile(occ.poly, filename = paste(o.path, paste0(lr.nm,".shp"), sep = "/" ), overwrite=TRUE)
  # }
  return(occ.poly)
}


#' Create occ polygon for several species
#'
#' @param spp.occ.list A named list of spatial data.frame of species occurence points
#' @inheritParams poly.c
#' @param plot logical. Plot results or not?
#' @return A named list of spatial polygons built using coordinates
#' @examples
#' Bvarieg.occ <- read.table(paste(system.file(package="dismo"), "/ex/bradypus.csv", sep=""), header=TRUE, sep=",")
#' colnames(Bvarieg.occ) <- c("SPEC", "LONG", "LAT")
#' spp.occ.list <- list(Bvarieg = Bvarieg.occ)
#' occ.polys <- poly.c.batch(spp.occ.list)
#' occ.polys <- poly.c.batch(spp.occ.list, convex=T, alpha=10)
#' @export
poly.c.batch <- function(spp.occ.list, convex=T, alpha=10, plot=T, crs.set = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"){ #, o.path=NULL
  occ.pgns <- vector("list", length(spp.occ.list)) # , names=
  lr.nm <- paste(names(spp.occ.list), "occ.poly", sep = ".")

  o.path.pts <- "1_sppData/occ.pts"
  if(dir.exists("1_sppData")==F) dir.create("1_sppData")
  if(dir.exists(o.path.pts)==F) dir.create(o.path.pts)
  #
  for(i in 1:length(spp.occ.list)){
    occ.spdf <- data.frame(spp.occ.list[[i]])
    sp::coordinates(occ.spdf) <- ~LONG+LAT
    raster::shapefile(occ.spdf , filename = paste(o.path.pts, paste0(paste(names(spp.occ.list), "occ.pts", sep = "."),".shp"), sep = "/" ), overwrite=TRUE)
    # raster::crs(occ.spdf) <- crs.set
    # if(!is.null(crs.set)){raster::projection(occ.spdf) <- crs.set}
    occ.pgns[[i]] <- poly.c(occ.spdf, lr.nm=lr.nm[i], convex=convex, alpha=alpha, crs.set=crs.set) # , o.path=o.path
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
bind.shp <- function(files, sp.nm="sp", crs.set = NULL){ # , o.path = "occ.poly"
  o.path <- "1_sppData/occ.poly"
  if(dir.exists("1_sppData")==F) dir.create("1_sppData")
  if(dir.exists(o.path)==F) dir.create(o.path)

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
  sp.nm <- paste(sp.nm, "occ.poly", sep = ".")
  filename <- paste(o.path, paste0(sp.nm,".shp"), sep = "/" )

  # writeOGR(poly.data, dsn=o.path, layer=paste0(sp.nm), overwrite_layer=T, driver="ESRI Shapefile")
  # return(rgdal::readOGR(paste(o.path, paste0(sp.nm, ".shp"), sep="/")) )
  raster::shapefile(poly.data, filename = filename, overwrite=TRUE)
  # return(raster::shapefile(x = filename))
  return(poly.data)
}

#' split a species occ polygon (whenever distribution seems disjoint) into K polygons and save in a single .shp
#'
#' @param spp.occ species occurence coordinates
#' @param k number of polygons to create based on coordinates
#' @inheritParams poly.c
#' @inheritParams bind.shp
#' @return spatial polygons built using coordinates
#' @examples
#' Bvarieg.occ <- read.table(paste(system.file(package="dismo"), "/ex/bradypus.csv", sep=""), header=TRUE, sep=",")
#' colnames(Bvarieg.occ) <- c("SPEC", "LONG", "LAT")
#' spp.occ.list <- list(Bvarieg = Bvarieg.occ)
#' occ.polys <- poly.c.batch(spp.occ.list)
#' occ.polys$Bvarieg <- poly.splt(spp.occ = spp.occ.list$Bvarieg, k=5)
#' @export
poly.splt <- function(spp.occ, k=2, convex=T, alpha=10, sp.nm = "sp1", crs.set = NULL){ # , o.path = "occ.poly"
  # o.path <- "1_sppData/occ.poly"
  # if(dir.exists("1_sppData")==F) dir.create("1_sppData")
  # if(dir.exists(o.path)==F) dir.create(o.path)

  hc <- stats::hclust(stats::dist(cbind(spp.occ$LONG, spp.occ$LAT)))
  # plot(hc)
  clust <- stats::cutree(hc, k)

  # create one polygon for each set of points
  spp.k.list <- lapply(1:k, function(i){spp.occ[clust==i,]})
  occ.polys.lst <- poly.c.batch(spp.k.list, convex=convex, alpha=alpha, plot=F)
  occ.polys.sp <- bind.shp(occ.polys.lst, sp.nm = sp.nm, crs.set = crs.set) # , o.path = o.path
  # raster::crs(occ.polys.sp) <- crs.set
  # if(!is.null(crs.set)){raster::projection(occ.polys.sp) <- crs.set}
  sp::plot(occ.polys.sp)
  spp.occ <- as.data.frame(spp.occ)
  sp::coordinates(spp.occ) <- ~LONG+LAT
  sp::plot(spp.occ, col="red", add=T)
  return(occ.polys.sp)
}





#' Create buffer based on species polygons
#'
#' @param occ.polys list of SpatialPolygons, usually obj returned from poly.c.batch()
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
#' occ.polys <- poly.c.batch(spp.occ.list)
#' occ.b <- bffr.batch(occ.polys, bffr.width=1.5)
#' @export
bffr.batch <- function(occ.polys, bffr.width = NULL, mult = .2, quadsegs = 100, crs.set = NULL, plot = T){ # , o.path = "occ.poly"
  o.path <- "1_sppData/occ.bffr"
  if(dir.exists("1_sppData")==F) dir.create("1_sppData")
  if(dir.exists(o.path)==F) dir.create(o.path)

  # https://gis.stackexchange.com/questions/194848/creating-outside-only-buffer-around-polygon-using-r
  occ.b <- vector("list", length(occ.polys))
  names(occ.b) <- names(occ.polys)
  if(dir.exists(o.path)==F) dir.create(o.path)
  bf.path <- paste(o.path,"bffr", sep = "/" )
  if(dir.exists(bf.path)==F) dir.create(bf.path)

  if(length(mult)==1){ mult <- rep(mult, length(occ.polys))}
  TF.b.w <- is.null(bffr.width)
  for(i in 1:length(occ.polys)){
    if(!is.null(crs.set) & is.null(raster::projection(occ.polys[[i]]))){raster::projection(occ.polys[[i]]) <- crs.set}
    if(TF.b.w){
      ext.proj <- raster::extent(occ.polys[[i]])
      bffr.width <- mean(c(abs(ext.proj[1]-ext.proj[2]), abs(ext.proj[3]-ext.proj[4])))*mult[i]
    }
    cat(c("Buffer width for", names(occ.b)[i], "is", bffr.width, "\n"))
    occ.b[[i]] <- rgeos::gBuffer(occ.polys[[i]], width=bffr.width, quadsegs=quadsegs)
    # raster::crs(occ.b[[i]]) <- crs.set
    # if(!is.null(crs.set)){raster::projection(occ.b[[i]]) <- crs.set}
    raster::shapefile(occ.b[[i]], filename = paste(o.path, paste0(names(occ.b)[i], ".bffr", ".shp"), sep = "/" ), overwrite=TRUE)
    occ.b[[i]] <- raster::shapefile(paste(o.path, paste0(names(occ.b)[i], ".bffr", ".shp"), sep = "/" ))

    if(plot == T){
      plot(occ.b[[i]], col="green", main=names(occ.b)[i])
      plot(occ.polys[[i]], border="red", add=T)
    }
  }
  return(occ.b)
}




#' Crop environmental variables for each species
#'
#' @param occ.b polygon, usually a buffer
#' @param env.uncut raster brick or stack to be cropped
#' @return list [for each species] of cropped environmental variables. Details in ?raster::crop
#' @examples
#' Bvarieg.occ <- read.table(paste(system.file(package="dismo"), "/ex/bradypus.csv", sep=""), header=TRUE, sep=",")
#' colnames(Bvarieg.occ) <- c("SPEC", "LONG", "LAT")
#' spp.occ.list <- list(Bvarieg = Bvarieg.occ)
#' occ.polys <- poly.c.batch(spp.occ.list)
#' occ.b <- bffr.batch(occ.polys, bffr.width=1.5)
#' env.uncut <- brick("path/to/env")
#' occ.b.env <- env.cut(occ.b, env.uncut)
#' for(i in 1:length(occ.b.env)){
#'    plot(occ.b.env[[i]][[1]])
#'    plot(occ.b[[i]], add=T)
#' }
#' @export
env.cut <- function(occ.b, env.uncut){
  path.env.out <- "2_envData"
  if(dir.exists(path.env.out)==F) dir.create(path.env.out)
  ## Clipping rasters for each species
  occ.b.env <- vector("list", length(occ.b))
  names(occ.b.env) <- names(occ.b)

  for(i in 1:length(occ.b)){
    cat(c("Cutting environmental variables of species", i, "of", length(occ.b), "\n"))
    # occ.b[[i]] <- sp::spTransform(occ.b[[i]], raster::crs(env.uncut))
    occ.b.env[[i]] <- raster::crop(env.uncut, raster::extent(occ.b[[i]]))
    occ.b.env[[i]] <- raster::mask(occ.b.env[[i]], occ.b[[i]])
    # raster::crs(occ.b.env[[i]]) <- raster::crs(env.uncut)
    # if(dir.exists(paste("2_envData", names(spp.occ.list)[i], sep = "/") )==F) dir.create(paste("2_envData", names(spp.occ.list)[i], sep = "/"))
    occ.b.env[[i]] <- raster::writeRaster(occ.b.env[[i]],
                                  filename = paste(path.env.out, paste("envData.", names(occ.b.env)[i], ".grd", sep=''), sep='/'),
                                  format = "raster", overwrite = T)
    # plot(occ.b.env[[i]])
  }

  return(occ.b.env)
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
#' @return Named list containing thinned datasets for each species. See ?thin of spThin package. Also, by default it saves log file and the first thinned dataset in the folder "occ.thinned.full".
#' @examples
#' thinned.dataset.batch <- thin.batch(loc.data.lst = spp.occ.list)
#' plotThin(thinned.dataset.batch[[1]])
#' length(thinned.dataset.batch[[1]])
#' @export
thin.batch <- function(loc.data.lst, lat.col = "LAT", long.col = "LONG", spec.col = "SPEC",
                       thin.par = 10, reps = 10, locs.thinned.list.return = TRUE,
                       write.files = TRUE, max.files = 1, write.log.file = TRUE) {

  out.dir <- "1_sppData/occ.thinned.full"
  if(dir.exists("1_sppData")==F) dir.create("1_sppData")
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
                 out.base = paste0(spp[i], ".occ.thinned"),
                 log.file = paste0(out.dir, "/", spp[i], ".occ.thinned.full.log.file.txt"),
                 write.log.file = write.log.file)
  }

  thinned.dataset.full <- vector(mode = "list", length = length(loc.data.lst))
  thinned.dataset.full <- lapply(1:length(loc.data.lst), t.loc, loc.data.lst=loc.data.lst, spp=spp )

  names(thinned.dataset.full) <- spp

  return(thinned.dataset.full)
}




#' Load "thin.batch" filtered occurrence data
#'
#' @param occ.list.thin named list returned from "thin.batch"
#' @param from.disk boolean. Read from disk or from one of thinned datasets stored in 'occ.list.thin' obj
#' @param wtd : which thinned dataset?
#' @return named list of thinned occurence data for each species in the list
#' @examples
#' occ.locs <- loadTocc(thinned.dataset.batch)
#' @export
loadTocc <- function(occ.list.thin, from.disk=F, wtd=1){

  if (from.disk){ # retrieve from disk
    out.dir <- "1_sppData/occ.thinned.full"

    occ.l <- vector("list", length(occ.list.thin))
    names(occ.l) <- names(occ.list.thin)
    for(i in 1:length(occ.list.thin)){
      occ.l[[i]] <- utils::read.csv(paste0(out.dir, "/", names(occ.list.thin)[i], ".occ.thinned", ".thin1.csv"),
                             header=TRUE, sep=',', stringsAsFactors=F)[2:3]
    }
  } else { # retrieve from thinned obj
    for(i in 1:length(occ.list.thin)){
      if(wtd > length(occ.list.thin[[i]])) stop(paste("There are only", length(occ.list.thin[[i]]), "thinned datasets. 'wtd' was", wtd))

      occ.l[[i]] <- occ.list.thin[[i]][[wtd]]
      colnames(occ.l[[i]]) <- c("LONG", "LAT")
    }
  }

  return(occ.l)
}
