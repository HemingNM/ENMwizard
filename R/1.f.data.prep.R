#####------- 1. Prepare environmental data

# 1. Prepare enviromental layers
#' Create occ polygon to crop rasters prior to modelling
#'
#' @param occ.spdf A spatial data.frame of coordinates, usually species occurence coordinates
#' @param o.path Output path
#' @param lr.nm Polygon output name
#' @param convex Concave or convex polygon (T or F)
#' @return spatial polygon built using coordinates
#' @examples
#' occ_poly <- f.poly(occ.spdf, o.path="occ_poly", lr.nm="occ_poly")
#' plot(occ_poly)
f.poly <- function(occ.spdf, o.path = NULL, lr.nm="occ_poly", convex=T, alpha=10, crs.set = NA ){
  if(convex==F){ # convex hulls to crop rasters
    # http://r.789695.n4.nabble.com/Concave-hull-td863710.html#a4688606
    # https://rpubs.com/geospacedman/alphasimple

    ch <- alphahull::ashape(unique(coordinates(occ.spdf)), alpha=10)
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
    ch <- grDevices::chull(coordinates(occ.spdf))
    coords <- sp::coordinates(occ.spdf)[c(ch, ch[1]),]
  }
  occ_poly <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(coords)), ID=1)))
  occ_poly <- sp::SpatialPolygonsDataFrame(occ_poly, data=data.frame(ID=1))
  crs(occ_poly) <- crs.set
  if(!is.null(o.path)){
    raster::shapefile(occ_poly, filename = paste(o.path, paste0(lr.nm,".shp"), sep = "/" ), overwrite=TRUE)
  }
  return(occ_poly)
}


#' Create occ polygon for several species
#'
#' @param spp.occ.list A named list of spatial data.frame of species occurence points
#' @param o.path Output path
#' @param convex Concave or convex polygon (T or F)
#' @return A named list of spatial polygons built using coordinates
#' @examples
#' occ_polys <- f.poly.batch(spp.occ.list, o.path="occ_poly", convex=T, alpha=10, crs.set = crs.set)
f.poly.batch <- function(spp.occ.list, o.path=NULL, crs.set = NA, convex=T, alpha=10){
  occ.pgns <- vector("list", length(spp.occ.list)) # , names=
  lr.nm <- paste(names(spp.occ.list), o.path, sep = ".")
  if(!is.null(o.path)) {if(dir.exists(o.path)==F) dir.create(o.path)}
  for(i in 1:length(spp.occ.list)){
    occ.spdf <- data.frame(spp.occ.list[[i]])
    sp::coordinates(occ.spdf) <- ~LONG+LAT
    raster::crs(occ.spdf) <- crs.set
    occ.pgns[[i]] <- f.poly(occ.spdf, o.path=o.path, lr.nm=lr.nm[i], convex=convex, alpha=alpha, crs.set=crs.set)
    plot(occ.pgns[[i]], main=names(spp.occ.list)[i])
    plot(occ.spdf, col="red", add=T)
  }
  names(occ.pgns) <- names(spp.occ.list)
  return(occ.pgns)
}


#' merge several shapefiles into a single one
#'
#' @param files
#' @param sp.nm
#' @param o.path Output path
#' @param crs.set
#' @return
f.bind.shp <- function(files, sp.nm="sp", o.path = "occ_poly", crs.set = NA ){
  # http://r-sig-geo.2731867.n2.nabble.com/merging-several-shapefiles-into-one-td6401613.html
  # Require packages: rgdal and maptool
  #-------------------------------------
  library(rgdal)
  library(maptools)

  # Get polygons and change IDs
  #-------------------------------------
  uid<-1
  poly.l <- vector("list", length(files))
  for (i in 1:length(files)) {
    temp.data <- files[[i]]
    n <- length(slot(temp.data, "polygons"))
    temp.data <- spChFIDs(temp.data, as.character(uid:(uid+n-1)))
    uid <- uid + n
    poly.l[[i]] <- temp.data
    # poly.data <- spRbind(poly.data,temp.data)
  }

  # mapunit polygoan: combin remaining  polygons with first polygoan
  #-----------------------------------------------------------------
  poly.data <- do.call(bind, poly.l)
  # names(poly.data)
  crs(poly.data) <- crs.set
  sp.nm <- paste0(sp.nm, ".occ_poly")
  if(!is.null(o.path)){
    shapefile(poly.data, filename = paste(o.path, paste0(sp.nm,".shp"), sep = "/" ), overwrite=TRUE)
    # writeOGR(poly.data, dsn=o.path, layer=paste0(sp.nm), overwrite_layer=T, driver="ESRI Shapefile")
    return(readOGR(paste(o.path, paste0(sp.nm, ".shp"), sep="/")) )
  } else {
    return(poly.data)
  }
}

#' create N polygons for a species (whenever distribution seems disjoint) and save in a single .shp
#'
#' @param spp.occ species occurence coordinates
#' @param k number of polygons to create based on coordinates
#' @param convex boolean (T or F). Create concave or convex polygons
#' @param alpha
#' @param sp.nm output name
#' @param o.path Output path
#' @param crs.set
#' @return spatial polygons built using coordinates
#' @examples
#' occ_polys$Bvarieg <- f.poly.splt(spp.occ = Bvarieg.occ, k=5, convex=T, alpha=10, sp.nm = "Bvarieg", o.path = "occ_poly", crs.set = crs.set)
f.poly.splt <- function(spp.occ, k=2, convex=T, alpha=10, sp.nm = "sp1", o.path = "occ_poly", crs.set = NA){
  hc <- hclust(dist(cbind(spp.occ$LONG, spp.occ$LAT)))
  # plot(hc)
  clust <- cutree(hc, k)

  # create one polygon for each set of points
  spp.k.list <- lapply(1:k, function(i){spp.occ[clust==i,]})
  occ_polys.lst <- f.poly.batch(spp.k.list, convex=convex, alpha=alpha)
  occ_polys.sp <- f.bind.shp(occ_polys.lst, sp.nm = sp.nm, o.path = o.path, crs.set = crs.set)
  crs(occ_polys.sp) <- crs.set
  plot(occ_polys.sp)
  spp.occ <- as.data.frame(spp.occ)
  coordinates(spp.occ) <- ~LONG+LAT
  plot(spp.occ, col="red", add=T)
  return(occ_polys.sp)
}




#### 1.2 creating buffer
#' Create buffer based on polygon
#'
#' @param occ_polys Spatial polygon
#' @param bffr.width Buffer width. See 'width' of ?gBuffer
#' @param mult How much expand bffr.width
#' @param quadsegs
#' @param o.path Output path
#' @param crs.set
#' @param plot Boolean, to draw plots or not
#' @return A named list of SpatialPolygons
#' @examples
#' occ_b <- f.bffr(occ_polys, bffr.width=1.5, crs.set=crs.set) #
f.bffr <- function(occ_polys, bffr.width=NULL, mult=.2, quadsegs=100, o.path = "occ_poly", crs.set=NULL, plot=T){
  # https://gis.stackexchange.com/questions/194848/creating-outside-only-buffer-around-polygon-using-r
  library(rgeos)
  library(raster)
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
    crs(occ_b[[i]]) <- crs.set
    shapefile(occ_b[[i]], filename = paste(o.path, "bffr", paste0(names(occ_b)[i], "_bffr", ".shp"), sep = "/" ), overwrite=TRUE)
    occ_b[[i]] <- raster::shapefile(paste(o.path, "bffr", paste0(names(occ_b)[i], "_bffr", ".shp"), sep = "/" ))

    if(plot == T){
      plot(occ_b[[i]], col="green", main=names(occ_b)[i])
      plot(occ_polys[[i]], border="red", add=T)
    }
  }
  return(occ_b)
}



# # path to environmental variables
# path.env <- "/Volumes/Samsung SSD/neanderh/Documents/CloudStation/ArcGIS-Virtual Machine/Mapas/Clima e Biosfera/ClimaticScenaries/PR/WorldClim/2_5min/bio_2-5m_bil"
# biovars <- paste0("bio", 1:17)#c("bio5", "bio8", "bio10", "bio13", "bio16") # "bio18" tem problemas para o cerrado
# pattern.env = 'asc'
# path.env.out <- "3_envData"
#
# # 1.3. Cut enviromental layers with M and save in hardrive.
# # Get uncut variables
# env_uncut <- list.files(path.env, pattern = pattern.env, full.names=TRUE)
# env_uncut <- env_uncut[grepl(paste(paste0(biovars, ".", pattern.env), collapse = "|"), env_uncut)]
# env_uncut <- stack(env_uncut) #predictors_uncut
# crs(env_uncut) <- crs.set
#
# # # Se as variáveis estiverem prontas:
# env_uncut <- brick(paste(path.env, "bio.grd", sep="/"))
#

#' Crop environmental variables for each species
#'
#' @param occ_b polygon, usually a buffer
#' @param env_uncut raster brick or stack to be cropped
#' @return
#' @examples
#' 1.3.2 crop environmental variables for each species
#' occ_b_env <- f.cut.env(occ_b, env_uncut)
#' for(i in 1:length(occ_b_env)){
#'    plot(occ_b_env[[i]][[1]])
#'    plot(occ_b[[i]], add=T)
#'    }
f.cut.env <- function(occ_b, env_uncut){
  path.env.out <- "3_envData"
  ## Clipping rasters for each species
  occ_b_env <- vector("list", length(occ_b))
  names(occ_b_env) <- names(occ_b)
  if(dir.exists(path.env.out)==F) dir.create(path.env.out)

  for(i in 1:length(occ_b)){
    cat(c("Cutting environmental variables of species", i, "of", length(occ_b), "\n"))
    occ_b[[i]] <- sp::spTransform(occ_b[[i]], raster::crs(env_uncut))
    occ_b_env[[i]] <- raster::crop(env_uncut, raster::extent(occ_b[[i]]))
    occ_b_env[[i]] <- raster::mask(occ_b_env[[i]], occ_b[[i]])
    raster::crs(occ_b_env[[i]]) <- raster::crs(env_uncut)
    # if(dir.exists(paste("3_envData", names(spp.occ.list)[i], sep = "/") )==F) dir.create(paste("3_envData", names(spp.occ.list)[i], sep = "/"))
    occ_b_env[[i]] <- raster::writeRaster(occ_b_env[[i]],
                                  filename = paste(path.env.out, paste("envData.", names(occ_b_env)[i], ".grd", sep=''), sep='/'),
                                  format = "raster", overwrite = T)
    # plot(occ_b_env[[i]])
  }

  return(occ_b_env)
}

# occ_b_env <- f.cut.env(occ_b, env_uncut)
#
#
