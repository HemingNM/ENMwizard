#' Create minimum convex polygon based on species occurence data
#'
#' This function will create minimum convex polygon based on coordinates of species occurence data.
#'
#' @param occ.spdf An object of class SpatialPoints or SpatialPointsDataFrame (e.g. species occurence coordinates)
#' @param sp.nm Species name, used on saving shapefile
#' @param convex Logical. Convex (T) or concave (F) polygon
#' @inheritParams alphahull::ashape
#' @param save Should save polygons on disk?
#' @param crs.set set the coordinate reference system (CRS) of the polygons
#'
#' @seealso \code{\link{poly.c.batch}}, \code{\link{poly.splt}}
#' @return An object of class SpatialPolygons or SpatialPolygonsDataFrame. Polygon built using coordinates of species occurence data
#' @examples
#' occ.poly <- poly.c(occ.spdf, sp.nm="occ.poly")
#' plot(occ.poly)
#' @keywords internal
#' @export
poly.c <- function(occ.spdf, sp.nm="sp.nm", convex=TRUE, alpha=10, save=TRUE, crs.set = "+proj=longlat +datum=WGS84"){ # , o.path = NULL
  o.path <- "1_sppData/occ.poly"
  if(dir.exists("1_sppData")==FALSE) dir.create("1_sppData")
  if(dir.exists(o.path)==FALSE) dir.create(o.path)

  u.pts <- as.data.frame(unique(sp::coordinates(occ.spdf)))

  if(convex==FALSE){ # convex hulls to crop rasters
    # http://r.789695.n4.nabble.com/Concave-hull-td863710.html#a4688606
    # https://rpubs.com/geospacedman/alphasimple

    ch <- alphahull::ashape(u.pts, alpha=alpha)
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
    ch <- grDevices::chull(u.pts)
    coords <- u.pts[c(ch, ch[1]),]
  }
  occ.poly <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(coords)), ID = 1)))
  occ.poly <- sp::SpatialPolygonsDataFrame(occ.poly, data = data.frame(ID = 1, row.names = 1)) # , sp.nm = sp.nm
  raster::crs(occ.poly) <- raster::crs(crs.set)
  # if(!is.null(crs.set)){raster::projection(occ.poly) <- crs.set}
  # if(!is.null(o.path)){
  sp.nm <- paste(sp.nm, "occ.poly", sep = ".")
  filename <- paste(o.path, paste0(sp.nm,".shp"), sep = "/" )
  if(save){
    raster::shapefile(occ.poly, filename = filename, overwrite = TRUE)
  }
  # }
  return(occ.poly)
}


#' Split a species occ polygon (whenever distribution seems disjoint) into K polygons and save in a single .shp
#'
#' Cluster points and create several small polygons. Implemented methods are 'Hierarchical Clustering' (when 'k'
#' number of clusters is defined 'k > 0'), 'Elbow' (c.m = "E"), 'Affinity Propagation' (c.m = AP), and several methods
#' implemented in function NbClust of NbClust package (c.m = "NB"). To use NbClust package, check arguments 'distance',
#' 'min.nc', 'max.nc', 'method', and 'index' in ?NbClust::NbClust.
#'
#' @param k number of polygons to create based on coordinates
#' @param c.m clustering method to find the best number of clusters (k). Currently E (Elbow) or (Affinity Propagation).
#' @param nm.col.dt "character". Name of a numeric column to use as grouping variable in addition to coordinates.
#' @inheritParams poly.c
#' @inheritParams bind.shp
#' @inheritParams apcluster::negDistMat
#' @inheritParams apcluster::apcluster
#' @inheritParams NbClust::NbClust
#'
#' @seealso \code{\link{poly.c}}, \code{\link{poly.c.batch}}, \code{\link[NbClust]{NbClust}}
#' @return spatial polygons built using coordinates
#' @examples
#' Bvarieg.occ <- read.table(paste(system.file(package="dismo"),
#'  "/ex/bradypus.csv", sep=""), header=TRUE, sep=",")
#' colnames(Bvarieg.occ) <- c("SPEC", "LONG", "LAT")
#' spp.occ.list <- list(Bvarieg = Bvarieg.occ)
#' occ.polys <- poly.c.batch(spp.occ.list)
#' occ.polys$Bvarieg <- poly.splt(occ.spdf = spp.occ.list$Bvarieg, k=5)
#' @keywords internal
#' @export
poly.splt <- function(occ.spdf, k=NULL, nm.col.dt=NULL, c.m = "NB", r = 2, q = 0.3,
                      distance = "euclidean", min.nc = 1, max.nc = 20,
                      method = "centroid", index = "trcovw", alphaBeale = 0.1,
                      convex=TRUE, alpha=10, sp.nm = "sp.nm", save = T,
                      crs.set = "+proj=longlat +datum=WGS84"){ # , o.path = "occ.poly"

  if(is.null(nm.col.dt)){
    u.pts <- as.data.frame(unique(sp::coordinates(occ.spdf)))
  } else {
    if(!is.numeric(occ.spdf@data[,nm.col.dt])){
      stop(paste(nm.col.dt, "must be numeric!"))
    }
    u.pts <- as.data.frame(unique(cbind(sp::coordinates(occ.spdf), occ.spdf@data[,nm.col.dt])))
  }

  if(k == 0 | is.null(k)){
    # http://www.sthda.com/english/articles/29-cluster-validation-essentials/96-determining-the-optimal-number-of-clusters-3-must-know-methods/
    # https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
    # for a package for further improvements see:
    # https://cran.r-project.org/web/packages/clValid/vignettes/clValid.pdf

    # if(c.m == "E"){ ## ELBOW method
    #   dist.obj <- stats::dist(u.pts)
    #   hclust.obj <- stats::hclust(dist.obj)
    #   css.obj <- GMD::css.hclust(dist.obj, hclust.obj)
    #   elbow.obj <- GMD::elbow.batch(css.obj)
    #   k <- elbow.obj$k
    #   clust <- stats::cutree(hclust.obj, k)
    # } else
    if (c.m == "AP") { # Affinity Propagation (AP)
      apclus <- apcluster::apcluster(apcluster::negDistMat(r=r), u.pts)
      # apclus <- apcluster::apcluster(apcluster::expSimMat(r=2, w=10), u.pts)
      apclus <- apcluster::apcluster(apclus@sim, q=q)
      length(apclus@clusters)
      k <- length(apclus@clusters)
      clust <- apclus@idx
      clust <- as.factor(clust)
      levels(clust) <- 1:k
      clust <- as.numeric(clust)
    } else if (c.m == "NB") {
      if(nrow(u.pts)<max.nc){
        max.nc <- nrow(u.pts)-1
      }
      nb <- NbClust::NbClust(u.pts, distance = distance, min.nc = min.nc, max.nc = max.nc,
                             method = method, index = index, alphaBeale = 0.1)
      clust <- nb$Best.partition
      k <- length(unique(clust))
    }
  } else { # Hierarchical Clustering
    # https://stackoverflow.com/questions/28672399/spatial-clustering-in-r-simple-example
    hclust.obj <- stats::hclust(stats::dist(u.pts))
    clust <- stats::cutree(hclust.obj, k)
  }

  # create one polygon for each set of points
  sp::coordinates(u.pts) <- cbind(u.pts[,1], u.pts[,2])

  spp.k.list <- lapply(1:k, function(i){u.pts[clust==i,]})
  names(spp.k.list) <- paste0(sp.nm, seq_along(spp.k.list))
  occ.polys.lst <- poly.c.batch(spp.k.list, k=1, convex=convex, alpha=alpha, plot=FALSE, save=FALSE)
  if(length(occ.polys.lst)>1){
    occ.polys.sp <- bind.shp(occ.polys.lst, sp.nm = sp.nm, save=save, crs.set = crs.set) # , o.path = o.path
  } else {
    occ.polys.sp <- bind.shp(occ.polys.lst, sp.nm = sp.nm, save=save, crs.set = crs.set) # , o.path = o.path
  }
  return(occ.polys.sp)
}


### for a package see:
# # https://cran.r-project.org/web/packages/clValid/vignettes/clValid.pdf
#
# # for methods
# # nb <- NbClust::NbClust(cbind(spp.occ$LONG, spp.occ$LAT), diss=NULL, distance = "euclidean",
# #                        min.nc=2, max.nc=15, method = "kmeans",
# #                        index = "ball", alphaBeale = 0.1)
# # k <- nb$Best.nc[1]
# # clust <- nb$Best.partition
#
# # elbow.k <- function(mydata){
# #   dist.obj <- dist(mydata)
# #   hclust.obj <- hclust(dist.obj)
# #   css.obj <- GMD::css.hclust(dist.obj, hclust.obj)
# #   elbow.obj <- GMD::elbow.batch(css.obj)
# #   k <- elbow.obj$k
# #   return(k)
# # }
# #
# # elbow.k(cbind(spp.occ.list$Bvarieg$LONG, spp.occ.list$Bvarieg$LAT))
#
#
# d <- cbind(spp.occ.list$Bvarieg$LONG, spp.occ.list$Bvarieg$LAT)
# apclus <- apcluster::apcluster(apcluster::negDistMat(r=2), d)
# plot(apclus, d)
# k <- length(apclus@clusters)
# clust <- apclus@idx
#
# # d <- cbind(spp.occ.list$Bvarieg$LONG, spp.occ.list$Bvarieg$LAT)
# # hclust.obj <- hclust(dist(d))
# # clust <- cutree(hclust.obj, k)
#
#
#
# # function to create N polygons for a species (whenever distribution seems disjoint) and save in a single .shp
#                       (spp.occ, k=2, convex=TRUE, alpha=10, sp.nm = "sp1", crs.set = NULL)
# poly.splt <- function(spp.occ, k=NULL, c.m="AP", convex=TRUE, alpha=10, sp.nm = "sp1", crs.set = NA){
#   # We need to create separated polygons, because we don't want such large area without points.
#
#   if(is.null(k)){
#     # http://www.sthda.com/english/articles/29-cluster-validation-essentials/96-determining-the-optimal-number-of-clusters-3-must-know-methods/
#     # https://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters
#     # for a package for further improvements see:
#     # https://cran.r-project.org/web/packages/clValid/vignettes/clValid.pdf
#     d <- cbind(spp.occ$LONG, spp.occ$LAT)
#     if(c.m = "E"){ ## ELBOW method
#       dist.obj <- dist(d)
#       hclust.obj <- hclust(dist.obj)
#       css.obj <- GMD::css.hclust(dist.obj, hclust.obj)
#       elbow.obj <- GMD::elbow.batch(css.obj)
#       k <- elbow.obj$k
#       clust <- cutree(hclust.obj, k)
#     } else if (c.m = "AP") { # Affinity Propagation (AP)
#       hclust.obj <- hclust(dist(d))
#       apclus <- apcluster::apcluster(apcluster::negDistMat(r=2), d)
#       k <- length(apclus@clusters)
#       clust <- apclus@idx
#       clust <- as.factor(clust)
#       levels(clust) <- 1:k
#       clust <- as.numeric(clust)
#     }
#   } else { # Hierarchical Clustering
#     # https://stackoverflow.com/questions/28672399/spatial-clustering-in-r-simple-example
#     d <- cbind(spp.occ$LONG, spp.occ$LAT)
#     hclust.obj <- hclust(dist(d))
#     clust <- cutree(hclust.obj, k)
#   }
#
#
#   # create one polygon for each set of points
# spp.k.list <- lapply(1:k, function(i){spp.occ[clust==i,]})
# names(spp.k.list) <- paste0(sp.nm, seq_along(spp.k.list))
# occ.polys.lst <- poly.c.batch(spp.k.list, convex=convex, alpha=alpha, plot=FALSE, save=FALSE)
# occ.polys.sp <- bind.shp(occ.polys.lst, sp.nm = sp.nm, crs.set = crs.set) # , o.path = o.path
# # raster::crs(occ.polys.sp) <- crs.set
# # if(!is.null(crs.set)){raster::projection(occ.polys.sp) <- crs.set}
# sp::plot(occ.polys.sp)
# spp.occ <- as.data.frame(spp.occ)
# sp::coordinates(spp.occ) <- ~LONG+LAT
# sp::plot(spp.occ, col="red", add=TRUE)
# return(occ.polys.sp)
# }
#


#' Create occ polygon for several species
#'
#' This function will use a list of coordinates of species occurence data and create minimum convex polygons
#' for each element in the list.
#' It is possible to create concave or convex polygons, create several small polygons based on clusters
#' of points.
#'
#' @param spp.occ.list A named list of species occurence points, either as "data.frame" or "SpatialPoints"/"SpatialPointsDataFrame"
#' @param plot logical. Plot results or not?
#' @param save.pts logical. Save each species' occurence points as shapefile?
#' @param numCores Number of cores to use for parallelization. If set to 1, no paralellization is performed
#' @inheritParams poly.c
#' @inheritParams poly.splt
# #' @inheritParams mxnt.c.batch
#'
#' @seealso \code{\link{poly.c}}, \code{\link{poly.splt}}, \code{\link[NbClust]{NbClust}}
#' #' @return A named list of spatial polygons built using coordinates
#' @examples
#' Bvarieg.occ <- read.table(paste(system.file(package="dismo"),
#'  "/ex/bradypus.csv", sep=""), header=TRUE, sep=",")
#' colnames(Bvarieg.occ) <- c("SPEC", "LONG", "LAT")
#' spp.occ.list <- list(Bvarieg = Bvarieg.occ)
#' occ.polys <- poly.c.batch(spp.occ.list)
#' occ.polys <- poly.c.batch(spp.occ.list, convex=TRUE, alpha=10)
#' @export
poly.c.batch <- function(spp.occ.list, k = 1, c.m = "AP", r = 2, q = .3,
                         distance = "euclidean", min.nc = 2, max.nc = 20,
                         method = "mcquitty", index = "all", alphaBeale = 0.1,
                         convex = T, alpha = 10,
                         plot = T, save = T, save.pts = F, numCores = 1,
                         crs.set = "+proj=longlat +datum=WGS84"){ #, o.path=NULL

  occ.pgns <- vector("list", length(spp.occ.list)) # , names=
  sp.nm <- names(spp.occ.list)
  # sp.nm2 <- paste(sp.nm, "occ.poly", sep = ".")


  o.path.pts <- "1_sppData/occ.pts"
  if(dir.exists("1_sppData")==FALSE) dir.create("1_sppData")
  if(save.pts){
    if(dir.exists(o.path.pts)==FALSE) dir.create(o.path.pts)
  }
  #
  f.poly <- function(i, spp.occ.list, o.path.pts, k,
                     sp.nm, convex, alpha, save, save.pts, crs.set,
                     c.m, r, q, distance, min.nc, max.nc,
                     method, index, alphaBeale, plot){
    occ.spdf <- spp.occ.list[[i]]
    if(!any(class(occ.spdf) %in% c("SpatialPoints", "SpatialPointsDataFrame"))){
      lon.col <- colnames(occ.spdf)[grep("^lon$|^long$|^longitude$", colnames(occ.spdf), ignore.case = T, fixed = F)][1]
      lat.col <- colnames(occ.spdf)[grep("^lat$|^latitude$", colnames(occ.spdf), ignore.case = T)][1]
      sp::coordinates(occ.spdf) <- c(lon.col, lat.col)
    }
    if(save.pts){
      raster::shapefile(occ.spdf, filename = paste(o.path.pts, paste0(paste(names(spp.occ.list)[i], "occ.pts", sep = "."),".shp"), sep = "/" ), overwrite=TRUE)
    }

    if(k == 1){
      resu <- poly.c(occ.spdf, sp.nm=sp.nm[i], convex=convex, alpha=alpha, save=save, crs.set=crs.set) # , o.path=o.path
    } else if (k != 1) {
      resu <-  poly.splt(occ.spdf, k = k, c.m = c.m, r = r, q = q,
                         distance = distance, min.nc = min.nc, max.nc = max.nc,
                         method = method, index = index, alphaBeale = alphaBeale,
                         convex=convex, alpha=alpha, sp.nm=sp.nm[i], save=save, crs.set=crs.set)
    }

    if(plot){
      sp::plot(resu, main=names(spp.occ.list)[i])
      sp::plot(occ.spdf, col="red", add=TRUE)
    }
    return(resu)
  }

  if(numCores>1){
    plot <- F

    cl<-parallel::makeCluster(numCores)

    occ.pgns <- parallel::clusterApply(cl, base::seq_along(spp.occ.list),
                                       function(i, spp.occ.list, o.path.pts, k,
                                                sp.nm, convex, alpha, save, save.pts, crs.set,
                                                c.m, r, q, distance, min.nc, max.nc,
                                                method, index, alphaBeale, plot){

                                         f.poly(i, spp.occ.list, o.path.pts, k,
                                                sp.nm, convex, alpha, save, save.pts, crs.set,
                                                c.m, r, q, distance, min.nc, max.nc,
                                                method, index, alphaBeale, plot)

                                       }, spp.occ.list, o.path.pts, k,
                                       sp.nm, convex, alpha, save, save.pts, crs.set,
                                       c.m, r, q, distance, min.nc, max.nc,
                                       method, index, alphaBeale, plot)

    parallel::stopCluster(cl)

  } else {
    occ.pgns <- lapply(base::seq_along(spp.occ.list),
                       function(i, spp.occ.list, o.path.pts, k,
                                sp.nm, convex, alpha, save, save.pts, crs.set,
                                c.m, r, q, distance, min.nc, max.nc,
                                method, index, alphaBeale, plot){

                         f.poly(i, spp.occ.list, o.path.pts, k,
                                sp.nm, convex, alpha, save, save.pts, crs.set,
                                c.m, r, q, distance, min.nc, max.nc,
                                method, index, alphaBeale, plot)

                       }, spp.occ.list, o.path.pts, k,
                       sp.nm, convex, alpha, save, save.pts, crs.set,
                       c.m, r, q, distance, min.nc, max.nc,
                       method, index, alphaBeale, plot)

  }

  names(occ.pgns) <- names(spp.occ.list)
  return(occ.pgns)
}



#' Bind list of SpatialPolygons into a single SpatialPolygons object
#'
#' This function will bind a list of splitted polygons (from poly.splt()) into a single SpatialPolygons object.
#'
#' @param occ.polys list of SpatialPolygons to bind
#' @inheritParams poly.c
#'
#' @seealso \code{\link{poly.c.batch}}, \code{\link{poly.c}}, \code{\link{poly.splt}}, \code{\link[NbClust]{NbClust}}
#' @return shapefile with binded polygons
#' @keywords internal
#' @export
bind.shp <- function(occ.polys, sp.nm="sp.nm", save=TRUE, crs.set = "+proj=longlat +datum=WGS84"){ # , o.path = "occ.poly"
  o.path <- "1_sppData/occ.poly"
  if(dir.exists("1_sppData")==FALSE) dir.create("1_sppData")
  if(dir.exists(o.path)==FALSE) dir.create(o.path)

  # http://r-sig-geo.2731867.n2.nabble.com/merging-several-shapefiles-into-one-td6401613.html
  # Get polygons and change IDs
  uid <- 1
  poly.l <- vector("list", length(occ.polys))
  for (i in 1:length(occ.polys)) {
    temp.data <- occ.polys[[i]]
    n <- length(methods::slot(temp.data, "polygons"))
    temp.data <- sp::spChFIDs(temp.data, as.character(uid:(uid+n-1)))
    uid <- uid + n
    poly.l[[i]] <- temp.data
  }

  # mapunit polygon: combin remaining  polygons with first polygoan
  p <- do.call(raster::bind, poly.l)

  { ## Convert to "SpatialPolygonsDataFrame"
    # poly.data@data$ID <- seq_along(poly.data@data$ID)
    # Extract polygon ID's
    pid <- sapply(methods::slot(p, "polygons"), function(x) methods::slot(x, "ID"))
    # Create dataframe with correct rownames
    p.df <- data.frame(ID = 1:length(p), row.names = pid) # , sp.nm = sp.nm
    # Coersion
    p <- sp::SpatialPolygonsDataFrame(p, p.df)
  }
  sp.nm <- paste(sp.nm, "occ.poly", sep = ".")
  filename <- paste(o.path, paste0(sp.nm,".shp"), sep = "/" )

  if(save){
    raster::shapefile(p, filename = filename, overwrite = TRUE)
  }
  return(p)
}


#' Create buffer based on species polygons
#'
#' This funcion will create a buffer using species polygons. Buffer width may be manually specified or
#' calculated based on the extent of the SpatialPolygons object (i.e. mean of latitudinal and longitudinal extent).
#' If width is calculated based on the extent of the SpatialPolygons object, it can be adjusted (enlarged or reduced)
#' using 'mult' argument.
#'
#' @param occ.polys list of SpatialPolygons, usually obj returned from poly.c.batch()
#' @param bffr.width Buffer width. See 'width' of ?rgeos::gBuffer
#' @param mult How much expand bffr.width
# #' @param plot Boolean, to draw plots or not
# #' @param quadsegs see ?rgeos::gBuffer
#' @inheritParams rgeos::gBuffer
#' @inheritParams poly.c
#' @inheritParams poly.c.batch
#'
#' @seealso \code{\link[rgeos]{gBuffer}}, \code{\link{poly.c.batch}}
#' @return A named list of SpatialPolygons
#' @examples
#' Bvarieg.occ <- read.table(paste(system.file(package="dismo"),
#'  "/ex/bradypus.csv", sep=""), header=TRUE, sep=",")
#' colnames(Bvarieg.occ) <- c("SPEC", "LONG", "LAT")
#' spp.occ.list <- list(Bvarieg = Bvarieg.occ)
#' occ.polys <- poly.c.batch(spp.occ.list)
#' occ.b <- bffr.batch(occ.polys, bffr.width=1.5)
#' @export
bffr.batch <- function(occ.polys, bffr.width = NULL, mult = .2, quadsegs = 100, numCores = 1, crs.set = NULL, plot = T){ # , o.path = "occ.poly"
  o.path <- "1_sppData/occ.bffr"
  if(dir.exists("1_sppData")==FALSE) dir.create("1_sppData")
  if(dir.exists(o.path)==FALSE) dir.create(o.path)

  # https://gis.stackexchange.com/questions/194848/creating-outside-only-buffer-around-polygon-using-r
  occ.b <- vector("list")
  if(length(mult)==1){ mult <- rep(mult, length(occ.polys))}
  TF.b.w <- is.null(bffr.width)

  f.bffr <- function(i, occ.polys, crs.set, TF.b.w, bffr.width,
                     quadsegs, mult, o.path){
    n.occp.i <- names(occ.polys)[i]
    occ.polys.i <- occ.polys[[i]]
    if(!is.null(crs.set) & is.null(raster::projection(occ.polys.i))){raster::projection(occ.polys.i) <- crs.set}
    if(TF.b.w){
      ext.proj <- raster::extent(occ.polys.i)
      bffr.width <- mean(c(abs(ext.proj[1]-ext.proj[2]), abs(ext.proj[3]-ext.proj[4])))*mult[i]
    }
    cat(c("Buffer width for", n.occp.i, "is", bffr.width, "\n"))
    occ.b.i <- rgeos::gBuffer(occ.polys.i, width=bffr.width, quadsegs=quadsegs)
    raster::shapefile(occ.b.i, filename = paste(o.path, paste0(n.occp.i, ".bffr", ".shp"), sep = "/" ), overwrite=TRUE)
    occ.b.i <- raster::shapefile(paste(o.path, paste0(n.occp.i, ".bffr", ".shp"), sep = "/" ))

    if(plot == T){
      sp::plot(occ.b.i, col="green", main=n.occp.i)
      sp::plot(occ.polys.i, border="red", add=TRUE)
    }
    return(occ.b.i)
  }

  if(numCores>1){
    plot <- F

    cl<-parallel::makeCluster(numCores)

    occ.b <- parallel::clusterApply(cl, base::seq_along(occ.polys),
                                       function(i, occ.polys, crs.set, TF.b.w, bffr.width,
                                                quadsegs, mult, o.path){

                                         f.bffr(i, occ.polys, crs.set, TF.b.w, bffr.width,
                                                quadsegs, mult, o.path)

                                       }, occ.polys, crs.set, TF.b.w, bffr.width,
                                    quadsegs, mult, o.path)

    parallel::stopCluster(cl)

  } else {
    occ.b <- lapply(base::seq_along(occ.polys),
                    function(i, occ.polys, crs.set, TF.b.w, bffr.width,
                             quadsegs, mult, o.path){

                      f.bffr(i, occ.polys, crs.set, TF.b.w, bffr.width,
                             quadsegs, mult, o.path)

                    }, occ.polys, crs.set, TF.b.w, bffr.width,
                    quadsegs, mult, o.path)

  }
  names(occ.b) <- names(occ.polys)
  return(occ.b)
}




#' Cut calibration area based on a list of SpatialPolygons
## #' Crop environmental variables for each species
#' Use a list of SpatialPolygons to crop environmental variables for each species.
#'
#' @param occ.b list of SpatialPolygons, usually returned from "bffr.batch" function
#' @param env.uncut raster brick or stack to be cropped
#' @inheritParams poly.c.batch
#'
#' @seealso \code{\link[raster]{crop}}, \code{\link{bffr.batch}}
#' @return list [one element for each species] of cropped environmental variables. Details in ?raster::crop
#' @examples
#' Bvarieg.occ <- read.table(paste(system.file(package="dismo"),
#'  "/ex/bradypus.csv", sep=""), header=TRUE, sep=",")
#' colnames(Bvarieg.occ) <- c("SPEC", "LONG", "LAT")
#' spp.occ.list <- list(Bvarieg = Bvarieg.occ)
#' occ.polys <- poly.c.batch(spp.occ.list)
#' occ.b <- bffr.batch(occ.polys, bffr.width=1.5)
#' env.uncut <- brick("path/to/env")
#' occ.b.env <- env.cut(occ.b, env.uncut)
#' for(i in 1:length(occ.b.env)){
#'    plot(occ.b.env[[i]][[1]])
#'    plot(occ.b[[i]], add=TRUE)
#' }
#' @export
env.cut <- function(occ.b, env.uncut, numCores = 1){
  path.env.out <- "2_envData/area.calib"
  if(dir.exists("2_envData")==FALSE) dir.create("2_envData")
  if(dir.exists(path.env.out)==FALSE) dir.create(path.env.out)
  ## Clipping rasters for each species
  occ.b.env <- vector("list")

  f.cut <- function(i, occ.b, env.uncut, path.env.out){
    n.occ.b.i <- names(occ.b)[i]
    occ.b.i <- occ.b[[i]]
    cat(c("Cutting environmental variables of species", i, "of", length(occ.b), "\n"))
    env.i <- raster::crop(env.uncut, raster::extent(occ.b.i))
    env.i <- raster::mask(env.i, occ.b.i)
    env.i <- raster::writeRaster(env.i,
                                 filename = paste(path.env.out, paste("envData.", n.occ.b.i, ".grd", sep=''), sep='/'),
                                 format = "raster", overwrite = T)
    return(env.i)
  }

  if(numCores>1){

    cl <- parallel::makeCluster(numCores)

    occ.b.env <- parallel::clusterApply(cl, base::seq_along(occ.b),
                                    function(i, occ.b, env.uncut, path.env.out){

                                      f.cut(i, occ.b, env.uncut, path.env.out)

                                    }, occ.b, env.uncut, path.env.out)

    parallel::stopCluster(cl)

  } else {
    occ.b.env <- lapply(base::seq_along(occ.b),
                        function(i, occ.b, env.uncut, path.env.out){

                          f.cut(i, occ.b, env.uncut, path.env.out)

                        }, occ.b, env.uncut, path.env.out)

  }

  names(occ.b.env) <- names(occ.b)

  return(occ.b.env)
}



#' Spatially thin a list of species occurrence data
#'
#' Will use spThin optimisation algorithm to subset the dataset such that
#' all occurrence locations are a minimum distance apart. This process helps
#' reduce the effect of biases in observation records on the predictive
#' performance of ecological niche models.
#'
#' Make sure coordinates are in decimal degrees. This function will use
#' great.circle.distance to thin the datasets
#'
#' @param loc.data.lst Named list containing data.frames/SpatialPoints/SpatialPointsDataFrame
#' of species occurence locations. Each data.frame can include several columnns, but must
#' include at minimum a column of latitude and a column of longitude values
#' @inheritParams spThin::thin
# #' @inheritParams spThin::spThin
#' @inheritParams poly.c.batch
#'
#' @seealso \code{\link[spThin]{thin}}, \code{\link{loadTocc}}
# #' @seealso \code{\link[spThin]{spThin}}, \code{\link{loadTocc}}
#' @return Named list containing thinned datasets for each species. See ?thin of spThin package.
# #'  Also, by default it saves log file and the first thinned dataset in the folder "occ.thinned.full".
#' @examples
#' thinned.dataset.batch <- thin.batch(loc.data.lst = spp.occ.list)
#' plotThin(thinned.dataset.batch[[1]])
#' length(thinned.dataset.batch[[1]])
#' @export
thin.batch <- function(loc.data.lst = list(),
                       lat.col = NULL, long.col = NULL,
                       spec.col = NULL,
                       thin.par = 10, reps = 10, # reps = 1000 thin.par 'é a distancia min (km) para considerar pontos distintos
                       locs.thinned.list.return = TRUE,
                       write.files = TRUE,
                       max.files = 1,
                       # out.dir = "3_occ.thinned.full",
                       write.log.file = TRUE) {

  thinned_dataset_full <- vector(mode = "list", length = length(loc.data.lst))
  out.dir <- "1_sppData/occ.thinned.full"
  if(dir.exists("1_sppData")==FALSE) dir.create("1_sppData")
  if(dir.exists(out.dir)==FALSE) dir.create(out.dir)

  spp <- names(loc.data.lst)


  t.loc <- function(i, loc.data.lst,  spp, ...){
    occ.spdf <- loc.data.lst[[i]]
    if(any(class(occ.spdf) %in% c("SpatialPoints", "SpatialPointsDataFrame"))){
      occ.spdf <- as.data.frame(occ.spdf)
    }
    if(is.null(lat.col) | is.null(long.col)) {
      long.col <- colnames(occ.spdf)[grep("^lon$|^long$|^longitude$", colnames(occ.spdf), ignore.case = T, fixed = F)][1]
      lat.col <- colnames(occ.spdf)[grep("^lat$|^latitude$", colnames(occ.spdf), ignore.case = T)][1]
    }
    if(is.null(spec.col)){
      spec.col <- colnames(occ.spdf)[grep("^spec$|^species$|^especie$", colnames(occ.spdf), ignore.case = T)][1]
    }


    th.ds <- spThin::thin(occ.spdf, # loc.data.lst[[i]],
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

    return(th.ds)
  }
  thinned_dataset_full <- lapply(1:length(loc.data.lst), t.loc, loc.data.lst=loc.data.lst, spp=spp )
  names(thinned_dataset_full) <- spp
  return(thinned_dataset_full)
}

# thin.batch <- function(loc.data.lst,
#                        dist = 10000, method = "heuristic", nrep = 10,
#                        great.circle.distance = FALSE,
#                        # lat.col = "lat", long.col = "lon", spec.col = "species",
#                        # thin.par = 10, reps = 10, locs.thinned.list.return = TRUE,
#                        # write.files = TRUE, max.files = 1, write.log.file = TRUE,
#                        # # numCores = 1,
#                        ...) {
#
#   out.dir <- "1_sppData/occ.thinned.full"
#   if(dir.exists("1_sppData")==FALSE) dir.create("1_sppData")
#   if(dir.exists(out.dir)==FALSE) dir.create(out.dir)
#
#   spp <- names(loc.data.lst)
#
#   if(utils::packageVersion("spThin")>= 1.0){
#     t.loc <- function(i, loc.data.lst,
#                       dist, method, nrep, great.circle.distance,
#                       spp, out.dir, ...
#                       # spp, lat.col, long.col, spec.col,
#                       # thin.par, reps, locs.thinned.list.return,
#                       # write.files, max.files, out.dir, write.log.file
#     ){
#
#       occ.spdf <- loc.data.lst[[i]]
#       if(!class(occ.spdf) %in% c("SpatialPoints", "SpatialPointsDataFrame")){
#         lon.col <- colnames(occ.spdf)[grep("^lon$|^long$|^longitude$", colnames(occ.spdf), ignore.case = T, fixed = F)][1]
#         lat.col <- colnames(occ.spdf)[grep("^lat$|^latitude$", colnames(occ.spdf), ignore.case = T)][1]
#         sp::coordinates(occ.spdf) <- c(lon.col, lat.col)
#       }
#
#       th.ds <- spThin::spThin(occ.spdf,
#                               dist = dist, method = method, nrep = nrep,
#                               great.circle.distance =  great.circle.distance, ...
#       )
#       wtd <- which.max(sapply(th.ds@samples, length))
#       utils::write.csv(as.data.frame(th.ds[[wtd]]), paste0(out.dir, "/", spp[i], ".occ.thinned.csv"))
#       return(th.ds)
#     }
#
#     thinned.dataset.full <- lapply(base::seq_along(loc.data.lst),
#                                    function(i, loc.data.lst, # spp,
#                                             dist, method, nrep, great.circle.distance,
#                                             spp, out.dir, ...
#                                             # lat.col, long.col, spec.col,
#                                             # thin.par, reps, locs.thinned.list.return,
#                                             # write.files, max.files, out.dir, write.log.file
#                                    ){
#
#                                      t.loc(i, loc.data.lst, # spp,
#                                            dist, method, nrep, great.circle.distance,
#                                            spp, out.dir, ...
#                                            # lat.col, long.col, spec.col,
#                                            # thin.par, reps, locs.thinned.list.return,
#                                            # write.files, max.files, out.dir, write.log.file
#                                      )
#
#                                    }, loc.data.lst, #spp,
#                                    dist, method, nrep, great.circle.distance,
#                                    spp, out.dir, ...
#                                    # lat.col, long.col, spec.col,
#                                    # thin.par, reps, locs.thinned.list.return,
#                                    # write.files, max.files, out.dir, write.log.file
#     )
#
#   } else {
#     t.loc <- function(i, loc.data.lst,
#                       # dist, method, nrep, great.circle.distance,
#                       # spp, out.dir
#                       spp, lat.col, long.col, spec.col,
#                       thin.par, reps, locs.thinned.list.return,
#                       write.files, max.files, out.dir, write.log.file
#     ){
#
#       occ.spdf <- loc.data.lst[[i]]
#       if(class(occ.spdf) %in% c("SpatialPoints", "SpatialPointsDataFrame")){
#         occ.spdf <- as.data.frame(occ.spdf)
#         lon.col <- colnames(occ.spdf)[grep("^lon$|^long$|^longitude$", colnames(occ.spdf), ignore.case = T, fixed = F)][1]
#         lat.col <- colnames(occ.spdf)[grep("^lat$|^latitude$", colnames(occ.spdf), ignore.case = T)][1]
#         # sp::coordinates(occ.spdf) <- c(lon.col, lat.col)
#       }
#
#       spThin::thin(as.data.frame(loc.data.lst[[i]]),
#                    lat.col = lat.col, long.col = long.col,
#                    spec.col = spec.col,
#                    thin.par = thin.par, reps = reps, # reps = 1000 thin.par 'é a distancia min (km) para considerar pontos distintos
#                    locs.thinned.list.return = locs.thinned.list.return,
#                    write.files = write.files,
#                    max.files = max.files,
#                    out.dir = out.dir,
#                    out.base = paste0(spp[i], ".occ.thinned"),
#                    log.file = paste0(out.dir, "/", spp[i], ".occ.thinned.full.log.file.txt"),
#                    write.log.file = write.log.file)
#     }
#
#     thinned.dataset.full <- lapply(base::seq_along(loc.data.lst),
#                                    function(i, loc.data.lst,
#                                             # dist, method, nrep, great.circle.distance,
#                                             # spp, out.dir
#                                             spp, lat.col, long.col, spec.col,
#                                             thin.par, reps, locs.thinned.list.return,
#                                             write.files, max.files, out.dir, write.log.file
#                                    ){
#
#                                      t.loc(i, loc.data.lst,
#                                            # dist, method, nrep, great.circle.distance,
#                                            # spp, out.dir
#                                            spp, lat.col, long.col, spec.col,
#                                            thin.par, reps, locs.thinned.list.return,
#                                            write.files, max.files, out.dir, write.log.file
#                                      )
#
#                                    }, loc.data.lst,
#                                    # dist, method, nrep, great.circle.distance,
#                                    # spp, out.dir
#                                    spp, lat.col, long.col, spec.col,
#                                    thin.par, reps, locs.thinned.list.return,
#                                    write.files, max.files, out.dir, write.log.file
#     )
#   }
#
#   # thinned.dataset.full <- vector(mode = "list", length = length(loc.data.lst))
#   # thinned.dataset.full <- lapply(1:length(loc.data.lst), t.loc, loc.data.lst=loc.data.lst, spp=spp )
#
#
#   names(thinned.dataset.full) <- spp
#   return(thinned.dataset.full)
# }




#' Load "thin.batch" filtered occurrence data
#'
#' @param occ.list.thin named list returned from "thin.batch"
#' @param from.disk boolean. Read from disk or from one of thinned datasets stored in 'occ.list.thin' obj
# #' @param wtd : which thinned dataset?
#'
#' @seealso \code{\link[spThin]{thin}}, \code{\link{thin.batch}}
# #' @seealso \code{\link[spThin]{spThin}}, \code{\link{thin.batch}}
#' @return named list of thinned occurence data for each species in the list
#' @examples
#' occ.locs <- loadTocc(thinned.dataset.batch)
#' @export
loadTocc <- function(occ.list.thin, from.disk=FALSE){ # , wtd=NULL
  occ.l <- vector("list", length(occ.list.thin))
  names(occ.l) <- names(occ.list.thin)

  if (from.disk){ # retrieve from disk
    out.dir <- "1_sppData/occ.thinned.full"
    for(i in 1:length(occ.list.thin)){
      occ.l[[i]] <- utils::read.csv(paste0(out.dir, "/", names(occ.list.thin)[i], ".occ.thinned.csv"),
                                    header=TRUE, sep=',', stringsAsFactors=FALSE)[2:3]
    }
  } else { # retrieve from thinned obj
    # occ.l <- lapply(occ.list.thin, function(x, wtd){
    # x <- x[[wtd]]
    # })
    for(i in 1:length(occ.list.thin)){
      # if(is.null(wtd)){
        # wtd <- which.max(sapply(occ.list.thin[[i]]@samples, length))
        wtd <- which.max(sapply(occ.list.thin[[i]], nrow))
      # }
      # if(wtd > length(occ.list.thin[[i]])) {
      #   stop(paste("There are only", length(occ.list.thin[[i]]), "thinned datasets. 'wtd' was", wtd))
      # }

      occ.l[[i]] <- as.data.frame(sp::coordinates(occ.list.thin[[i]][[wtd]]))
      # occ.l[[i]] <- as.data.frame(occ.list.thin[[i]][[wtd]])
      # occ.l[[i]] <- occ.list.thin[[i]][[wtd]]
    }
  }

  return(occ.l)
}
