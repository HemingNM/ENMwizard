#####------- 1. Prepare environmental data

# 1. Prepare enviromental layers
##### 1.1 create occ polygon to crop rasters prior to modelling

f.poly <- function(occ.spdf, o.path = NULL, lr.nm="occ_poly", convex=T, alpha=10, crs.set = NA ){
  if(convex==F){ # convex hulls to crop rasters
    # http://r.789695.n4.nabble.com/Concave-hull-td863710.html#a4688606
    # https://rpubs.com/geospacedman/alphasimple
    library(alphahull)
    library(igraph)
    library(sp)

    ch <- ashape(unique(coordinates(occ.spdf)), alpha=10)
    chg = graph.edgelist(cbind(as.character(ch$edges[, "ind1"]),
                               as.character(ch$edges[, "ind2"])), directed = FALSE)
    if (!is.connected(chg)) {
      stop("Graph not connected")
    } else if (any(degree(chg) != 2)) {
      stop("Graph not circular")
    } else if (clusters(chg)$no > 1) {
      stop("Graph composed of more than one circle")
    }
    # In the next step, we delete one edge of the circle to create a chain,
    # and find the path from one end of the chain to the other.
    cutg <- chg - E(chg)[1]
    # find chain end points
    ends <- names(which(degree(cutg) == 1))
    path <- get.shortest.paths(cutg, ends[1], ends[2])[[1]]
    # this is an index into the points
    pathX <- as.numeric(V(chg)[path[[1]]]$name)
    # join the ends
    pathX <- c(pathX, pathX[1])
    # now show the alpha shape plot with our poly on top
    # plot(ch, lwd = 10, col = "gray")
    # get the points from the ashape object
    # lines(ch$x[pathX, ], lwd = 2)
    coords <- ch$x[pathX, ]
  } else {
    ch <- chull(coordinates(occ.spdf))
    coords <- coordinates(occ.spdf)[c(ch, ch[1]),]
  }
  occ_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)), ID=1)))
  occ_poly <- SpatialPolygonsDataFrame(occ_poly, data=data.frame(ID=1))
  crs(occ_poly) <- crs.set
  if(!is.null(o.path)){
    shapefile(occ_poly, filename = paste(o.path, paste0(lr.nm,".shp"), sep = "/" ), overwrite=TRUE)
  }
  return(occ_poly)
}
# rm(ch, chg, cutg, ends, path, pathX, coords)
# occ_poly <- f.poly(occ.spdf, o.path="occ_poly", lr.nm="occ_poly", convex=T, alpha=10)
# plot(occ_poly)

##### 1.1 create occ polygon for several spp
f.poly.batch <- function(spp.occ.list, o.path=NULL, crs.set = NA, convex=T, alpha=10){
  occ.pgns <- vector("list", length(spp.occ.list)) # , names=
  lr.nm <- paste(names(spp.occ.list), o.path, sep = ".")
  if(!is.null(o.path)) {if(dir.exists(o.path)==F) dir.create(o.path)}
  for(i in 1:length(spp.occ.list)){
    occ.spdf <- data.frame(spp.occ.list[[i]])
    coordinates(occ.spdf) <- ~LONG+LAT
    crs(occ.spdf) <- crs.set
    occ.pgns[[i]] <- f.poly(occ.spdf, o.path=o.path, lr.nm=lr.nm[i], convex=convex, alpha=alpha, crs.set=crs.set)
    plot(occ.pgns[[i]], main=names(spp.occ.list)[i])
    plot(occ.spdf, col="red", add=T)
  }
  names(occ.pgns) <- names(spp.occ.list)
  return(occ.pgns)
}

occ_polys <- f.poly.batch(spp.occ.list, o.path="occ_poly", convex=T, alpha=10, crs.set = crs.set)



#### 1.2 creating buffer
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
    if(!is.null(crs.set) & is.null(projection(occ_polys[[i]]))){projection(occ_polys[[i]]) <- crs.set}
    if(TF.b.w){
      ext_proj <- extent(occ_polys[[i]])
      bffr.width <- mean(c(abs(ext_proj[1]-ext_proj[2]), abs(ext_proj[3]-ext_proj[4])))*mult[i]
    }
    cat(c("Buffer width for", names(occ_b)[i], "is", bffr.width, "\n"))
    occ_b[[i]] <- gBuffer(occ_polys[[i]], width=bffr.width, quadsegs=quadsegs)
    crs(occ_b[[i]]) <- crs.set
    shapefile(occ_b[[i]], filename = paste(o.path, "bffr", paste0(names(occ_b)[i], "_bffr", ".shp"), sep = "/" ), overwrite=TRUE)
    occ_b[[i]] <- shapefile(paste(o.path, "bffr", paste0(names(occ_b)[i], "_bffr", ".shp"), sep = "/" ))

    if(plot == T){
      plot(occ_b[[i]], col="green", main=names(occ_b)[i])
      plot(occ_polys[[i]], border="red", add=T)
      # occ.spdf <- data.frame(spp.occ.list[[i]])
      # coordinates(occ.spdf) <- ~LONG+LAT
      # crs(occ.spdf) <- crs.set
      # plot(occ.spdf, col="red", add=T)
      # plot(NewWorld, add=T)
    }
  }
  return(occ_b)
}

occ_b <- f.bffr(occ_polys, bffr.width=1.5, crs.set=crs.set) #


# path to environmental variables
path.env <- "/Volumes/Samsung SSD/neanderh/Documents/CloudStation/ArcGIS-Virtual Machine/Mapas/Clima e Biosfera/ClimaticScenaries/PR/WorldClim/2_5min/bio_2-5m_bil"
biovars <- paste0("bio", 1:17)#c("bio5", "bio8", "bio10", "bio13", "bio16") # "bio18" tem problemas para o cerrado
pattern.env = 'asc'
path.env.out <- "3_envData"

# 1.3. Cut enviromental layers with M and save in hardrive.
# Get uncut variables
env_uncut <- list.files(path.env, pattern = pattern.env, full.names=TRUE)
env_uncut <- env_uncut[grepl(paste(paste0(biovars, ".", pattern.env), collapse = "|"), env_uncut)]
env_uncut <- stack(env_uncut) #predictors_uncut
crs(env_uncut) <- crs.set

# # Se as variáveis estiverem prontas:
env_uncut <- brick(paste(path.env, "bio.grd", sep="/"))


#### Function to crop environmental variables for each species
f.cut.env <- function(occ_b, env_uncut){
  library(raster)
  path.env.out <- "3_envData"
  ## Clipping rasters for each species
  occ_b_env <- vector("list", length(occ_b))
  names(occ_b_env) <- names(occ_b)
  if(dir.exists(path.env.out)==F) dir.create(path.env.out)

  for(i in 1:length(occ_b)){
    cat(c("Cutting environmental variables of species", i, "of", length(occ_b), "\n"))
    occ_b[[i]] <- spTransform(occ_b[[i]], crs(env_uncut))
    occ_b_env[[i]] <- crop(env_uncut, extent(occ_b[[i]]))
    occ_b_env[[i]] <- mask(occ_b_env[[i]], occ_b[[i]])
    crs(occ_b_env[[i]]) <- crs(env_uncut)
    # if(dir.exists(paste("3_envData", names(spp.occ.list)[i], sep = "/") )==F) dir.create(paste("3_envData", names(spp.occ.list)[i], sep = "/"))
    occ_b_env[[i]] <- writeRaster(occ_b_env[[i]],
                                  filename = paste(path.env.out, paste("envData.", names(occ_b_env)[i], ".grd", sep=''), sep='/'),
                                  format = "raster", overwrite = T)
    # plot(occ_b_env[[i]])
  }

  return(occ_b_env)
}

# 1.3.2 crop environmental variables for each species
occ_b_env <- f.cut.env(occ_b, env_uncut)

for(i in 1:length(occ_b_env)){
  plot(occ_b_env[[i]][[1]])
}
