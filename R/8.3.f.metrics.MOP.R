#' Detecting values out of the environmental range of M
#' @keywords internal
mopout <- function(M, G){
  d1 <- raster::nlayers(M)
  MRange <- raster::cellStats(M, "range")

  out <- which(raster::getValues(sum(raster::stack(
    lapply(1:d1, function(i, G, MRange){
      G[[i]] < MRange[1,i] | G[[i]] > MRange[2,i]
    }, G, MRange))))>0)

  return(sort(out))
}
# mopout(M, G)


# .messi3
#' low level MOP function
#' @keywords internal
mopi3 <- function(x, probs, reff){
  di <- fields::rdist(matrix(x, nrow=1), reff)
  qdi <- stats::quantile(di, probs, na.rm = T)
  di <-  di[di <= qdi]
  return(mean(di))
}



##
#' Extrapolation risk analysis
#'
#' This function will compute the omission rate (OR) for a species' AICc Averaged Model
#' from a 'mcmp' object, based on the selected threshold value.
#'
#' @param M RasterStack of environmental variables from calibration area
#' @param G RasterStack of environmental variables from projection area
#' @param p Percent of values, sampled from calibration area, used as
#' reference to calculate the MOP. Must be >0 and <=1.
#' @param q Quantile. Proportion of closest points in M is to be compared
#' with G to calculate the MOP. Must be >0 and <=1.
#' @param min.M.sz Threshold value to be used to compute OR
#' @inheritParams mxntCalib
#' @inheritParams mxntProj
#' @seealso \code{\link{mop_b}}
#' @return Data frame with average and variance of OR values across partition groups of data
# #' @examples
#' @export
mop <- function(M, G, p=0.1, q=0.1, min.M.sz=100, filename=NULL, scn.nm="", numCores=1, ...) {

  if(any(names(M) != names(G))) {
    stop("M and G must contain the same environmental variables")
  }

  # path to save
  if(is.null(filename)){
    path.mop <- "2_envData/MOP"
    if(dir.exists("2_envData")==FALSE) dir.create("2_envData")
    if(dir.exists(path.mop)==FALSE) dir.create(path.mop)
    filename <- paste0(path.mop, "/", scn.nm, "_MOP.grd")
  }

  # remove NAs from reference (calibration) area
  Mmat <- raster::getValues(M)
  Mmat <- stats::na.omit(Mmat)

  if(nrow(Mmat)<min.M.sz) {
    min.M.sz <- nrow(Mmat)
    # sample.size <- nrow(Mmat)
    # cat("All", nrow(Mmat), "pixel values taken from M \n")
  }
  sample.size <- round(nrow(Mmat)*p)
  if(nrow(Mmat)>sample.size){
    if(sample.size<min.M.sz) {
      sample.size <- ifelse(nrow(Mmat)<min.M.sz, nrow(Mmat), min.M.sz)
    }
    cat(sample.size, "pixels sampled from M \n")
    Mmat <- Mmat[sort(sample(1:nrow(Mmat), sample.size)),]
  } else {
    cat("Used all", nrow(Mmat), "pixels from M \n")
  }

  # find pixels with values outside range of calibration area variables
  outIDs <- mopout(M, G)
  # remove extrapolation pixels from MOP computation
  G[outIDs] <- NA

  # mop.r <- raster(ext=extent(G), res=res(G), crs=crs(G))

  # compute MOP
  if(numCores==1){
    if (raster::canProcessInMemory(G)) {
      mop.r <- raster::raster(G)
      G <- raster::getValues(G)
      G <- apply(G, 1, mopi3, reff=Mmat, probs=q)
      names(mop.r) <- "MOP"
      mop.r <- raster::setValues(mop.r, G)
      # return(mop.r)
    }
    else {
      mop.r <- raster::calc(G,
                    fun = function(x) {
                      apply(x, 1, mopi3, probs=q, reff=Mmat)
                    })
      # fun=function(x, probs=q, reff=Mmat){
      #   di <- fields::rdist(matrix(x, nrow=1), reff)
      #   qdi <- stats::quantile(di, probs, na.rm = T)
      #   ii <-  which(di <= qdi)
      #   return(mean(di[ii]))
      # }) #, filename = filename, overwrite=T, progress='text')
      #
      # mop.r[] <- apply(getValues(G), 1, ff, probs=q, reff=Mmat)
      # plot(mop.r)
    }
  }
  else {
    cl <- parallel::makeCluster(numCores)
    # mop.r[] <- parApply(cl, getValues(G), 1, ff, probs=q, reff=Mmat)
    parallel::clusterExport(cl, varlist=c("q", "Mmat"), envir=environment())

    mop.r <- raster::calc(G,
                  fun=function(x){
                    parallel::parApply(cl, x, 1, mopi3, probs=q, reff=Mmat)
                  }) #, q, Mmat , filename = filename, overwrite=T, progress='text')

    parallel::stopCluster(cl)
  }

  # scale MOP
  max.mop <- raster::cellStats(mop.r, "max")
  mop.r <- 1 - (mop.r/max.mop)

  # remove pixels with values outside range of variables at calibration area
  mop.r[outIDs] <- 0 # set all extrapolation areas to zero
  # set layer name
  names(mop.r) <- "MOP"

  # write raster
  # cat("MOP computed, saving file in: \n", filename, "\n")
  mop.r <- raster::writeRaster(mop.r, filename = filename, overwrite=T, ...)
  return(mop.r)
}


#' Extrapolation risk analysis for a list of species
#'
#' This function will compute the omission rate (OR) for each species' AICc Averaged Model
#' from a 'mcmp.l' object, based on the selected threshold value.
#'
#' @inheritParams mop
#' @inheritParams mxntCalibB
#' @inheritParams mxntProjB
#' @inheritParams plotScnDiff
#' @seealso \code{\link{mop}}
#' @return Data frame with average and variance of OR values across partition groups of data
# #' @examples
#' @export
mop_b <- function(a.calib.l, proj.area.l,
                  p=0.1, q=0.1, min.M.sz=100, ref.scn="current", format = "raster", numCores=1){
    {
      path.res <- "2_envData"
      if (dir.exists(path.res) == FALSE)
        dir.create(path.res)
      path.MOP <- paste0(path.res, "/MOPr")
      if (dir.exists(path.MOP) == FALSE)
        dir.create(path.MOP)
    }

    sp.MOP <- stats::setNames(vector("list", length(a.calib.l)), names(a.calib.l))
    for (sp in names(a.calib.l)) {
      cat(c("\n", "Species: ", sp, "\n"))

      proj.area.spi <- proj.area.l[[sp]]
      a.calib.spi <- a.calib.l[[sp]]

      var.nms <- names(a.calib.spi)
      n.scn <- names(proj.area.spi)
      n.scn <- n.scn[!n.scn %in% ref.scn]

      mop.spi <- vector("list", length(n.scn))
      for (g in n.scn) {
        cat(c(g, " - "))
        G <- proj.area.spi[[g]][[var.nms]]
        mop.spi[[g]] <- mop(M=a.calib.spi[[var.nms]], G, p, q, min.M.sz, filename=NULL, scn.nm=g, numCores)
      }
      mop.spi <- raster::stack(unlist(mop.spi))
      names(mop.spi) <- paste0("MOP_", n.scn)

      flnm <- paste0(path.MOP, "/", sp, "_MOP")
      mop.spi <- raster::writeRaster(mop.spi, flnm,
                             format = format, overwrite=T)
      cat("MOP computed for all climatic scenarios of", sp,". File saved in: \n", flnm, "\n")

      sp.MOP[[sp]] <- mop.spi
    }
    unlink("2_envData/MOP", T)
    return(sp.MOP)
}


