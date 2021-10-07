# own functions

####
#' Re-scale a variable to a desired range.
#'
#' @param x vector of values
#' @param x.min minimum values to be considered
#' @param x.max maximu values to be considered
#' @param new.min new minimum values
#' @param new.max new maximum values
#' @export
rescale <- function(x, x.min = NULL, x.max = NULL, new.min = 0, new.max = 1) {
  if(is.null(x.min)) x.min = min(x)
  if(is.null(x.max)) x.max = max(x)
  return(new.min + (x - x.min) * ((new.max - new.min) / (x.max - x.min)))
}

####
#' Thinning of occurrence records using environmental information
#'
#' The function thins spatial points to minimize biases in the environmental space.
#' It divides de range of variables into a desired number of @param bins and selects points
#' closer to the bin center.
#'
#'
#' @return
#' A vector containing the index (row numbers) of selected records.
#'
#' @param data A matrix or data.frame of environmental data at occurrence records.
#' @param bins Number of bins to divide the environmental range of
#' each variable.
#' @export
e_thin_algorithm <- function(data, bins=20){
  ## size of group
  grp_size <- rowSums(apply(data, 2,
                            function(x, bins){
                              class.size <- table(sort(cut(x,
                                                           qunif(seq(0, 1, length.out = bins), min(x), max(x)),
                                                           labels=F, include.lowest=T)))
                              class.size
                              rep(class.size, times=class.size)-1
                            }, bins = bins))

  ## id of group for each variable
  grp_ids <- apply(apply(data, 2,
                         function(x, bins){
                           breaks <- seq(0, 1, length.out = bins)
                           cut(x,
                               stats::qunif(breaks, min(x), max(x)),
                               labels=F, include.lowest=T)
                         }, bins = bins), 1,
                   function(x){
                     paste(x, collapse = "_")
                   })

  ## distance from group center
  grp_centerdist <- rowSums( apply(data, 2,
                                   function(x, bins){
                                     breaks <- seq(0, 1, length.out = bins)
                                     meds <- numeric(length(breaks)-1)
                                     for(i in 1:(length(breaks)-1)){
                                       meds[i] <- mean(breaks[i:(i+1)])
                                     }
                                     xr <- rescale(x)
                                     central.dists <- outer(xr, meds,
                                                            function(x, y){
                                                              sqrt((y-x)^2)
                                                            })
                                     apply(central.dists, 1, min)
                                   }, bins = bins) )


  #### Standard solution: Combine all groups
  grps <- unique(grp_ids)
  sel.rec <- numeric(length(grps))
  for(i in seq_along(grps)){
    gr <- grp_ids == grps[i]
    sel.rec[i] <- which(gr)[which.min(grp_centerdist[gr])]
  }

  return(sort(sel.rec))
}


#' Thinning of occurrence records using environmental information
#'
#' The function thins spatial points to minimize biases in the environmental space.
#' It divides de range of variables into a desired number of @param bins and selects points
#' closer to the bin center.
#'
#' @return
#' A data.frame, matrix, sp or sf object containing the selected records.
#'
#' @param predictors \code{\link[raster]{Raster-class}} object of environmental
#' predictor variables.
#' @param p Two column matrix or data.frame with point coordinates
#' or \code{\link[sp]{SpatialPoints}} of occurrence records.
#' @param bins Number of bins to divide each environmental variable.
#' @param plot Logical. Should results be plotted?
#'
#' @examples \dontrun{
#' file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
#' bradypus <- read.table(file, header=TRUE, sep=',')
#' coord <- bradypus [,2:3]
#'
#' # Download data for present
#' dir.create("rasters")
#' predictors <- getData('worldclim', var='bio', res=10, path="rasters")
#'
#' #cor(data)
#' #vars <- ENMwizard::select_vars(cutoff = 0.9, corr.mat = cor(data), names.only = T)[[1]]
#' #length(vars)
#' env_thin(coord, predictors)
#' }
#'
#' @export
env_thin <- function(p, predictors, bins=20, plot=F){
  data <- raster::extract(predictors, p)
  selRec <- e_thin_algorithm(data, bins=bins)

  if(plot){
    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 0, 0.5))
    plot(data[,1:2], col="gray", pch=19)
    graphics::points(data[selRec, 1:2], col="brown", pch=19)
    plot(p, col="gray", pch=19)
    # maps::map(add = T)
    graphics::points(p[selRec,], col="brown", pch=19)
  }
  return(p[selRec,])
}


#' Thinning of occurrence records using environmental information
#'
#' The function thins spatial points to minimize biases in the environmental space.
#' It divides de range of variables into a desired number of @param bins and selects points
#' closer to the bin center.
#'
#' @return
#' A data.frame, matrix, sp or sf object containing the selected records.
#'
#' @param predictors.lst List of \code{\link[raster]{Raster-class}} object of environmental
#' predictor variables.
#' @param p.lst  List of two column matrix or data.frame with point coordinates
#' or \code{\link[sp]{SpatialPoints}} of occurrence records.
#' @param bins Number of bins to divide each environmental variable.
#' @param plot Logical. Should results be plotted?
#'
#' @examples \dontrun{
#' file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
#' bradypus <- read.table(file, header=TRUE, sep=',')
#' coord <- bradypus [,2:3]
#'
#' # Download data for present
#' dir.create("rasters")
#' predictors <- getData('worldclim', var='bio', res=10, path="rasters")
#'
#' #cor(data)
#' #vars <- ENMwizard::select_vars(cutoff = 0.9, corr.mat = cor(data), names.only = T)[[1]]
#' #length(vars)
#' env_thin_b(list(coord), list(predictors))
#' }
#'
#' @export
env_thin_b <- function(p.lst, predictors.lst, bins=20, plot=F){
  spp <- names(p.lst)

  if(class(predictors.lst) != "list"){
    f_thin <- function(i, p.lst, predictors.lst, ...){
      env_thin(p.lst[[i]], predictors.lst, ...)
    }
  } else {
    f_thin <- function(i, p.lst, predictors.lst, ...){
      env_thin(p.lst[[i]], predictors.lst[[i]], ...)
    }
  }

  thinned_dataset_full <- lapply(1:length(p.lst), f_thin, p.lst=p.lst, predictors.lst=predictors.lst )

  names(thinned_dataset_full) <- spp

  return(thinned_dataset_full)
}
