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
  if(sum(is.na(data))>0){ stop("data cannot contain NAs") }
  ## size of group
  grp_size <- rowSums(apply(data, 2,
                            function(x, bins){
                              breaks <- seq(0, 1, length.out = bins)
                              class.size <- table(sort(cut(x,
                                                           stats::qunif(breaks, min(x)-abs(min(x))*.001, max(x)+abs(max(x))*.001),
                                                           labels=F, include.lowest=T)))
                              # class.size
                              rep(class.size, times=class.size)-1
                            }, bins = bins))

  ## id of group for each variable
  grp_ids <- apply(apply(data, 2,
                         function(x, bins){
                           breaks <- seq(0, 1, length.out = bins)
                           cut(x,
                               stats::qunif(breaks, min(x)-abs(min(x))*.001, max(x)+abs(max(x))*.001),
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
  sel.vec <- vector("logical",nrow(data))
  sel.vec[sel.rec] <- T
  return(sel.vec)
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
#' @param lat.col,long.col Name of columns that contain coordinates (latitude and longitude)
#' @param bins Number of bins to divide each environmental variable.
#' @param file either a character string naming a file or a connection open for writing.
#' @param plot Logical. Should results be plotted?
#' @param verbose Logical. Should information be printed on screen?
#'
#' @seealso \code{\link{env_thin_b}}, \code{\link{load_env_thin_occ}}
#'
#' @examples \dontrun{
#' file <- paste(system.file(package="dismo"), "/ex/bradypus.csv", sep="")
#' bradypus <- read.table(file, header=TRUE, sep=',')
#' coord <- bradypus [,2:3]
#'
#' # Download data for present
#' dir.create("rasters")
#' library(raster)
#' predictors <- getData('worldclim', var='bio', res=10, path="rasters")
#'
#' #cor(data)
#' #vars <- ENMwizard::select_vars(cutoff = 0.9, corr.mat = cor(data), names.only = T)[[1]]
#' #length(vars)
#' env_thin(coord, predictors)
#' }
#'
#' @export
env_thin <- function(p, predictors, long.col=NULL, lat.col=NULL, bins=20, file=NULL, plot=F, verbose=F){
  if(inherits(p, c("data.frame", "matrix"))){
    if(inherits(p, "matrix")){
      p <- as.data.frame(p)
    }
    if(is.null(lat.col) | is.null(long.col)) {
      long.col <- colnames(p)[grep("^lon$|^long$|^longitude$", colnames(p), ignore.case = T, fixed = F)][1]
      lat.col <- colnames(p)[grep("^lat$|^latitude$", colnames(p), ignore.case = T)][1]
    }
    coords <- p[,c(long.col, lat.col)]
  } else { # if(class(p) %in% c("SpatialPoints", "SpatialPointsDataFrame")){
    coords <- p
  }

  data <- raster::extract(predictors, coords)
  nas <- rowSums(is.na(data))>0
  if(sum(nas)>0) warning(paste(sum(nas), "record(s) with NA for at least one predictor"))
  selRec <- e_thin_algorithm(data[!nas,], bins=bins)

  if(plot){
    graphics::par(mfrow = c(1, 2), mar = c(7.2, 4, .5, .5))
    plot(data[,1:2], col="white", pch=19)
    graphics::points(data[!nas,][selRec, 1:2], col=grDevices::rgb(34,139,34, 100, maxColorValue=255) , pch=19)
    graphics::points(data[!nas,][!selRec,1:2], col=grDevices::rgb(255,0,0, 130, maxColorValue=255), pch=19, cex=.6)
    # graphics::points(data[nas,1:2], col=grDevices::rgb(255,0,0, 130, maxColorValue=255) , pch=19)
    plot(coords, col="white", pch=19)
    # maps::map(add = T)
    graphics::points(coords[!nas,][selRec,], col=grDevices::rgb(34,139,34, 100, maxColorValue=255) , pch=19)
    graphics::points(coords[!nas,][!selRec,], col=grDevices::rgb(255,0,0, 130, maxColorValue=255), pch=19, cex=.6)
    graphics::points(coords[nas,], col="black", pch=19, cex=.6)
    graphics::par(mfrow = c(1, 1), mar = c(1, .5, .5, .5), xpd=TRUE)
    graphics::legend("bottom", inset=c(.20,-.90), horiz = F, ncol=2,
           title = "Occurrences",
           legend=c("retained", "with NAs", "removed"), pch=19,
           col=c(grDevices::rgb(34,139,34, 150, maxColorValue=255),
                 "black",
                 grDevices::rgb(255,0,0, 150, maxColorValue=255)),
           cex=c(.6,.6,.6))
  }

  m <- paste("\n**********************************************\n
        Environmental thinning\n
        Environmental variables used for thinning:", paste(colnames(data), collapse=", "), "\n
        Number of bins:", bins, "\n
        Original number of records:", nrow(p), "\n
        Number of retained records:", sum(selRec), "\n
        Number of records with NA for at least one variable:", sum(nas), "\n
        Thinning done on:",  as.character(Sys.time()), "\n
        File saved as:", file, "\n
        **********************************************")
  if(verbose){
    cat(m, "\n")
  }
  if(!is.null(file)){
    p$selRec <- F
    p$selRec[nas] <- NA
    p$selRec[!nas][selRec] <- T
    utils::write.csv(as.data.frame(p), file = file)
    sink(gsub(".csv", ".log.txt", file))
    cat(m)
    sink()
  }
  return(coords[!nas,][selRec,])
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
#' @param predictors.lst A \code{\link[raster]{Raster-class}} or list of \code{\link[raster]{Raster-class}}
#' object containing the predictor variables.
#' @param p.lst  List of two column matrix or data.frame with point coordinates
#' or \code{\link[sp]{SpatialPoints}} of occurrence records.
#' @inheritParams env_thin
# #' @param bins Number of bins to divide each environmental variable.
# #' @param plot Logical. Should results be plotted?
#'
#' @seealso \code{\link{env_thin}}, \code{\link{load_env_thin_occ}}
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
env_thin_b <- function(p.lst, predictors.lst, long.col=NULL, lat.col=NULL, bins=20, plot=F, verbose=F){
  spp <- names(p.lst)

  if(!inherits(predictors.lst, "list")){
    f_thin <- function(i, p.lst, predictors.lst, long.col, lat.col, bins, spp, ...){
      file=paste0("1_sppData/occ.thinned.full/",
                  spp[i],"_occ_env_thinned.csv")
      env_thin(p.lst[[i]], predictors.lst, long.col, lat.col, bins, file, ...)
    }
  } else {
    f_thin <- function(i, p.lst, predictors.lst, long.col, lat.col, bins, spp, ...){
      file=paste0("1_sppData/occ.thinned.full/",
                  spp[i],"_occ_env_thinned.csv")
      env_thin(p.lst[[i]], predictors.lst[[i]], long.col, lat.col, bins, file, ...)
    }
  }
  if(!dir.exists("1_sppData/occ.thinned.full"))dir.create("1_sppData/occ.thinned.full", recursive = T)
  thinned_dataset_full <- lapply(1:length(p.lst), f_thin,
                                 p.lst=p.lst, predictors.lst=predictors.lst,
                                 long.col=long.col, lat.col=lat.col, bins=bins,
                                 spp=spp, plot=plot, verbose=verbose)

  names(thinned_dataset_full) <- spp

  return(thinned_dataset_full)
}


#' Load filtered occurrence data
#'
#' Load filtered occurrence data from object returned by \code{\link{env_thin_b}}
#'
#' @param p.lst named list returned from \code{\link{env_thin_b}}
#' @inheritParams env_thin
#'
#' @seealso \code{\link{env_thin}}, \code{\link{env_thin_b}}
#'
#' @return named list of thinned occurence data for each species in the list
#' @examples
#'\dontrun{
#' occ.locs <- load_env_thin_occ(thinned.dataset.batch)
#' }
#' @export
load_env_thin_occ <- function(p.lst, long.col=NULL, lat.col=NULL){
  lapply(p.lst, function(p, ...){
    if(inherits(p, c("data.frame", "matrix"))){
      if(is.null(lat.col) | is.null(long.col)) {
        long.col <- colnames(p)[grep("^lon$|^long$|^longitude$", colnames(p), ignore.case = T, fixed = F)][1]
        lat.col <- colnames(p)[grep("^lat$|^latitude$", colnames(p), ignore.case = T)][1]
      }
      p <- p[,c(long.col, lat.col)]
    } else {
      p <- sp::coordinates(p)
    }
    return(p)
  }, long.col=NULL, lat.col=NULL)
}
