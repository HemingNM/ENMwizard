##### Community metric functions

#' Group taxa by trait
#'
#' This function will group taxa by trait (trait.df) and return a list of
#' presence/absence (binary distribution maps) stacks for each trait
#'
#' @param trait.df Data frame containing columns with taxon names (taxon),
#' and categorical traits
#' @param trait.cols Character. Name of columns containing traits
#' @param threshold Name of threshold to be found within mtp.l
#' object (returned by \code{\link{thrshld_b}})
#' @param model Name of model to be found within mtp.l
#' object (returned by \code{\link{thrshld_b}})
#' @inheritParams plot_mdl_diff_b
#' @export
comm_gr_trait <- function(mtp.l, trait.df, trait.cols, threshold=NULL, model=NULL){
  # tests before stack layers
  {
    if(is.null(threshold) | !threshold %in% names(mtp.l[[1]][[1]]$binary)) {
      stop(paste("Select one of available thresholds:",
                 paste(names(mtp.l[[1]][[1]]$binary), sep = ", ")))
    }

    if(length(names(mtp.l[[1]][[1]]$binary[[1]]))>1){
      mdls <- sapply(names(mtp.l[[1]][[1]]$binary[[threshold]]), function(x, split){
        strsplit(x, split)[[1]][1]
      }, split= "[.]")

      if(is.null(model) | !model %in% mdls){
        stop(paste("Select one of available models:", mdls))
      } else {
        model <- grep(model, mdls)
      }
    } else {
      model <- 1
    }

    notINmdls <- !trait.df$taxon %in% names(mtp.l)
    notINtbl <- !names(mtp.l) %in% trait.df$taxon
    if(sum(!notINmdls)==0) stop("none of modeled taxa is present in table")
    if(sum(notINmdls)>0) cat(paste("\n", sum(notINmdls), "taxa are not present in modeled list"))
    if(sum(notINtbl)>0) cat(paste("\n", sum(notINtbl), "modeled taxa are not present in table"))
    trait.df <- trait.df[trait.df$taxon %in% names(mtp.l),]
    notINcols <- !trait.cols %in% colnames(trait.df)
    if(sum(notINcols)>0) cat(paste("\n", sum(notINcols), "trait colunms are not present in table:",
                                   paste(colnames(trait.df)[notINcols], sep = ", ")))
    trait.cols <- trait.cols[!notINcols]
  }

  scn.nm <- names(mtp.l[[1]])

  trait.gr <- vector("list", length(trait.cols))
  names(trait.gr) <- trait.cols
  for(tc in trait.cols){
    cat("\n", "Group by", tc)
    traits <- unique(trait.df[,tc])

    trait.comm <- vector("list", length(traits))
    names(trait.comm) <- traits
    for(tr in traits){
      cat("\n", tr)
      FH <- trait.df$taxon[grep(tr, trait.df[,tc])]

      scn.l <- vector("list", length(scn.nm))
      names(scn.l) <- scn.nm
      for(s in scn.nm){
        scn.l[[s]] <- raster::stack(sapply(mtp.l[FH],
                                   function(x){
                                     x[[s]]$binary[[threshold]][[model]]
                                   }))
      }
      trait.comm[[tr]] <- scn.l
    }
    trait.gr[[tc]] <- trait.comm
  }
  return(trait.gr)
}




#' Compute spatial metrics for each trait grouped taxa
#'
#' This function will compute user defined spatial metrics
#' for each trait grouped taxa
#'
#' @param trait.gr Object returned by \code{\link{comm_gr_trait}}
#' @param fun function to be passed to \code{\link[raster]{calc}}
#' @param user.fun User specified function. Can be used to compute focal/window metrics
#' @param result.name Name of directory where results will be saved
#' @inheritParams calib_mdl
#' @inheritParams raster::writeRaster
#' @export
comm_spatial_metric <- function(trait.gr, fun=sum, user.fun=NULL, result.name="richness",
                                numCores=1, format="raster", ...){
  res.dir <- paste0("4.results/", result.name, "/")
  if(dir.exists(res.dir)==F) dir.create(res.dir, recursive = T)

  trait.cols <- names(trait.gr)

  for(tc in trait.cols){
    cat("\n", "Group by", tc)
    traits <- names(trait.gr[[tc]])
    for(tr in traits){
      cat("\n", tr)
      scn.nm <- names(trait.gr[[tc]][[tr]])
      for(s in scn.nm){
        filename <- paste(file.path("4.results", result.name, tc), tr, s, sep = "_")
        if(!is.null(user.fun)){
          trait.gr[[tc]][[tr]][[s]] <- user.fun(trait.gr[[tc]][[tr]][[s]],
                                                  format=format, overwrite=T,
                                                  filename=filename,
                                                  ...)
        } else {
          if(numCores>1){
            cl <- parallel::makeCluster(numCores)
            # parallel::clusterExport(cl)
            trait.gr[[tc]][[tr]][[s]] <- raster::calc(trait.gr[[tc]][[tr]][[s]],
                                                      fun=function(x){ parallel::parApply(cl, x, 1, fun)},
                                                      format=format, overwrite=T,
                                                      filename=filename,
                                                      ... )
            # set flag that cluster is available again
            parallel::stopCluster(cl)
          } else {
            trait.gr[[tc]][[tr]][[s]] <- raster::calc(trait.gr[[tc]][[tr]][[s]],
                                                      fun=fun,
                                                      format=format, overwrite=T,
                                                      filename=filename,
                                                      ...)
          }
        }
      }
    }
  }
  return(trait.gr)
}

