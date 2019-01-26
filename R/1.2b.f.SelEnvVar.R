#' Find, optionally remove, highly correlated variables from a raster brick/stack
#'
#' This function creates a correlation matrix for the layers of a raster brick/stack
#' and returns a brick containing the least correlated variables below the cutoff value.
#' @details
#' This function creates a correlation matrix for the layers of a raster brick/stack using
#' raster::layerStats function.
#' Then, through caret::findCorrelation function, it searches the correlation matrix
#' for pair-wise correlations above the cutoff value and, for each pair of correlated variables,
#' it removes the variable with the largest mean absolute correlation. At the end it returns
#' a raster brick containing only the least correlated variables with correlation below the cutoff
#' value.
#'
#' @param env raster brick/stack
#' @param corr_matrix Correlation matrix from which variables will be selected. If the correlation
#' matrix was already computed from env, you can just input here and choose other cutoff values for
#' selecting variable layers.
#' @param names.only Logical. Return only the names of selected variables (T) or return
#' the raster brick containing the least correlated variables (F).
#' @param plot.dend Logical. Plot dendrogram of correlation distance between variables,
#' showing cutoff limit and selected (black) and discarded (red) variables
#' @param rm.old Logical. Remove (T) old env variables from folder?
#' @inheritParams mxntCalib
#' @inheritParams caret::findCorrelation
#' @inheritParams raster::writeRaster
#' @export
selEnv <- function(env=NULL, corr_matrix=NULL, cutoff=.9, names.only=F, plot.dend=T, rm.old=F, sp.nm="sp", filename=NULL){
  if(is.null(corr_matrix)){
    lStats <- raster::layerStats(env, 'pearson', na.rm=T)
    corr_matrix <- lStats[['pearson correlation coefficient']]
  }
  to.rm <- caret::findCorrelation(corr_matrix, cutoff)
  if(length(to.rm)==0){
    sel.nms <- sort(colnames(corr_matrix))
  } else {
    sel.nms <- sort(colnames(corr_matrix)[-to.rm])
  }

  ### plot dendrogram with selected and discarded variables
  if(plot.dend){
    dist_matrix <- stats::as.dist(1 - base::abs(corr_matrix))
    dend <- stats::as.dendrogram(stats::hclust(dist_matrix)) # as.dendrogram
    ## function to set label color
    labelCol <- function(x, sel.nms) {
      if (stats::is.leaf(x)) {
        ## fetch label
        label <- base::attr(x, "label")
        ## set label color to red for A and B, to blue otherwise
        base::attr(x, "nodePar") <- base::list(lab.col=ifelse(label %in% sel.nms, "black", "firebrick"))
      }
      return(x)
    }

    ## apply labelCol on all nodes of the dendrogram
    dend <- stats::dendrapply(dend, labelCol, sel.nms)
    # graphics::plot(dend, main=sp.nm, ylab = "1 - absolute correlation", xlab = "", sub = "")
    graphics::plot(dend, main=sp.nm, axes=F, ylab = "Absolute correlation", xlab = "", sub = "")
    axis(2, at = seq(0,1,.2), labels=rev(seq(0,1,.2)), ylab = "Absolute correlation")
    graphics::abline(h = 1 - cutoff, col = "red")
  }

  ### return names of selected variables only
  if(names.only){
    return(list(sel_vars=sel.nms, corr_matrix=corr_matrix))
    # return(sel.nms)
  } else { ### return brick with selected variables
    if(is.null(corr_matrix) | is.null(env) |
       nrow(corr_matrix) != raster::nlayers(env) |
       any(!names(env) %in% colnames(corr_matrix)) |
       any(!colnames(corr_matrix) %in% names(env))){
      stop("corr_matrix does not match environmental variables layers")
    }
    path.env.out <- "2_envData/area.calib"
    cat("Selected layers: ", sel.nms)
    env <- env[[-to.rm]]
    env <- raster::writeRaster(env,
                       filename = ifelse(is.null(filename),
                                         paste("2_envData/area.calib", paste0("envDataSel.", sp.nm, ".grd"), sep = "/"),
                                         filename),
                       format = "raster", overwrite=T)
    if(rm.old & is.null(filename)){
      unlink(list.files(path.env.out, pattern = paste0("envData.", sp.nm), full.names=T), recursive = T)
    }
    return(env)
  }
}



#' Find, optionally remove, highly correlated variables from a list (several species) of raster brick/stack
#'
#' This function is a wrapper for ENMwizard::env.sel. See ?ENMwizard::env.sel for details
#' It works with a named list of environmental variables (raster brick/stack).
#' It creates a correlation matrix for the layers of a raster brick/stack for each species
#' (item in the list)
#' and returns a brick containing the least correlated variables below the cutoff value.
#' @param env.l List of raster brick/stack.
#' @param corr_matrix.l List of correlation matrices from which variables will be selected. If the correlation
#' matrix was already computed from env, you can just input here and choose other cutoff values for
#' selecting variable layers.
#' @examples
#' env.sel.b(occ.b.env, .9, names.only=T)
#' occ.b.env <- env.sel.b(occ.b.env, .9, names.only=F, rm.old=F)
#' @inheritParams selEnv
#' @export
selEnvB <- function(env.l, corr_matrix.l=NULL, cutoff=.9, names.only = F, plot.dend = T, rm.old=F, filename=NULL){
  if(is.null(filename)){
    path.env.out <- "2_envData/area.calib"
    if (dir.exists("2_envData") == FALSE) {
      dir.create("2_envData")
    }
    if (dir.exists(path.env.out) == FALSE){
      dir.create(path.env.out)
    # filename <- path.env.out
    }
  }

  env.l.sel <- base::lapply(base::seq_along(env.l), function(i, x, cutoff, corr_matrix.l, names.only, rm.old, filename){ # , n.env.l
    sp.nm <- names(x)[i]
    selEnv(x[[i]], cutoff=cutoff, corr_matrix=corr_matrix.l[[i]]$corr_matrix, names.only = names.only, plot.dend = plot.dend,
            rm.old = rm.old, sp.nm = sp.nm, filename = filename)
  }, x=env.l, cutoff, corr_matrix.l, names.only, rm.old, filename=filename) # , n.env.l

  names(env.l.sel) <- names(env.l) # n.env.l
  return(env.l.sel)
}
