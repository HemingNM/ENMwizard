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
#' @param corr.mat Correlation matrix from which variables will be selected. If the correlation
#' matrix was already computed from env, you can just input here and choose other cutoff values for
#' selecting variable layers.
#' @param sample.prop Numeric. Proportion of cells of each layer to be sampled and used for computing
#' correlation. Values must be > 0 and <= 1. If 1, all cells are used to compute correlation
#' between variables.
#' @param names.only Logical. Return only the names of selected variables (T) or return
#' the raster brick containing the least correlated variables (F).
#' @param plot.dend Logical. Plot dendrogram of correlation distance between variables,
#' showing cutoff limit and selected (black) and discarded (red) variables
#' @param rm.old Logical. Remove (T) old env variables from folder?
#' @inheritParams calib_mdl
#' @inheritParams caret::findCorrelation
#' @inheritParams raster::writeRaster
#' @seealso \code{\link{select_vars_b}}, \code{\link[caret]{findCorrelation}}
#' @export
select_vars <- function(env = NULL, cutoff = .9, corr.mat = NULL, sample.prop = 0.1,
                        names.only = F, plot.dend = T, rm.old = F, sp.nm = "sp",
                        filename = NULL){
  if(is.null(corr.mat)){
    if(sample.prop<1){
      # corr.mat <- cor(raster::sampleRandom(env, round(ncell(env)*sample.prop)), method = "pearson")
      envs <- raster::sampleRandom(env, round(ncell(env)*sample.prop), asRaster=T)
      lStatss <- raster::layerStats(envs, 'pearson', na.rm=T)
      corr.mat <- lStats[['pearson correlation coefficient']]
    } else {
      lStats <- raster::layerStats(env, 'pearson', na.rm=T)
      corr.mat <- lStats[['pearson correlation coefficient']]
    }
  }
  to.rm <- caret::findCorrelation(corr.mat, cutoff=cutoff)
  if(length(to.rm)==0){
    sel.nms <- sort(colnames(corr.mat))
  } else {
    sel.nms <- sort(colnames(corr.mat)[-to.rm])
  }

  ### plot dendrogram with selected and discarded variables
  if(plot.dend){
	  # modified from rafalib::myplclust
	  myplclust <- function (hclust, labels = hclust$labels, lab.col = rep(1, length(hclust$labels)),
	                         lab.face = rep(1, length(hclust$labels)),
	                         hang = 0.1, xlab = NA, sub = NA, axes=F, ...) {
	    y <- rep(hclust$height, 2)
	    x <- as.numeric(hclust$merge)
	    y <- y[which(x < 0)]
	    x <- x[which(x < 0)]
	    x <- abs(x)
	    y <- y[order(x)]
	    x <- x[order(x)]
	    graphics::plot(hclust, labels = FALSE, hang = hang, xlab = xlab, sub = sub, axes = axes, ...)
	    graphics::text(x = x, y =
	                     if(hang > 0){
	                       (y[hclust$order] - .04 - max(hclust$height) * hang)
	                     } else if(hang == 0){
	                       y[hclust$order] - .04  # y[hclust$order] - .1
	                     } else {
	                       - .04 #(mean(hclust$height) * hang)
	                     },
	                   labels = labels[hclust$order], col = lab.col[hclust$order],
	                   font = lab.face[hclust$order],
	                   srt = 90, adj = c(1, 0.5), xpd = NA, ...)
	  }
	  dist_matrix <- stats::as.dist(1 - base::abs(corr.mat))
	  hcd <- stats::hclust(dist_matrix)
	  lab.col <- ifelse(rownames(corr.mat) %in% sel.nms, "black", "gray30")
	  lab.face <- ifelse(rownames(corr.mat) %in% sel.nms, 2, 1)
	  myplclust(hcd, hang=-.1, axes=F, xlab = "Variables",
	            lab.col = lab.col, lab.face=lab.face, ylab = "Absolute correlation")

	  # # dist_matrix <- stats::dist(corr.mat)
	  # dist_matrix <- stats::as.dist(1 - base::abs(corr.mat))
	  # dend <- stats::as.dendrogram(stats::hclust(dist_matrix)) # as.dendrogram
	  # ## function to set label color
	  # labelCol <- function(x, sel.nms) {
	  #  if (stats::is.leaf(x)) {
	  #    ## fetch label
	  #    label <- base::attr(x, "label")
	  #    ## set label color to red for A and B, to blue otherwise
	  #    base::attr(x, "nodePar") <- base::list(lab.col=ifelse(label %in% sel.nms, "black", "firebrick"))
	  #  }
	  #  return(x)
	  # }
	  #
	  # ## apply labelCol on all nodes of the dendrogram
	  # dend <- stats::dendrapply(dend, labelCol, sel.nms)
	  # # graphics::plot(dend, main=sp.nm, ylab = "1 - absolute correlation", xlab = "", sub = "")
	  # graphics::plot(dend, main=sp.nm, axes=F, ylab = "Absolute correlation", xlab = "", sub = "")
	  graphics::abline(h = 1 - cutoff, col = "firebrick", lwd=1.5)
	  graphics::axis(2, at = seq(0,1,.2), labels=rev(seq(0,1,.2)), ylab = "Absolute correlation")
	  graphics::legend("topright", horiz=F, # title="Variables:",
	                   legend=c("selected vars","removed vars", "cutoff"), text.col=c("black", "gray30", "firebrick"), text.font=c(2, 1,1),
	                   col=c(NA, NA, "firebrick"), lty=c(0,0,1), lwd=1.5, seg.len = 1,
	                   xpd=T, cex=.7)
  }

  ### return names of selected variables only
  if(names.only){
    return(list(sel_vars=sel.nms, corr.mat=corr.mat))
    # return(sel.nms)
  } else { ### return brick with selected variables
    if(is.null(corr.mat) | is.null(env) |
       nrow(corr.mat) != raster::nlayers(env) |
       any(!names(env) %in% colnames(corr.mat)) |
       any(!colnames(corr.mat) %in% names(env))){
      stop("corr.mat does not match environmental variables layers")
    }
    path.env.out <- "2_envData/area.calib"
    cat("Selected variables: ", sel.nms, "\n")
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
#' @param corr.mat.l List of correlation matrices from which variables will be selected. If the correlation
#' matrix was already computed from env, you can just input here and choose other cutoff values for
#' selecting variable layers.
#' @seealso \code{\link{select_vars}}, \code{\link[caret]{findCorrelation}}
#' @examples
#'\dontrun{
#' select_vars_b(occ.b.env, .9, names.only=T)
#' occ.b.env <- select_vars_b(occ.b.env, .9, names.only=F, rm.old=F)
#' }
#' @inheritParams select_vars
#' @export
select_vars_b <- function(env.l, cutoff = .9, corr.mat.l = NULL, sample.prop = 0.1, names.only = F, plot.dend = T, rm.old = F, filename = NULL){
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

  env.l.sel <- base::lapply(base::seq_along(env.l), function(i, x, cutoff, corr.mat.l, sample.prop, names.only, rm.old, filename){ # , n.env.l
    sp.nm <- names(x)[i]
    select_vars(x[[i]], cutoff=cutoff, corr.mat=corr.mat.l[[i]]$corr.mat, sample.prop=sample.prop, names.only = names.only, plot.dend = plot.dend,
            rm.old = rm.old, sp.nm = sp.nm, filename = filename)
  }, x=env.l, cutoff, corr.mat.l, names.only, rm.old, filename=filename) # , n.env.l

  names(env.l.sel) <- names(env.l) # n.env.l
  return(env.l.sel)
}
