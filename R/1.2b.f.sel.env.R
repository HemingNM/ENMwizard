#' Find highly correlated variables
#'
#' This function creates a correlation matrix for the layers of a raster brick/stack
#' and returns a brick containing the least correlated variables below the cutoff value.
#' @details
#' The function finds variables with correlation above the cutoff and
#' sucessively picks up the variable with the largest number of pair-wise
#' correlations above the cutoff.
#' At each step, the variable is assigned to a group containing all variables
#' with at least one correlated variable in common. The variable
#' with the largest within group mean correlation is discarded.
#' @return
#' A vector of variable names (when names = TRUE) or variable/column indexes
#' (when names = F). If there are no correlations above the cutoff, all variable
#' names (or indexes) is returned
#' @param corr.mat A correlation matrix
#' @param cutoff A numeric value for the pair-wise absolute correlation cutoff.
#' @param names Logical. Should variable names or indexes be returned?
#' @seealso \code{\link{select_vars_b}}, \code{\link{select_vars}}, \code{\link[caret]{findCorrelation}}
#' @export
correlated <- function(corr.mat, cutoff=0.9, names=F){
  corr.mat <- abs(corr.mat)
  xcut <- corr.mat > cutoff
  ncorr <- sort(rowSums(xcut), decreasing = T)
  kept <- names(ncorr[ncorr==1])
  rmv <- character(0)
  ncorr <- ncorr[ncorr!=1]

  for(i in names(ncorr)){
    if(i %in% c(rmv,kept)) next
    col.check <- rownames(corr.mat)[!(rownames(corr.mat) %in% c(rmv,kept))]

    tosum <- xcut[col.check, col.check[xcut[col.check,i]] ]
    if(is.null(nrow(tosum))) next
    group <- rowSums(tosum)>0
    group <- names(group[group])
    if(length(group)>2){
      avg.wtin.corr <- (rowSums(corr.mat[group,group])-1)/(length(group)-1)
      rmv <- c(rmv, names(which.max(avg.wtin.corr)))
    } else {
      avg.wtin.corr <- (rowSums(corr.mat[group,])-1)/(length(group)-1)
      rmv <- c(rmv, names(which.max(avg.wtin.corr)))
    }
  }
  # kept <- sort(c(kept, names(ncorr)[!(names(ncorr) %in% rmv)]))
  index <- which(row.names(corr.mat) %in% rmv)
  if(names){
    return(row.names(corr.mat)[index])
  } else {
    return(index)
  }
}


#' Find, optionally remove, highly correlated variables from a raster brick/stack
#'
#' This function creates a correlation matrix for the layers of a raster brick/stack
#' and returns a brick containing the least correlated variables below the cutoff value.
#' @details
#' This function creates a correlation matrix for the layers of a raster brick/stack,
#' then searches for pair-wise correlations above the cutoff, eliminates variables
#' with most pair-wise correlations and retains the largest number of variables with
#' pair-wise correlation below the cutoff
#' Optionally, a raster brick containing only the variables with correlation below
#' the cutoff is returned.
#'
#' @param env raster brick/stack
#' @param corr.mat Correlation matrix from which variables will be selected. If the correlation
#' matrix was already computed from env, you can just input here and choose other cutoff values for
#' selecting variable layers.
#' @param sample.size Numeric. Number of cells of each layer to be sampled and used for computing
#' correlation. If NULL (default), all cells are used to compute correlation.
#' between variables.
#' @param names.only Logical. Return only the names of selected variables (T) or return
#' the raster brick containing the least correlated variables (F).
#' @param plot.dend Logical. Plot dendrogram of correlation distance between variables,
#' showing cutoff limit and selected (black) and discarded (red) variables
#' @param rm.old Logical. Remove (T) old env variables from folder?
#' @inheritParams calib_mdl
#' @inheritParams caret::findCorrelation
#' @inheritParams raster::writeRaster
#' @seealso \code{\link{select_vars_b}}, \code{\link{correlated}}, \code{\link[caret]{findCorrelation}}
#' @export
select_vars <- function(env = NULL, cutoff = .9, corr.mat = NULL, sample.size = NULL,
                        names.only = F, plot.dend = T, rm.old = F, sp.nm = "sp",
                        filename = NULL){
  if(is.null(corr.mat) & is.null(env)){
    stop("Cannot select variables without environmental variables or correlation matrix")
  }
  if(!is.null(corr.mat) & !is.null(env)){
    if(grepl("raster", class(env), ignore.case = T)){
      n.env <- raster::nlayers(env)
    } else if (class(env) %in% c("data.frame", "matrix")){
      n.env <- nrow(env)
    }
    if(nrow(corr.mat) != n.env |
       any(!names(env) %in% colnames(corr.mat)) |
       any(!colnames(corr.mat) %in% names(env))){
      stop("Variable names in corr.mat does not match names of environmental variables")
    }
  }
  if(is.null(corr.mat)){
    if(!is.null(sample.size)){
      if(grepl("raster", class(env), ignore.case = T)){
        # ns <- round(raster::ncell(env)*sample.prop)
        # message(paste("Computing correlation with", ns, "sampled cells"))
        sample.size <- ifelse(sample.size >= raster::ncell(env), raster::ncell(env), sample.size)

        n <- ceiling(raster::nlayers(env)*(sample.size/raster::ncell(env)))
        if(raster::canProcessInMemory(env, n=n)){
          corr.mat <- stats::cor(raster::sampleRandom(env, sample.size), method = "pearson")
        } else {
          lStats <- raster::layerStats(raster::sampleRandom(env, sample.size, asRaster=T), 'pearson', na.rm=T)
          corr.mat <- lStats[['pearson correlation coefficient']]
        }
      } else if (class(env) %in% c("data.frame", "matrix")){
        sample.size <- ifelse(sample.size >= nrow(env), nrow(env), sample.size)
        sample(1:nrow(env), sample.size)
        corr.mat <- stats::cor(env, method = "pearson")
      }
    } else {
      if(grepl("raster", class(env), ignore.case = T)){
        if(raster::canProcessInMemory(env, n=raster::nlayers(env))){
          corr.mat <- stats::cor(raster::sampleRandom(env, raster::ncell(env)), method = "pearson")
        } else {
          lStats <- raster::layerStats(env, 'pearson', na.rm=T)
          corr.mat <- lStats[['pearson correlation coefficient']]
        }
      } else if (class(env) %in% c("data.frame", "matrix")){
        corr.mat <- stats::cor(env, method = "pearson")
      }
    }
  }
  to.rm <- correlated(corr.mat, cutoff=cutoff)
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
    # standardize corr.mat range between min and 1
    # corr.mat <- (corr.mat-min(corr.mat))/(max(corr.mat)-min(corr.mat))
    corr.mat2 <- corr.mat/base::max(base::abs(corr.mat))
    diag(corr.mat2) <- 1
    dist_matrix <- stats::as.dist(1 - base::abs(corr.mat2))
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
    cat("Selected variables: ", sel.nms, "\n")
    if(grepl("raster", class(env), ignore.case = T)){
      path.env.out <- "2_envData/area.calib"
      env <- env[[-to.rm]]
      if(dir.exists(path.env.out)){
        env <- raster::writeRaster(env,
                                   filename = ifelse(is.null(filename),
                                                     paste("2_envData/area.calib", paste0("envDataSel.", sp.nm, ".grd"), sep = "/"),
                                                     filename),
                                   format = "raster", overwrite=T)
      }
      if(rm.old & is.null(filename)){
        unlink(list.files(path.env.out, pattern = paste0("envData.", sp.nm), full.names=T), recursive = T)
      }
    } else if (class(env) %in% c("data.frame", "matrix")){
      env <- env[-to.rm]
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
#' @seealso \code{\link{select_vars}}, \code{\link{correlated}}, \code{\link[caret]{findCorrelation}}
#' @examples
#'\dontrun{
#' select_vars_b(occ.b.env, .9, names.only=T)
#' occ.b.env <- select_vars_b(occ.b.env, .9, names.only=F, rm.old=F)
#' }
#' @inheritParams select_vars
#' @export
select_vars_b <- function(env.l, cutoff = .9, corr.mat.l = NULL, sample.size = NULL, names.only = F, plot.dend = T, rm.old = F, filename = NULL){
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

  env.l.sel <- base::lapply(base::seq_along(env.l), function(i, x, cutoff, corr.mat.l, sample.size, names.only, rm.old, filename){ # , n.env.l
    sp.nm <- names(x)[i]
    select_vars(x[[i]], cutoff=cutoff, corr.mat=corr.mat.l[[i]]$corr.mat, sample.size=sample.size, names.only = names.only, plot.dend = plot.dend,
            rm.old = rm.old, sp.nm = sp.nm, filename = filename)
  }, x=env.l, cutoff, corr.mat.l, sample.size, names.only, rm.old, filename=filename) # , n.env.l

  names(env.l.sel) <- names(env.l) # n.env.l
  return(env.l.sel)
}
