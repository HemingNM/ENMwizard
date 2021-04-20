
### ensemble GCMs
#' Group climate scenarios for consensual projections
#'
#' This function groups climate scenarios by supplied groups (e.g. GCMs,
#' RCPs, SSPs) for consensual projections.
#' User need to supply a list containing vectors of names for grouping projections
#' and names of climate scenarios.
#' Projections will be grouped by matching the character vectors in the list
#' against projection names
#'
#' @param groups list containing vectors of names for grouping projections.
#' @param clim.scn.nms Vector with names of climate scenarios
#' @seealso \code{\link{consensus_scn}}, \code{\link{consensus_scn_b}}
#' @return A data.frame with three columns: climate scenario names, group names,
#' and group numbers
#' @examples
#' \dontrun{
#' # see projection names
#' names(mxnt.mdls.preds.cf[[1]]$mxnt.preds)
#' # vector with projection names
#' clim.scn.nms <- c("CCSM4.2050.RCP45",  "MIROC.ESM.2050.RCP45", "MPI.ESM.LR.2050.RCP45",
#'                   "CCSM4.2070.RCP45",  "MIROC.ESM.2070.RCP45", "MPI.ESM.LR.2070.RCP45",
#'                   "CCSM4.2050.RCP85",  "MIROC.ESM.2050.RCP85", "MPI.ESM.LR.2050.RCP85",
#'                   "CCSM4.2070.RCP85",  "MIROC.ESM.2070.RCP85", "MPI.ESM.LR.2070.RCP85")
#' # create two vectors containing grouping codes
#' yr <- c(2050, 2070)
#' rcp <- c("RCP45", "RCP85")
#' # run
#' consensus_gr(groups = list(yr, rcp), clim.scn.nms)
#' }
#' @export
consensus_gr <- function(groups, clim.scn.nms){
  grps <- apply(as.data.frame(expand.grid(groups)), 1,
                function(x){
                  g <-  paste0(x, collapse=".")
                  nbr <- suppressWarnings(as.numeric(substr(g, 1,1)))
                  if(!is.na(nbr)){
                    gsub(paste0("^", nbr,"{1}"), paste0("X", nbr), g)
                  } else {
                    g
                  }
                })

  comb <- apply(as.data.frame(expand.grid(groups)), 1,
                function(x){
                  paste0("(?=.*", paste0(x, collapse=")(?=.*"), ")")
                })
  # create a data frame with groups
  out <- data.frame(clim.scn=clim.scn.nms, consensus.nm=NA)
  # levels(out$consensus.nm) <- grps #unique(c(levels(out$consensus.nm), grps))
  for(i in seq_along(comb)){
    out[grep(comb[i], clim.scn.nms, perl=T), 2] <- grps[i]
  }
  out$consensus.nm <- factor(out$consensus.nm)
  out$consensus.group <- as.numeric(out$consensus.nm)

  return(out[order(out$consensus.group),])
}

#' Create consensual projections of climatic scenarios
#'
#' This function creates consensual projections of climatic scenarios
#' (e.g. GCMs, RCPs, SSPs).
#' User need to supply a 'mcmp' object (returned by \code{\link{proj_mdl}}),
#' and a list containing vectors of names for grouping projections.
#' Projections will be grouped by matching the character vectors in the list
#' against projection names
#'
#' @param ref Character. Name of reference projection (i.e. the one used for
#' calibration and that will not be averaged with any other)
#' @param save Logical. TRUE to save consensual rasters.
#' @inheritParams thrshld
#' @inheritParams consensus_gr
#' @seealso \code{\link{consensus_scn_b}}, \code{\link{consensus_gr}}
#' @return A 'mcmp.l' object. An object returned from function \code{\link{proj_mdl_b}},
#'  containing the consensual projections for each element (species) of the list
# #' @examples
# #' \dontrun{
# #' # see projection names
# #' names(mxnt.mdls.preds.cf[[1]]$mxnt.preds)
# #' # create two vectors containing grouping codes
# #' yr <- c(2050, 2070)
# #' rcp <- c("RCP45", "RCP85")
# #' # run
# #' mxnt.mdls.preds <- consensus_scn(mcmp=mxnt.mdls.preds, groups = list(yr, rcp))
# #' }
#' @export
consensus_scn <- function(mcmp, groups, ref=NULL, sp.nm="species", save=T){
  pred.args <- mcmp$pred.args
  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

  mdl <- names(mcmp$mxnt.preds[[1]])# gsub("Mod.","", )

  grps <- consensus_gr(groups, names(mcmp$mxnt.preds)[!names(mcmp$mxnt.preds) %in% ref])
  comb <- levels(grps$consensus.nm)
  cnssl <- stats::setNames(vector("list", length(comb)), comb)
  cnssl.sd <- cnssl
  for(lr.nm in comb){
    lrs <- names(mcmp$mxnt.preds) %in% grps[grps$consensus.nm == lr.nm,1] # grep(lr.nm, names(mcmp$mxnt.preds))
    cnss.m <- raster::stack()
    cnss.sd <- cnss.m
    for(m in seq_along(mdl)){
      stm <- raster::stack(lapply(mcmp$mxnt.preds[lrs],
                                  function(r, m){
                                    r[[m]] ##
                                  }, m=m))
      cnss.m <- raster::stack(cnss.m, raster::calc(stm, base::mean))
      cnss.sd <- raster::stack(cnss.sd, raster::calc(stm, stats::sd))
    } # for mdl
    names(cnss.m) <- paste0(mdl)
    names(cnss.sd) <- paste0(mdl, "_", "sd")
    if(save){
      dr <- paste0("3_out.MaxEnt/Mdls.", sp.nm, "/", outpt, "/Mdls.consensus")
      if(dir.exists(dr)==F) dir.create(dr, recursive = T)
      cnss.m <- raster::writeRaster(cnss.m,
                                    paste0(dr, "/scn.consensus_", lr.nm, "_mean", ".grd"),
                                    overwrite=T)
      cnss.sd <- raster::writeRaster(cnss.sd,
                                     paste0(dr, "/scn.consensus_", lr.nm, "_sd", ".grd"),
                                     overwrite=T)
    }
    cnssl[[lr.nm]] <- cnss.m
    cnssl.sd[[lr.nm]] <- cnss.sd
  }

  if(!is.null(ref)){
    if(!ref %in% names(mcmp$mxnt.preds)){
      warning("'ref' is not among predictions. Insert correct climatic scenario name.")
      mcmp$scn.consensus <- cnssl
    } else {
      ref.l <- list(mcmp$mxnt.preds[[ref]])
      names(ref.l) <- ref
      mcmp$scn.consensus <- c(ref.l, cnssl)
    }
  } else {
    mcmp$scn.consensus <- cnssl
  }
  mcmp$scn.consensus.sd <- cnssl.sd
  return(mcmp)
}

#' Create consensual projections of climatic scenarios for multiple species
#'
#' This function creates consensual projections of climatic scenarios
#' (e.g. GCMs, RCPs, SSPs) for multiple species.
#' User need to supply a 'mcmp.l' object (returned by \code{\link{proj_mdl_b}}),
#' and a list containing vectors of names for grouping projections.
#' Projections will be grouped by matching the character vectors in the list
#' against projection names
#'
#' @inheritParams consensus_scn
#' @inheritParams thrshld_b
#' @seealso \code{\link{consensus_scn}}, \code{\link{consensus_gr}}
#' @return A 'mcmp.l' object. An object returned from function \code{\link{proj_mdl_b}},
#'  containing the consensual projections for each element (species) of the list
#' @examples
#' \dontrun{
#' # see projection names
#' names(mxnt.mdls.preds.cf[[1]]$mxnt.preds)
#' # create two vectors containing grouping codes
#' yr <- c(2050, 2070)
#' rcp <- c("RCP45", "RCP85")
#' # run
#' mxnt.mdls.preds.cf <- consensus_scn_b(mcmp.l=mxnt.mdls.preds.cf, groups = list(yr, rcp))
#' }
#' @export
consensus_scn_b <- function(mcmp.l, groups, ref=NULL, save=T){
  cnss.l <- lapply(names(mcmp.l),
                   function(i, x, groups, ref, save){
                     consensus_scn(x[[i]], groups, ref, sp.nm=i, save)
                   }
                   , x=mcmp.l, groups=groups, ref, save=save)
  names(cnss.l) <- names(mcmp.l)
  cnss.l
}

# names(mxnt.mdls.preds.cf$Srobustus$mxnt.preds)
# mxnt.mdls.preds.cf2 <- consensus_scn_b(mxnt.mdls.preds.cf, groups = list(yr, rcp), ref="current", save=T)
# names(mxnt.mdls.preds.cf2$Srobustus$scn.consensus)
# plot(mxnt.mdls.preds.cf2$Srobustus$scn.consensus$X2050.RCP45$Mod.EBPM)
# plot(mxnt.mdls.preds.cf2$Srobustus$scn.consensus.sd$X2050.RCP45$Mod.EBPM_sd)
# plot(mxnt.mdls.preds.cf2$Srobustus$scn.consensus$X2050.RCP45$Mod.LowAIC)
# plot(mxnt.mdls.preds.cf2$Srobustus$scn.consensus.sd$X2050.RCP45$Mod.LowAIC_sd)
#
# thrshld_cons <- thrshld_b(mxnt.mdls.preds.cf2)
