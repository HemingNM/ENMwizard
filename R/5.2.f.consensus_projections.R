
### ensemble GCMs

#' Create consensual projections of climatic scenarios
#'
#' This function creates consensual projections of climatic scenarios
#' (e.g. GCMs, RCPs, SSPs).
#' User need to supply a 'mcmp' object (returned by \code{\link{proj_mdl}}),
#' and a list containing vectors of names for grouping projections.
#' Projections will be grouped by matching the character vectors in the list
#' against projection names
#'
#' @param ref Object returned by \code{\link{proj_mdl_b}}, containing a list of calibrated models
#' and model projections for each species.
#' @inheritParams thrshld
#' @seealso \code{\link{consensus_scn_b}}
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
consensus_scn <- function(mcm, groups, ref=NULL, sp.nm="species", save=T){
  pred.args <- mcm$pred.args
  outpt <- ifelse(grep('cloglog', pred.args)==1, 'cloglog',
                  ifelse(grep("logistic", pred.args)==1, 'logistic',
                         ifelse(grep("raw", pred.args)==1, 'raw', "cumulative")))

  mdl <- names(mcm$mxnt.preds[[1]])# gsub("Mod.","", )
  comb <- apply(as.data.frame(expand.grid(groups)), 1,
                function(x){
                  paste0("((.*", paste0(x, collapse=")(.*"), "))")
                })
  cnssl <- list()
  cnssl.sd <- cnssl
  for(lr.nm in comb){
    c.nm <- gsub("^\\.", "",
                 gsub("[()]", "",
                      gsub(")(", ".",
                           gsub(paste0(ref,"|\\||[.*]"), "", lr.nm), fixed=T))) # c.nm <- gsub("[&]|[*]|[.]|", "", gsub(paste0(ref,"|\\||[.*]"), "", lr.nm))
    # insert "X" before numbers
    nbr <- suppressWarnings(as.numeric(substr(c.nm, 1,1)))
    if(!is.na(nbr)){
      c.nmx <- gsub(paste0("^", nbr,"{1}"), paste0("X", nbr), c.nm)
    } else {
      c.nmx <- c.nm
    }

    lrs <- grep(lr.nm, names(mcm$mxnt.preds))
    cnss.m <- raster::stack()
    cnss.sd <- cnss.m
    for(m in seq_along(mdl)){
      stm <- raster::stack(lapply(mcm$mxnt.preds[lrs],
                          function(r, m){
                            r[[m]] ##
                          }, m=m))
      cnss.m <- raster::stack(cnss.m, base::mean(stm))
      cnss.sd <- raster::stack(cnss.sd, raster::calc(stm, stats::sd))
    } # for mdl
      names(cnss.m) <- paste0(mdl)
      names(cnss.sd) <- paste0(mdl, "_", "sd")
    if(save){
      dr <- paste0("3_out.MaxEnt/Mdls.", sp.nm, "/", outpt, "/Mdls.consensus")
      if(dir.exists(dr)==F) dir.create(dr, recursive = T)
      cnss.m <- raster::writeRaster(cnss.m,
                                 paste0(dr, "/scn.consensus_", c.nm, "_mean", ".grd"),
                                 overwrite=T)
      cnss.sd <- raster::writeRaster(cnss.sd,
                            paste0(dr, "/scn.consensus_", c.nm, "_sd", ".grd"),
                            overwrite=T)
      }
    cnssl[[c.nmx]] <- cnss.m
    cnssl.sd[[c.nmx]] <- cnss.sd
  }

  if(!is.null(ref)){
    if(!ref %in% names(mcm$mxnt.preds)){
      warning("'ref' is not among predictions. Insert correct climatic scenario name.")
      mcm$scn.consensus <- cnssl
    } else {
      ref.l <- list(mcm$mxnt.preds[[ref]])
      names(ref.l) <- ref
      mcm$scn.consensus <- c(ref.l, cnssl)
    }
  } else {
    mcm$scn.consensus <- cnssl
  }
  mcm$scn.consensus.sd <- cnssl.sd
  return(mcm)
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
#' @seealso \code{\link{consensus_scn}}
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
