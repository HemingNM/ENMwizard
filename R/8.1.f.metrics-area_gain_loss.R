####  Geographical shifts in suitable areas

#' Map shifts in suitable areas
#'
#' Map shifts in suitable areas between a selected climatic scenario
#' and all other climatic scenarios.
#'
#' @inheritParams mop_b
#' @inheritParams thrshld_b
#' @inheritParams plot_mdl_diff
#' @inheritParams calib_mdl
# #' @param ref.scn Selected climatic scenario to compare with all others. Usually "ncurrent".
#' @seealso \code{\link{range_shift_b}}
#' @return A list containing rasters of differences between predictions of climatic scenarios
#' @examples
#' \dontrun{
#' range_shift(mcmp.l=mxnt.mdls.preds.lst, mtp.l=mods.thrshld.lst)
#' }
#' @export
range_shift <- function(mtp, ref.scn="ncurrent",
                        sp.nm="species", numCores=1, format="raster"){
  if(length(mtp)==1){
    stop("Only one climatic scenario. Nothing to compare to.")
  }

  if(dir.exists("4.results/range_shift")==F) dir.create("4.results/range_shift", recursive = T)
  comb.plots <- utils::combn(length(mtp), 2)
  cli.scn.pres <- which(names(mtp) == ref.scn)
  sel.col <- apply(comb.plots == cli.scn.pres, 2, sum)==TRUE # rowsum(comb.plots == cli.scn.pres)==TRUE # comb.plots[, ]
  comb.plots <- matrix(comb.plots[, sel.col], nrow = 2)
  comb.plots[, comb.plots[2,] == cli.scn.pres] <- comb.plots[c(2,1), comb.plots[2,] == cli.scn.pres]

  thrshld_crit <- names(mtp[[comb.plots[1,1]]]$binary)
  m.nms <- names(mtp[[comb.plots[1,1]]]$binary[[1]])
  mSel <- gsub("_$", "", gsub(paste0(".", ref.scn, ".",thrshld_crit, collapse = "|"), "_", m.nms))

  n.scn <- ncol(comb.plots)
  scn.nms <- names(mtp)

  res.thr <- stats::setNames(vector("list", length(thrshld_crit)), thrshld_crit)
  res.scn <- sapply(mSel, function(i, x){x}, x=res.thr, simplify = F)

  for(m in seq_along(mSel)){ # model selection criteria

    for(tc in thrshld_crit){ # threshold criteria
      res.scn[[m]][[tc]] <-
        raster::stack(lapply(1:n.scn, # climatic scenario
                             function(j, mtp, comb.plots, tc, m, scn.nms){
                               stats::setNames(raster::overlay(mtp[[comb.plots[1,j]]]$binary[[tc]][[m]],
                                                        mtp[[comb.plots[2,j]]]$binary[[tc]][[m]],
                                                        fun=function(r1, r2) {
                                                          ifelse((r2+r1)==0, NA, r2-r1)
                                                        }),
                                        paste(scn.nms[comb.plots[1,j]], scn.nms[comb.plots[2,j]], sep = "_to_"))
                             }, mtp, comb.plots, tc, m, scn.nms))

      res.scn[[m]][[tc]] <- raster::writeRaster(res.scn[[m]][[tc]],
                                        paste(paste0("4.results/range_shift/", sp.nm), mSel[m], tc, ref.scn, sep = "_"),
                                        format=format, overwrite=T)

    } # threshold criteria
  } # model selection criteria
  return(res.scn)
}


#' Map shifts in suitable areas
#'
#' Map shifts in suitable areas between a selected climatic scenario
#' and all other climatic scenarios for each species.
#'
#' @inheritParams thrshld_b
#' @inheritParams range_shift
#' @inheritParams plot_mdl_diff_b
#' @inheritParams calib_mdl
#' @seealso \code{\link{range_shift}}
#' @return Return a nested list, containig rasters with differences between
#' a selected climatic scenario and all other climatic scenarios.
#' @examples
#' \dontrun{
#' spp_diff <- range_shift_b(mods.thrshld.lst, ref.scn = "ncurrent")
#'
#' col <- c("red", "gray", "blue") #RColorBrewer::brewer.pal(9,"YlOrRd")
#' breaks <- round(seq(from=-1, to=1, .666), 2)
#' colors <- colorRampPalette(col)(length(breaks)-1)
#'
#' plot(spp_diff[[1]][[1]][[1]], col=colors, breaks=breaks)
#' }
#' @export
range_shift_b <- function(mtp.l, ref.scn="ncurrent", numCores=1, format="raster"){
  sp.nm.l <- names(mtp.l)
  res.list <- stats::setNames(vector("list", length(mtp.l)), sp.nm.l)
  for(sp in 1:length(mtp.l)){ # species
    mtp <- mtp.l[[sp]]
    sp.nm <- sp.nm.l[sp]

    res.list[[sp]] <- range_shift(mtp, ref.scn, sp.nm, numCores, format)
  } # fecha # species

  return(res.list)
}



#' Reshape array to data.frame
#' @param a array
#' @param clim.scen Climate scenario
#' @param threshold threshold
#' @param model model
#' @param location location
#' @keywords internal
array2df <- function(a, clim.scen, threshold, model, location){
  # https://stackoverflow.com/questions/40921426/converting-array-to-matrix-in-r
  data.frame(expand.grid(Clim.scen = clim.scen, # model criteria
                         Threshold = threshold, # threshold criteria
                         Model = model), # scenario diff
             Location = location, #rep((masks), length(ar.mods.t.p)/rep.mdl),
             TotSuitArea = array(aperm(a, c(1,3,2))))

}

##
#' Compute changes in total suitable area
#'
#' Compute changes in total suitable area at multiple climatic scenarios, threshold and model criteria.
#'
#' @param mrs List containing stack or brick of range shifts. See \code{\link{range_shift}}
#' @inheritParams thrshld
#' @inheritParams get_tsa
#' @inheritParams base::round
#' @seealso \code{\link[raster]{area}}, \code{\link{get_tsa_b}}, \code{\link{get_cont_permimport}}, \code{\link{get_fpa}},
#' \code{\link{get_cont_permimport_b}}, \code{\link{get_fpa_b}}
#' @return An object of class "data.table" "data.frame" containing calculated areas of
#'  range shift between a selected climatic scenario and all other climatic scenarios
#'  for threshold, model and masked locations (if projections were masked, see \code{\link{mask_thr_projs_mscn_b}}).
#' @examples
#' \dontrun{
#' range_diff <- range_shift(mods.thrshld.lst$haddadus_binotatus)
#' spp_diff_a <- get_rsa(mrs=range_diff)
#' }
#' @export
get_rsa <- function(mrs, area.raster=NULL, digits=2){ # species, areas
  thrshld.nms <- paste(paste0(".", tnm), collapse = "|")
  c.nms <- gsub(paste0("Mod\\.|", gsub("\\.", "\\\\.", thrshld.nms)), "", names(mrs))
  c.nms2 <- vector("character", length(c.nms))
  s.nms <- c("LowAIC", "ORmtp", "OR10", "AUCmtp", "AUC10", "^AUC", "^AvgAIC", "^EBPM", "^WAAUC", "^ESORIC")
  c.nms2 <- unlist(sapply(seq_along(s.nms), function(i, x, y, z){
    si <- grepl(s.nms[i], c.nms)
    if(sum(si)>0){
      c.nms2[si] <- gsub("\\^|^\\.", "", paste(c.nms2[si], s.nms[i], sep = "."))
    }
    return(c.nms2[si])
  }, c.nms, s.nms, c.nms2))
  rep.nm <- find_repeated_characters(gsub(paste(unique(c.nms2), collapse = "|"), "", c.nms))
  masks <- gsub(paste(c(rep.nm, paste(unique(c.nms2), collapse = "|")), collapse = "|"), "", c.nms)
  masks[masks==""] <- "all"
  rep.mdl <- length(c.nms2)/length(unique(c.nms2))

  scn.nm <- names(mrs[[1]][[1]])
  areas <- array(dim=c(raster::nlayers(mrs[[1]][[1]]), # sheet (3rd dim) for scn diff
                       length(mrs[[1]]), # cols for threshold criteria
                       length(mrs)), # rows for model criteria
                 dimnames = list(Clim.scen.change = scn.nm, # scn diff
                                 threshold = names(mrs[[1]]), # threshold criteria
                                 Model = names(mrs))) # model criteria
  areas.g <- areas.l <- areas
  thrshld.crit <- names(mrs[[1]])

  zref <- matrix(c(-1,0,1, 0,0,0), nrow = 3, dimnames = list(NULL,c("zone", "sum")))

  for(m in seq_along(mrs)){
    for(t in seq_along(mrs[[m]])){
     for(sc in 1:raster::nlayers(mrs[[m]][[t]])){
        ar <- mrs[[m]][[t]][[sc]]

        if(is.null(area.raster)){
          area.raster <- raster::area(ar, na.rm=TRUE)
        }
        if(isTRUE(all.equal(raster::extent(ar), raster::extent(area.raster)))){
          area.raster <- raster::crop(area.raster, ar)
        }

        zstat <- raster::zonal(area.raster, ar, "sum", digits=digits)

        missrow <- !zref[,1] %in% zstat[,1]
        zstat <- rbind(zref[missrow,],
              zstat[,])
        zstat[] <- zstat[order(zstat[,1]),]

        areas.l[sc,t,m] <- empty2zero(zstat[zstat[,1]==-1, 2])
        areas[sc,t,m] <- empty2zero(zstat[zstat[,1]==0, 2])
        areas.g[sc,t,m] <- empty2zero(zstat[zstat[,1]==1, 2])
     }
    }
  }

  # result <- data.table::rbindlist(lapply(1:3, function(i, a, Change, loc){
  #   cbind(reshape2::melt(a[[i]], value.name = "SuitArea"),
  #         Location = loc,
  #         Change = Change[[i]])
  # }, a = list(areas.l, areas, areas.g), Change = c("loss", "unchanged", "gain"),
  # loc = rep(unique(masks), each = length(areas)/rep.mdl)))
  result <- data.table::rbindlist(lapply(1:3, function(i, a, Change, loc){
    cbind(data.frame(expand.grid(Clim.scen = scn.nm, # pred.scenario
                                 Threshold = names(mrs[[1]]), # threshold criteria
                                 Model = rep(unique(c.nms2), rep.mdl)), # model criteria
                     Location = loc,
                     SuitArea = array(a[[i]])),
          Change = Change[[i]]) # check order of scn.nm and c.nms2
  }, a = list(areas.l, areas, areas.g),
  Change = c("loss", "unchanged", "gain"),
  loc = rep(unique(masks), each=length(areas)/rep.mdl)))
  result$Model <- gsub(paste0("_", masks, collapse = "|"), "", result$Model)

  return(result)
}


#' Compute changes in total suitable area for multiple species
#'
#' Compute changes in total suitable area at multiple climatic scenarios,
#' threshold and model criteria for each species.
#'
#' @param mrs.l List containing stack or brick of range shifts. See \code{\link{range_shift_b}}
#' @inheritParams thrshld_b
#' @inheritParams get_rsa
#' @seealso \code{\link{range_shift}}, \code{\link{get_rsa}}
#' @return An object of class "data.table" "data.frame" containing calculated areas of
#'  range shift for each species between a selected climatic scenario and all other climatic scenarios
#'  for threshold, model and masked locations (if projections were masked,
#'  see \code{\link{mask_thr_projs_mscn_b}}).
#' @examples
#' \dontrun{
#' spp_diff_a <- get_rsa_b(mods.thrshld.lst)
#' }
#' @export
get_rsa_b <- function(mrs.l, area.raster=NULL, digits=2, numCores=1){

  area.occ.spp <- lapply(seq_along(mrs.l), function(i, mrs.l, area.raster, digits){
    get_rsa(mrs.l[[i]], area.raster, digits)
  }, mrs.l, area.raster, digits)

  names(area.occ.spp) <- names(mrs.l)

  area.occ.spp.c <- data.table::rbindlist(area.occ.spp, idcol = "Taxon")
  utils::write.csv(area.occ.spp.c, paste0("3_out.MaxEnt/metric.rangeShiftArea.csv")) # reorder ds

  return(area.occ.spp.c)
}
