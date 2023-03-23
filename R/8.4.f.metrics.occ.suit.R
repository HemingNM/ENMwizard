#' Extract suitability values from all projections for known species occurrences
#' of a 'mcmp' object
#'
#' This function will extract the suitability values from all projections of a 'mcmp'
#' object using either the species occurrences used for niche modeling or an user
#' supplied object with species occurrences.
#' @inheritParams thrshld
#' @inheritParams calib_mdl
#' @inheritParams consensus_scn
#' @examples
#' \dontrun{
#' get_occ_suit(mcmp=mxnt.mdls.preds)
#' }
#' @export
get_occ_suit <- function(mcmp, occ=NULL, ref="current", t.all=F){
  if(is.null(occ)){
    occ <- mcmp$occ.pts
  }
  npts <- nrow(as.data.frame(occ))
  ID <- 1:npts
  suit <- data.frame(matrix(ncol = 4, nrow = npts))
  colnames(suit) <- c("Clim.scen", "Model", "Suitability", "ID")
  suit.f <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(suit.f) <- colnames(suit)
  ###- get projections to apply thresholds
  # choose between consensus, mxnt.preds, or both
  if(!is.null(mcmp$scn.consensus)){
    if(t.all){ # individual and consensus projections
      scn.nms <- c(names(mcmp$mxnt.preds), names(mcmp$scn.consensus))
    } else { # consensus projections only
      scn.nms <- names(mcmp$scn.consensus)
    }
  } else { # individual projections only
    scn.nms <- names(mcmp$mxnt.preds)
  }
  # scn.nms <- names(mcmp$scn.consensus)
  for(sc in scn.nms){
    # choose between consensus, mxnt.preds, or both
    if(!is.null(mcmp$scn.consensus)){
      if(t.all){ # individual and consensus projections
        if(sc <= length(mcmp$mxnt.preds)){
          pred.r <- mcmp$mxnt.preds[[sc]] # mcmp[[match(pred.nm, names(mcmp))]] # , fixed=TRUE # [pred.i]
          mod <- names(mcmp$mxnt.preds[[sc]])
          } else {
          pred.r <- mcmp$scn.consensus[[sc]]
          mod <- names(mcmp$scn.consensus[[sc]])
        }
      } else { # consensus projections only
        pred.r <- mcmp$scn.consensus[[sc]]
        mod <- names(mcmp$scn.consensus[[sc]])
      }
    } else { # individual projections only
      pred.r <- mcmp$mxnt.preds[[sc]]
      mod <- names(mcmp$mxnt.preds[[sc]])
    }

    for(m in mod){
      suit[,"Suitability"] <- raster::extract(pred.r[[m]], occ)
      suit[,"Model"] <- m
      suit[,"Clim.scen"] <- sc
      suit[,"ID"] <- ID
      suit.f <- rbind(suit.f, suit)
    } # for mod
  } # for scs
  suit.f[ref] <- suit.f$Suitability[suit.f$Clim.scen==ref]
  suit.f$Direction <- "increase"
  suit.f$Direction[suit.f[ref] == suit.f$Suitability] <- "stable"
  suit.f$Direction[suit.f[ref] > suit.f$Suitability] <- "decrease"
  return(suit.f)
}


#' Extract suitability values from all projections for known species occurrences from
#' a 'mcmp.l' object
#'
#' This function will extract the suitability values from all projections of each species
#' in a 'mcmp.l' object using either the species occurrences used for niche modeling or
#' a user supplied list of objects containing species occurrences.
#' @inheritParams thrshld_b
#' @inheritParams calib_mdl_b
#' @inheritParams get_occ_suit
#' @examples
#' \dontrun{
#' get_occ_suit_b(mcmp.l=mxnt.mdls.preds.cf)
#' }
#' @export
get_occ_suit_b <- function(mcmp.l, occ.l=NULL, ref="current", t.all=F){

  suit.occ.spp <- lapply(seq_along(mcmp.l), function(i, mcmp.l, occ.l, ref, t.all){
    get_occ_suit(mcmp.l[[i]], occ.l[[i]], ref=ref, t.all=t.all)
  }, mcmp.l, occ.l, ref, t.all)
  names(suit.occ.spp) <- names(mcmp.l)

  suit.occ.spp.c <- data.table::rbindlist(suit.occ.spp, idcol = "Taxon")

  utils::write.csv(suit.occ.spp.c, paste0("3_out.MaxEnt/metric.occ.suit.csv")) # reorder ds
  return(suit.occ.spp.c)
}
