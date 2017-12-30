# #### 5.3 comparar as distribuições (rasteres de adequabilidade) geradas por diferentes critérios de
# ## seleção de modelo (AICAvg, AICLow, Mean.ORmin, Mean.OR10, Mean.AUCmin, Mean.AUC10),
# ## usando função ENMTools::raster.overlap(r1, r2).
#' Raster overlap between models selected using different criteria
#'
#' Measures overlap between two ENMs. Used to compare differences among model selection criteria.
#' See ?ENMTools::raster.overlap for details.
#' @inheritParams f.thr.batch
#' @param scn.nm name of climatic scenarios to compute overlap.
#' @param model.compare Reference model to be compared (AvgAICc, LowAICc, Mean.ORmin,
#'  Mean.OR10, Mean.AUCmin, Mean.AUC10)
#' @return A list of matrices containing the three metrics (I, D, and Spearman rank correlation)
#' for each climatic scenario
#' @examples
#' mxnt.mdls.ovlp <- f.raster.overlap.mscn(mxnt.mdls.preds.pf, scn.nm="current", 1)
#' #' mxnt.mdls.ovlp <- f.raster.overlap.mscn(mxnt.mdls.preds.pf, scn.nm=c("futAC5085", "futAC7085"), 3)
#' @export
f.raster.overlap.mscn <- function(mmp.spl, scn.nm="current", model.compare=1){
  clim.scn <- grep(paste(scn.nm, collapse = "|"), names(mmp.spl[[1]]$mxnt.preds))
  comb.plots <- utils::combn(raster::nlayers(mmp.spl[[1]]$mxnt.preds[[scn.nm[1]]]), 2)
  comb.plots <- comb.plots[,comb.plots[1,]==model.compare| comb.plots[2,]==model.compare]
  comb.plots[,comb.plots[2,]==model.compare] <- comb.plots[c(2,1),comb.plots[2,]==model.compare]
  # comb.plots[,comb.plots[2,]==model.compare]<- rev(comb.plots[,comb.plots[2,]==model.compare])

  # ovr.spp <- vector("list", length = length(mmp.spl))
  # names(ovr.spp) <- names(mmp.spl)
  names(area.occ.spp) <- names(mmp.spl)
  ovr.vl <- array(dim= c(length(clim.scn), # rows for pred.scenario
                         ncol(comb.plots), # cols for combinations among model selection criteria
                         3),  # sheets for metrics
                  dimnames = list(names(mmp.spl[[1]]$mxnt.preds[clim.scn]),
                                  gsub("Mod.", "", paste0(names(mmp.spl[[1]]$mxnt.preds[[clim.scn[1]]])[comb.plots[1,]],".vs.",
                                                          names(mmp.spl[[1]]$mxnt.preds[[clim.scn[1]]])[comb.plots[2,]])),
                                  c("D", "I", "rank.cor"))) #, # sheet (3rd dim) for model criteria

  ovr.spp <- lapply(names(mmp.spl), function(sp, mmp.spl, comb.plots, ovr.vl){
    # for(sp in names(mmp.spl)){ # species
    # ovr.spp[[sp]] <- ovr.vl

    r.o.l <- lapply(scn.nm, function(sc, mmp.spl, comb.plots, ovr.vl, sp){
      # for(sc in scn.nm){ # climatic scenari0 # names(mmp.spl[[sp]]$mxnt.preds[clim.scn])
      mods.comp <- mmp.spl[[sp]]$mxnt.preds[[sc]]

      ## TODO use mclapply
      r.o.s <- sapply(1:ncol(comb.plots), function(j, mods.comp, comb.plots, sp,sc){
        # for(j in 1:ncol(comb.plots)){ # ncol(comb.plots)
        r.o <- ENMTools::raster.overlap(mods.comp[[comb.plots[1,j]]], mods.comp[[comb.plots[2,j]]])
        # ovr.spp[[sp]][sc,j,"D"] <- r.o$D # Schoener - climatic
        # ovr.spp[[sp]][sc,j,"I"] <- r.o$I # Hellinger - geographic
        # ovr.spp[[sp]][sc,j,"rank.cor"] <- r.o$rank.cor
        # } # ncol(comb.plots)
        return(r.o)}, mods.comp, comb.plots, sp,sc)

      # }  # climatic scenario
      return(r.o.s)}, mmp.spl, comb.plots, ovr.vl, sp)

    r.o.l <- array(simplify2array(r.o.l), dim= dim(ovr.vl), dimnames = dimnames(ovr.vl))


    xlsx::write.xlsx(r.o.l[,,"D"], paste0("4_ENMeval.results/Mdls.", sp, "/Raster.ovrlp.", sp, ".xlsx"), sheetName="D")
    xlsx::write.xlsx(r.o.l[,,"I"], paste0("4_ENMeval.results/Mdls.", sp, "/Raster.ovrlp.", sp, ".xlsx"), sheetName="I", append=T)
    xlsx::write.xlsx(r.o.l[,,"rank.cor"], paste0("4_ENMeval.results/Mdls.", sp, "/Raster.ovrlp.", sp, ".xlsx"), sheetName="rank.cor", append=T)
    # } # species
    return(r.o.l)}, mmp.spl, comb.plots, ovr.vl)

  names(ovr.spp) <- names(mmp.spl)
  return(ovr.spp)
}

# f.raster.overlap.mscn(mmp.spl, scn.nm=c("futAC5085", "futAC7085"), 3)

