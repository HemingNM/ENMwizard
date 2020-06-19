# #### 5.3 comparar as distribuições (rasteres de adequabilidade) geradas por diferentes critérios de
# ## seleção de modelo (AvgAIC, LowAIC, avg.test.orMTP, avg.test.or10pct, avg.test.AUC.MTP, avg.test.AUC10pct),
# ## usando função ENMTools::raster.overlap(r1, r2).
#' Raster overlap between models selected using different criteria
#'
#' Measures overlap between two ENMs. Used to compare differences among model selection criteria.
#' See ?ENMTools::raster.overlap for details.
#'
#' @inheritParams thrshld_b
#' @param scn.nm name of climatic scenarios to compute overlap.
#' @param model.compare Reference model to be compared (AvgAIC, LowAIC, avg.test.orMTP,
#'  avg.test.or10pct, avg.test.AUC.MTP, avg.test.AUC10pct)
# #' @seealso \code{\link[ENMTools]{raster.overlap}}, \code{\link{cSArea}}, \code{\link{cVarCI}}, \code{\link{cOR}}, \code{\link{cFPA}}
#' @return A list of matrices containing the three metrics (I, D, and Spearman rank correlation)
#' for each climatic scenario
# #' @examples
# #' mxnt.mdls.ovlp <- raster_overlap_b(mxnt.mdls.preds.pf, scn.nm="current", 1)
# #' mxnt.mdls.ovlp <- raster_overlap_b(mxnt.mdls.preds.pf,
# #' scn.nm=c("futAC5085", "futAC7085"), 3)
# #' @export
raster_overlap_b <- function(mcmp.l, scn.nm="current", model.compare=1){print("function disabled until ENMtools is fixed")}
# raster_overlap_b <- function(mcmp.l, scn.nm="current", model.compare=1){
#
#   ovr.spp <- lapply(names(mcmp.l), function(sp, mcmp.l, comb.plots){ # , ovr.vl
#     ovr.vl <- array(dim= c(length(clim.scn), # rows for pred.scenario
#                            ncol(comb.plots), # cols for combinations among model selection criteria
#                            3),  # sheets for metrics
#                     dimnames = list(names(mcmp.l[[sp]]$mxnt.preds[clim.scn]),
#                                     gsub("Mod.", "", paste0(names(mcmp.l[[sp]]$mxnt.preds[[clim.scn[1]]])[comb.plots[1,]],".vs.",
#                                                             names(mcmp.l[[sp]]$mxnt.preds[[clim.scn[1]]])[comb.plots[2,]])),
#                                     c("D", "I", "rank.cor"))) #, # sheet (3rd dim) for model criteria
#
#     clim.scn <- grep(paste(scn.nm, collapse = "|"), names(mcmp.l[[sp]]$mxnt.preds))
#     comb.plots <- utils::combn(raster::nlayers(mcmp.l[[sp]]$mxnt.preds[[scn.nm[1]]]), 2)
#     comb.plots <- comb.plots[,comb.plots[1,]==model.compare| comb.plots[2,]==model.compare]
#     comb.plots[,comb.plots[2,]==model.compare] <- comb.plots[c(2,1),comb.plots[2,]==model.compare]
#
#     r.o.l <- lapply(scn.nm, function(sc, mcmp.l, comb.plots, ovr.vl, sp){
#       mods.comp <- mcmp.l[[sp]]$mxnt.preds[[sc]]
#
#       r.o.s <- sapply(1:ncol(comb.plots), function(j, mods.comp, comb.plots, sp,sc){
#         r.o <- ENMTools::raster.overlap(mods.comp[[comb.plots[1,j]]], mods.comp[[comb.plots[2,j]]])
#         return(r.o)}, mods.comp, comb.plots, sp,sc)
#       return(r.o.s)}, mcmp.l, comb.plots, ovr.vl, sp)
#     r.o.l <- array(simplify2array(r.o.l), dim= dim(ovr.vl), dimnames = dimnames(ovr.vl))
#     utils::write.csv(r.o.l[,,"D"], paste0("3_out.MaxEnt/Mdls.", sp, "/Raster.ovrlp.D", sp, ".csv"))
#     utils::write.csv(r.o.l[,,"I"], paste0("3_out.MaxEnt/Mdls.", sp, "/Raster.ovrlp.I", sp, ".csv"))
#     utils::write.csv(r.o.l[,,"rank.cor"], paste0("3_out.MaxEnt/Mdls.", sp, "/Raster.ovrlp.rank.cor", sp, ".csv"))
#     return(r.o.l)}, mcmp.l, comb.plots)
#
#   names(ovr.spp) <- names(mcmp.l)
#   return(ovr.spp)
# }

# raster_overlap_b(mcmp.l, scn.nm=c("futAC5085", "futAC7085"), 3)

