#' Deprecated & Defunct Functions
#'
#' These functions are Deprecated or Defunct in this release of ENMwizard.
#' Deprecated functions will be marked as Defunct and removed in a future version.
#' Defunct functions were already removed from this version.
#' @name ENMwizard-deprecated-defunct
#' @keywords internal
NULL

##### 1.1.
#' Create minimum convex polygon based on species occurence data
#'
#' 'poly.c' is defunct. Use 'set_calibarea' instead.
#'
#' @param ... additional arguments
#' @rdname ENMwizard-deprecated-defunct
#' @export
poly.c <- function(...){
  .Defunct("'set_calibarea'")
}

#' Create minimum convex polygon based on coordinates of species occurence data for several species
#'
#' 'poly.c.batch' is defunct. Use 'set_calibarea_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
poly.c.batch <- function(...){
  .Defunct("'set_calibarea_b'")
}




#' Create buffer based on species polygons
#'
#'  'bffr.batch' is defunct. Use 'buffer_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
bffr.batch <- function(...){
  .Defunct("'buffer_b'")
}



#' Cut calibration area based on a list of SpatialPolygons
#'
#'  'env.cut' is defunct. Use 'cut_calibarea_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
env.cut <- function(...){
  .Defunct("'cut_calibarea_b'")
}


#' Spatially thin a list of species occurrence data
#'
#'  'thin.batch' is defunct. Use 'thin_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
thin.batch <- function(...){
  .Defunct("'thin_b'")
}


#' Load filtered occurrence data
#'
#'  'loadTocc' is defunct. Use 'load_thin_occ' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
loadTocc <- function(...){
  .Defunct("'load_thin_occ'")
}



##### 1.2.

#' Find, optionally remove, highly correlated variables from a list (several species) of raster brick/stack
#'
#'  'sel_env_b' is defunct. Use 'select_vars_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
sel_env_b <- function(...){
  .Defunct("'select_vars_b'")
}

#' Find, optionally remove, highly correlated variables from a raster brick/stack
#'
#' 'sel_env' is defunct. Use 'select_vars' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
sel_env <- function(...){
  .Defunct("'select_vars'")
}


##### 2.
#' Select area for projection based on the extent of occ points
#'
#' 'pred.a.poly' is defunct. Use 'set_calibarea' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
pred.a.poly <- function(...){
  .Defunct("'set_projarea'")
}


#' Select area for projection based on the extent of occ points for multiple species
#'
#' 'pred.a.poly.batch' is defunct. Use 'set_projarea_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
pred.a.poly.batch <- function(...){
  .Defunct("'set_projarea_b'")
}



#' Cut area for projection based on a list of SpatialPolygons
#'
#' 'pred.a' is defunct. Use 'cut_projarea_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
pred.a <- function(...){
  .Defunct("'cut_projarea_b'")
}

#' Cut area for projection based on a list of SpatialPolygons
#'
#' 'pred.a.batch' is defunct. Use 'cut_projarea_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
pred.a.batch <- function(...){
  .Defunct("'cut_projarea_b'")
}

#' Cut multiple projection areas (climatic scenarios) for multiple species (list of SpatialPolygons)
#'
#' 'pred.a.batch.mscn' is defunct. Use 'cut_projarea_mscn_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
pred.a.batch.mscn <- function(...){
  .Defunct("'cut_projarea_mscn_b'")
}

#' Cut a projection area based on a SpatialPolygon (e.g. Ecoregion)
#'
#' 'pred.a.rst' is defunct. Use 'cut_projarea_rst' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
pred.a.rst <- function(...){
  .Defunct("'cut_projarea_rst'")
}

#' Cut projection areas of multiple species based on a single SpatialPolygon object (e.g. Ecoregion)
#'
#' 'pred.a.batch.rst' is defunct. Use 'cut_projarea_rst_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
pred.a.batch.rst <- function(...){
  .Defunct("'cut_projarea_rst_b'")
}

#' Cut multiple projection areas of multiple species based on a single SpatialPolygon object (e.g. Ecoregion)
#'
#' 'pred.a.batch.rst.mscn' is defunct. Use 'cut_projarea_rst_mscn_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
pred.a.batch.rst.mscn <- function(...){
  .Defunct("'cut_projarea_rst_mscn_b'")
}



##### 3.
#' Tuning and evaluation of ENMs with Maxent for several species using ENMeval
#'
#'  'ENMevaluate.batch' is defunct. Use 'ENMevaluate_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
ENMevaluate.batch <- function(...){
  .Defunct("'ENMevaluate_b'")
}



#' Optimize the size of ENMevaluate_b objects
#'
#' 'ENMevaluate.l.opt' is defunct. Use 'optENMevalObjL' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
ENMevaluate.l.opt <- function(...){
  .Defunct("'optENMevalObjL'")
}



##### 4.
## 4.2 Generate final models for occ
#' Model selection and creation of MaxEnt arguments for selected models
#'
#' 'f.args' is defunct. Use 'mod_sel' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
f.args <- function(...){
  .Defunct("'mod_sel'")
}

#### 4.3 Run top corresponding models and save predictions
#### 4.3.1 save maxent best models and predictions for each model
# "f.mxnt.mdl.pred" renamed to "calib_mdl" and now to "calib_mdl"
#' Calibrate MaxEnt models based on model selection criteria
#'
#' 'mxnt.c' is defunct. Use 'calib_mdl' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
mxnt.c <- function(...){
  .Defunct("'calib_mdl'")
}

#' Calibrate MaxEnt models based on model selection criteria for several species
#'
#' 'mxnt.c.batch' is defunct. Use 'calib_mdl_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
mxnt.c.batch <- function(...){
  .Defunct("'calib_mdl_b'")
}


##### 5.
### 4.8 predictions for future and past
### functions to predict areas based on fitted models

#' Project calibrated MaxEnt models
#'
#' 'mxnt.p' is defunct. Use 'proj_mdl' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
mxnt.p <- function(...){
  .Defunct("'proj_mdl'")
}

#' Project calibrated MaxEnt models
#'
#' 'mxnt.p.batch' is defunct. Use 'proj_mdl' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
mxnt.p.batch <- function(...){
  .Defunct("'proj_mdl'")
}

#' Project calibrated MaxEnt models for several species onto multiple environmental scenarios
#'
#' 'mxnt.p.batch.mscn' is defunct. Use 'proj_mdl_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
mxnt.p.batch.mscn <- function(...){
  .Defunct("'proj_mdl_b'")
}




##### 6.
#### 4.3.3 aplicar threshold
# name of arg "mxnt.mdls.preds.sp[...]" shortened to "mcmp"

#' Apply threshold for MaxEnt projections of a species
#'
#' 'f.thr' is defunct. Use 'thrshld' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
f.thr <- function(...){
  .Defunct("'thrshld'")
}


#### 4.8.5 threshold for past and future pred
#' Apply threshold for MaxEnt projections for multiple species
#'
#' 'f.thr.batch' is defunct. Use 'thrshld_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
f.thr.batch <- function(...){
  .Defunct("'thrshld_b'")
}


##### 7.
#### 4.4 plot differences among model selection criteria predictions
#### 4.8.6 plot prediction diff between models
#' Plot differences in suitable areas between models selected using distinct criteria for multiple species
#'
#' 'plotMdlDiff' is defunct. Use 'plot_mdl_diff_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
plotMdlDiff <- function(...){
  .Defunct("'plot_mdl_diff_b'")
}





#### 4.8.6 plot prediction diff between models
#' Plot differences in suitable areas between climatic scenarios for multiple species
#'
#' 'plotScnDiff' is defunct. Use 'plot_scn_diff_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
plotScnDiff <- function(...){
  .Defunct("'plot_scn_diff_b'")
}



##### 8.1.

#' Compute species' total suitable area
#'
#' 'f.area.occ.mscn' is defunct. Use 'get_tsa_b' or 'get_tsa' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
f.area.occ.mscn <- function(...){
  .Defunct("'get_tsa_b' or 'get_tsa'")
}



# #### 4.7 extract model results
# ### 4.7.1 variable contribution and importance

#' Compute variable contribution and permutation importance
#'
#' 'f.var.ci' is defunct. Use 'get_cont_permimport_b' or 'get_cont_permimport' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
f.var.ci <- function(...){
  .Defunct("'get_cont_permimport_b' or 'get_cont_permimport'")
}


#' Compute "Fractional predicted area" ('n of occupied pixels'/n)
#'
#' 'f.FPA' is defunct. Use 'get_fpa_b' or 'get_fpa' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
f.FPA <- function(...){
  .Defunct("'get_fpa_b' or 'get_fpa'")
}


#' Compute "Omission Rate"
#'
#' 'f.OR' is defunct. Use 'get_OR' or 'get_or_ensemble' or 'get_or_ensemble_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
f.OR <- function(...){
  .Defunct("'get_OR' or 'get_or_ensemble' or 'get_or_ensemble_b'")
}


##### 8.2.

#########	Omission Rate for AICc Averaged Model #############

##### 8.3.



##### 9.
# #### 5.3 comparar as distribuições (rasteres de adequabilidade) geradas por diferentes critérios de
# ## seleção de modelo (AvgAIC, LowAIC, avg.test.orMTP, avg.test.or10pct, avg.test.AUC.MTP, avg.test.AUC10pct),
# ## usando função ENMTools::raster.overlap(r1, r2).

#' Raster overlap between models selected using different criteria
#'
#' 'f.raster.overlap.mscn' is defunct. Use 'raster_overlap_b' instead.
#' @inheritParams poly.c
#' @rdname ENMwizard-deprecated-defunct
#' @export
f.raster.overlap.mscn <- function(...){
  .Defunct("'raster_overlap_b'")
}
