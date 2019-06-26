##### 1.1.
#' Create minimum convex polygon based on species occurence data
#'
#' Deprecated function. Use 'set_calibarea'
#' This function will create minimum convex polygon based on coordinates of species occurence data.
#'
#' @param ... additional arguments
#' @export
poly.c <- function(...){
  stop("Deprecated function. Use 'set_calibarea'")
}

#' Create minimum convex polygon based on coordinates of species occurence data for several species
#'
#' Deprecated function. Use 'set_calibarea_b'
#' This function will use a list of coordinates of species occurence data and create minimum convex polygons
#' for each element in the list.
#' It is possible to create concave or convex polygons, create several small polygons based on clusters
#' of points.
#'
#' @inheritParams poly.c
#' @export
poly.c.batch <- function(...){
  stop("Deprecated function. Use 'set_calibarea_b'")
}




#' Create buffer based on species polygons
#'
#' Deprecated function. Use 'buffer_b'
#' This funcion will create a buffer using species polygons. Buffer width may be manually specified or
#' calculated based on the extent of the SpatialPolygons object (i.e. mean of latitudinal and longitudinal extent).
#' If width is calculated based on the extent of the SpatialPolygons object, it can be adjusted (enlarged or reduced)
#' using 'mult' argument.
#'
#' @inheritParams poly.c
#' @export
bffr.batch <- function(...){
  stop("Deprecated function. Use 'buffer_b'")
}



#' Cut calibration area based on a list of SpatialPolygons
#'
#' Deprecated function. Use 'cut_calibarea_b'
#' Use a list of SpatialPolygons to crop environmental variables for each species.
#'
#' @inheritParams poly.c
#' @export
env.cut <- function(...){
  stop("Deprecated function. Use 'cut_calibarea_b'")
}


#' Spatially thin a list of species occurrence data
#'
#' Deprecated function. Use 'thin_b'
#' Will use spThin optimisation algorithm to subset the dataset such that
#' all occurrence locations are a minimum distance apart. This process helps
#' reduce the effect of biases in observation records on the predictive
#' performance of ecological niche models.
#'
#' Make sure coordinates are in decimal degrees. This function will use
#' great.circle.distance to thin the datasets
#'
#' @inheritParams poly.c
#' @export
thin.batch <- function(...){
  stop("Deprecated function. Use 'thin_b'")
}


#' Load filtered occurrence data
#'
#' Deprecated function. Use 'load_thin_occ'
#' Load filtered occurrence data from object returned by "thin_b" function
#'
#' @inheritParams poly.c
#' @export
loadTocc <- function(...){
  stop("Deprecated function. Use 'load_thin_occ'")
}



##### 1.2.

#' Find, optionally remove, highly correlated variables from a list (several species) of raster brick/stack
#'
#' Deprecated function. Use 'select_vars_b'
#' @inheritParams poly.c
#' @export
sel_env_b <- function(...){
  stop("Deprecated function. Use 'select_vars_b'")
}

#' Find, optionally remove, highly correlated variables from a raster brick/stack
#'
#' Deprecated function. Use 'select_vars'
#' @inheritParams poly.c
#' @export
sel_env <- function(...){
  stop("Deprecated function. Use 'select_vars'")
}


##### 2.
#' Select area for projection based on the extent of occ points
#'
#' Deprecated function. Use 'set_calibarea'
#' This function will create SpatialPolygons that will be used to crop/mask raster/brick objects to be used on model projections. It has several options
#' The user can crop a squared area with an extent larger than the extent of occ.poly. By default, the "extent increase"
#' is the maximum of latitudinal and longitudinal extent "max(abs(ext.proj[1]-ext.proj[2]), abs(ext.proj[3]-ext.proj[4]))".
#' The result is added to each side of the occ.poly extent. This may be changed by setting "mult", which will be multiplied
#' by the "extent increase". Latitudinal and longitudinal increase may also vary independently by setting "same=FALSE".
#'
#' The user can also set a buffer around occ.poly to cut the raster/brick object. Buffer value is, by default,
#' calculated in the same way as "extent increase", using "max(abs(ext.proj[1]-ext.proj[2]), abs(ext.proj[3]-ext.proj[4]))".
#' But the exact value can be defined through "deg.incr" and "mult". This method takes longer to run.
#'
#' @inheritParams poly.c
#' @export
pred.a.poly <- function(...){
  stop("Deprecated function. Use 'set_projarea'")
}


#' Select area for projection based on the extent of occ points for multiple species
#'
#' Deprecated function. Use 'set_projarea_b'
#' This function is a wrapper for "cut_projarea". See ?cut_projarea
#' It works with a named list of occ.polys to delimit the projection area for each of the species.
#'
#' @inheritParams poly.c
#' @export
pred.a.poly.batch <- function(...){
  stop("Deprecated function. Use 'set_projarea_b'")
}



#' Cut area for projection based on a list of SpatialPolygons
#'
#' Deprecated function. Use 'cut_projarea_b'
#' This function is a wrapper for "cut_projarea". See ?cut_projarea
#' It works with a named list of pred.poly to delimit the projection area for each of the species.
#'
#' @inheritParams poly.c
#' @export
pred.a <- function(...){
  stop("Deprecated function. Use 'cut_projarea_b'")
}

#' Cut area for projection based on a list of SpatialPolygons
#'
#' Deprecated function. Use 'cut_projarea_b'
#' This function is a wrapper for "cut_projarea". See ?cut_projarea
#' It works with a named list of pred.poly to delimit the projection area for each of the species.
#'
#' @inheritParams poly.c
#' @export
pred.a.batch <- function(...){
  stop("Deprecated function. Use 'cut_projarea_b'")
}

#' Cut multiple projection areas (climatic scenarios) for multiple species (list of SpatialPolygons)
#'
#' Deprecated function. Use 'cut_projarea_mscn_b'
#' This function is a wrapper for "cut_projarea". See ?cut_projarea. This function delimits the projection area for
#' each of the species contained in the pred.polys named list and crops
#' multiple rasters/bricks (i.e. representing distinct climatic scenaries) based on the same criteria for each species.
#' @inheritParams poly.c
#' @export
pred.a.batch.mscn <- function(...){
  stop("Deprecated function. Use 'cut_projarea_mscn_b'")
}

#' Cut a projection area based on a SpatialPolygon (e.g. Ecoregion)
#'
#' Deprecated function. Use 'cut_projarea_rst'
#' This function will use a single SpatialPolygon to crop/mask raster/brick objects to be used on model projections.
#'
#' @inheritParams poly.c
#' @export
pred.a.rst <- function(...){
  stop("Deprecated function. Use 'cut_projarea_rst'")
}

#' Cut projection areas of multiple species based on a single SpatialPolygon object (e.g. Ecoregion)
#'
#' Deprecated function. Use 'cut_projarea_rst_b'
#' This function will use a single SpatialPolygon to crop/mask raster/brick objects to be used on model projections.
#'
#' @inheritParams poly.c
#' @export
pred.a.batch.rst <- function(...){
  stop("Deprecated function. Use 'cut_projarea_rst_b'")
}

#' Cut multiple projection areas of multiple species based on a single SpatialPolygon object (e.g. Ecoregion)
#'
#' Deprecated function. Use 'cut_projarea_rst_mscn_b'
#' This function will use a single SpatialPolygon to crop/mask multiple raster/brick objects
#' to be used on model projections.
#'
#' @inheritParams poly.c
#' @export
pred.a.batch.rst.mscn <- function(...){
  stop("Deprecated function. Use 'cut_projarea_rst_mscn_b'")
}



##### 3.
#' Tuning and evaluation of ENMs with Maxent for several species using ENMeval
#'
#' Deprecated function. Use 'ENMevaluate_b'
#' This function is a wrapper for ENMeval::ENMevaluate. See ?ENMeval::ENMevaluate for details
#' It works with a named list of species occurrence data (occ.l) and a list of
#' cropped environmental variables (a.calib.l) for model tuning.
#'
#' @inheritParams poly.c
#' @export
ENMevaluate.batch <- function(...){
  stop("Deprecated function. Use 'ENMevaluate_b'")
}



#' Optimize the size of ENMevaluate_b objects
#'
#' Deprecated function. Use 'optENMevalObjL'
#' This function will set to NULL (erase) the largest slots ('predictions', 'models', 'occ.grp', and 'bg.grp')
#' of ENMeval::ENMevaluate objects. Only results, occ.pts, and bg.pts are returned.
#' Use with care. It will not be possible to check MaxEnte models, predictions and grouping of occ and bg points.
#' Can be used to optimize allocated RAM memory when 'ENMevaluate' objects are too large.
#'
#' @inheritParams poly.c
#' @export
ENMevaluate.l.opt <- function(...){
  stop("Deprecated function. Use 'optENMevalObjL'")
}



##### 4.
## 4.2 Generate final models for occ
#' Model selection and creation of MaxEnt arguments for selected models
#'
#' Deprecated function. Use 'mod_sel'
#' This function will read an object of class ENMevaluation (See ?ENMeval::ENMevaluate for details) and
#' return the results table with models selected by the chosen criteria. It can also return the
#' necessary arguments for final model calibration and predictions.
#'
#' @inheritParams poly.c
#' @export
f.args <- function(...){
  stop("Deprecated function. Use 'mod_sel'")
}

#### 4.3 Run top corresponding models and save predictions
#### 4.3.1 save maxent best models and predictions for each model
# "f.mxnt.mdl.pred" renamed to "calib_mdl" and now to "calib_mdl"
#' Calibrate MaxEnt models based on model selection criteria
#'
#' Deprecated function. Use 'calib_mdl'
#' This function will read an object of class ENMevaluation (See ?ENMeval::ENMevaluate for details) and
#' calibrate the selected maxent models.
#'
#' @inheritParams poly.c
#' @export
mxnt.c <- function(...){
  stop("Deprecated function. Use 'calib_mdl'")
}

# "f.mxnt.mdl.pred.batch" renamed to "calib_mdl_b" and now to "calib_mdl_b"


#' Calibrate MaxEnt models based on model selection criteria for several species
#'
#' Deprecated function. Use 'calib_mdl_b'
#' This function will read a list of objects of class ENMevaluation (See ?ENMeval::ENMevaluate for details) and
#' return selected maxent model calibrations and predictions. Each element on the list is usually a species.
#' a.proj.l, a.calib.l, occ.l are lists with occurence data, projection and calibration/predictor data.
#' Species in these lists must all be in the same order of species in ENMeval.o.
#'
#' @inheritParams poly.c
#' @export
mxnt.c.batch <- function(...){
  stop("Deprecated function. Use 'calib_mdl_b'")
}


##### 5.
### 4.8 predictions for future and past
### functions to predict areas based on fitted models

#' Project calibrated MaxEnt models
#'
#' Deprecated function. Use 'proj_mdl'
#' This function will read an object returned by "calib_mdl", read the calibrated models and project into
#' new areas/climatic scenarios. These new projections will be returned together with (appended to)
#' the original object.
#'
#' @inheritParams poly.c
#' @export
mxnt.p <- function(...){
  stop("Deprecated function. Use 'proj_mdl'")
}

#' Project calibrated MaxEnt models
#'
#' Deprecated function. Use 'proj_mdl'
#' This function will read an object returned by "calib_mdl", read the calibrated models and project into
#' new areas/climatic scenarios. These new projections will be returned together with (appended to)
#' the original object.
#'
#' @inheritParams poly.c
#' @export
mxnt.p.batch <- function(...){
  stop("Deprecated function. Use 'proj_mdl'")
}

#' Project calibrated MaxEnt models for several species onto multiple environmental scenarios
#'
#' Deprecated function. Use 'proj_mdl_b'
#' This function will read an object returned by "calib_mdl_b", read the calibrated models and project into
#' several environmental (areas/climatic) scenarios (specified in a.proj.l). These new projections will be returned together with (appended to)
#' each element (species) the original object.
#'
#' @inheritParams poly.c
#' @export
mxnt.p.batch.mscn <- function(...){
  stop("Deprecated function. Use 'proj_mdl_b'")
}




##### 6.
#### 4.3.3 aplicar threshold
# name of arg "mxnt.mdls.preds.sp[...]" shortened to "mcmp"

#' Apply threshold for MaxEnt projections of a species
#'
#' Deprecated function. Use 'thrshld'
#' This function will apply the selected threshold criterias to MaxEnt model projection(s) of a 'mcmp' object
#' and save on the folder "3_out.MaxEnt/Mdls.[species name]/Mdls.thrshld". For each projection (species and climatic
#' scenario), two layers will be generated, one with suitability above the threshold value and another with presence/absence only.
#'
#' @inheritParams poly.c
#' @export
f.thr <- function(...){
  stop("Deprecated function. Use 'thrshld'")
}


#### 4.8.5 threshold for past and future pred
#' Apply threshold for MaxEnt projections for multiple species
#'
#' Deprecated function. Use 'thrshld_b'
#' This function will apply the selected threshold criterias to MaxEnt model projection(s) of a 'mcmp.l' object
#' and save on the folder "3_out.MaxEnt/Mdls.[species name]/Mdls.thrshld". For each projection (species and climatic
#' scenario), two layers will be generated, one with suitability above the threshold value and another with presence/absence only.
#'
#' @inheritParams poly.c
#' @export
f.thr.batch <- function(...){
  stop("Deprecated function. Use 'thrshld_b'")
}


##### 7.
#### 4.4 plot differences among model selection criteria predictions
#### 4.8.6 plot prediction diff between models
#' Plot differences in suitable areas between models selected using distinct criteria for multiple species
#'
#' Deprecated function. Use 'plot_mdl_diff_b'
#' Plot differences between predictions of models selected using distinct
#' criteria (e.g. "AvgAIC", "LowAIC", "OR", "AUC") for multiple climatic scenarios and multiple species
#'
# #' @inheritParams f.plot.mxnt.preds
#' @inheritParams poly.c
#' @export
plotMdlDiff <- function(...){
  stop("Deprecated function. Use 'plot_mdl_diff_b'")
}





#### 4.8.6 plot prediction diff between models
#' Plot differences in suitable areas between models selected using distinct criteria
#'
#' Deprecated function. Use 'set_calibarea'
#' Plot differences between predictions of models selected using distinct
#' criteria (e.g. "AvgAIC", "LowAIC", "OR", "AUC") for multiple climatic scenarios
#'
# #' @inheritParams f.plot.mxnt.preds
# #' @export
# plotMdlDiff <- function(...){
#   stop("Deprecated function. Use 'plot_mdl_diff'")
# }




####
#' Plot differences in suitable areas between climatic scenarios for multiple species
#'
#' Deprecated function. Use 'plot_scn_diff_b'
#' Plot differences between a selected climatic scenario and all other climatic scenarios for each species.
#' This function will plota and (optionally) save the figures on pdf files in the folder "Mdls.thrshld/figs".
#' @inheritParams poly.c
#' @export
plotScnDiff <- function(...){
  stop("Deprecated function. Use 'plot_scn_diff_b'")
}



#' Plot differences in suitable areas between climatic scenarios
#'
#' Deprecated function. Use 'set_calibarea'
#' Plot differences between a selected climatic scenario and all other climatic scenarios for each species.
#' This function will plota and (optionally) save the figures on pdf files in the folder "Mdls.thrshld/figs".
#'
# #' @export
# plot_scn_diff <- function(...){
#   stop("Deprecated function. Use 'set_projarea'")
# }




##### 8.1.

#' Compute species' total suitable area
#'
#' Deprecated function. Use 'get_tsa_b' or 'get_tsa'
#' Compute total suitable area at multiple climatic scenario, threshold and model criteria.
#'
#' @inheritParams poly.c
#' @export
f.area.occ.mscn <- function(...){
  stop("Deprecated function. Use 'get_tsa_b' or 'get_tsa'")
}


#' Compute species' total suitable area
#'
#' Deprecated function. Use 'set_calibarea'
#' Compute total suitable area at multiple climatic scenario, threshold and model criteria.
#'
# #' @export
# get_tsa <- function(...){
#   stop("Deprecated function. Use 'set_projarea'")
# }


# #### 4.7 extract model results
# ### 4.7.1 variable contribution and importance

#' Compute variable contribution and permutation importance
#'
#' Deprecated function. Use 'get_cont_permimport_b' or 'get_cont_permimport'
#' Compute variable contribution and importance for each model
#'
#' @inheritParams poly.c
#' @export
f.var.ci <- function(...){
  stop("Deprecated function. Use 'get_cont_permimport_b' or 'get_cont_permimport'")
}


#' Compute variable contribution and permutation importance
#'
#' Deprecated function. Use 'set_calibarea'
#' Compute variable contribution and importance for each model
#'
# #' @export
# get_cont_permimport <- function(...){
#   stop("Deprecated function. Use 'set_projarea'")
# }


#' Compute "Fractional predicted area" ('n of occupied pixels'/n)
#'
#' Deprecated function. Use 'get_fpa_b' or 'get_fpa'
#' Compute "Fractional predicted area" ('n of occupied pixels'/total n) or ('area of occupied pixels'/total area)
#'
#' @inheritParams poly.c
#' @export
f.FPA <- function(...){
  stop("Deprecated function. Use 'get_fpa_b' or 'get_fpa'")
}

#' Compute "Fractional predicted area" ('n of occupied pixels'/n)
#'
#' Compute "Fractional predicted area" ('n of occupied pixels'/total n) or ('area of occupied pixels'/total area)
#'
# #' @export
# get_fpa <- function(...){
#   stop("Deprecated function. Use 'set_projarea'")
# }


#' Compute "Omission Rate"
#'
#' Deprecated function. Use 'get_OR' or 'get_or_ensemble' or 'get_or_ensemble_b'
#' Compute "Omission Rate" of species occurence points for a climatic scenario (usually "current")
#'
#' @inheritParams poly.c
#'@export
f.OR <- function(...){
  stop("Deprecated function. Use 'get_OR' or 'get_or_ensemble' or 'get_or_ensemble_b'")
}


##### 8.2.

#########	Omission Rate for AICc Averaged Model #############
#' Compute Omission Rate for a species' ensembled model
#'
#' This function will compute the omission rate (OR) for a species' ensembled model
#' from a 'mcmp' object, based on the selected threshold value.
#'
# #' @export
# get_or_ensemble <- function(...){
#   stop("Deprecated function. Use 'set_projarea'")
# }



#' Compute Omission Rate for a list of species' ensembled models
#'
#' This function will compute the omission rate (OR) for each species' ensembled model
#' from a 'mcmp.l' object, based on the selected threshold value.
#'
# #' @export
# get_or_ensemble_b <- function(...){
#   stop("Deprecated function. Use 'set_projarea'")
# }


##### 8.3.


##
#' Extrapolation risk analysis
#'
#' This function will compute the omission rate (OR) for a species' AICc Averaged Model
#' from a 'mcmp' object, based on the selected threshold value.
#'
# #' @export
# mop <- function(...){
#   stop("Deprecated function. Use 'set_projarea'")
# }

#' Extrapolation risk analysis for a list of species
#'
#' This function will compute the omission rate (OR) for each species' AICc Averaged Model
#' from a 'mcmp.l' object, based on the selected threshold value.
#'
# #' @export
# mop_b <- function(...){
#   stop("Deprecated function. Use 'set_projarea'")
# }

##### 9.
# #### 5.3 comparar as distribuições (rasteres de adequabilidade) geradas por diferentes critérios de
# ## seleção de modelo (AvgAIC, LowAIC, avg.test.orMTP, avg.test.or10pct, avg.test.AUC.MTP, avg.test.AUC10pct),
# ## usando função ENMTools::raster.overlap(r1, r2).

#' Raster overlap between models selected using different criteria
#'
#' Deprecated function. Use 'raster_overlap_b'
#' Measures overlap between two ENMs. Used to compare differences among model selection criteria.
#' See ?ENMTools::raster.overlap for details.
#'
#' @inheritParams poly.c
#' @export
f.raster.overlap.mscn <- function(...){
  stop("Deprecated function. Use 'raster_overlap_b'")
}
