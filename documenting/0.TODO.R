#
# # https://stackoverflow.com/questions/26799630/r-s4-best-practice-slot-with-vector-of-s4-objects
# #  An alternative to implementing your own type-checked list (which is really the other alternative)
# # is to re-use the infrastructure from Biocdonductor's S4Vectors class
#
# X = setClass("X", representation(x="numeric"))
# .XList = setClass("XList", contains="SimpleList",
#                   prototype=prototype(elementType="X"))
#
#
# # And in action
# xl = .XList(listData=list(.X(x=1), .X(x=2)))
# xl
# # XList of length 2
# xl[[2]]
# # An object of class "X"
# # Slot "x":
# #   [1] 2



############# TODO
# 1.f.data.prep
# use species names as data frames when creating spatialPolygonsDataFrames in functions:
# poly.c (lines #59-60), bind.shp (lines #424)

# 4.f.tunning.final
# lines 216 and 232: 'resu <- dismo::maxent'
# Add option to use 'maxnet::maxnet', based on 'algorithm' [line 175] returned from [ENMeval 0.3.0], which allows use of maxnet
# instead of 'dismo::maxent'

# 10. ensemble of selected models
#
#

# 7.f.plot
# A. create function to map uncertainties across models. This function was made for ENMavg paper
# this function might compare projections among models (tuning), thresholds, climatic scenarios
# B. function to create consensual models across climatic scenarios
# C.a function to get thresholds, climatic scenarios, model selection criteria;
# C.b function to plot projections using 'C.a'

## 3.f.tuning
# create class "ENMevaluation.list" # see 'documenting/all.classes.R'

## 4.f.tuning.final
# create class "MaxEnt.list" (list of "MaxEnt" objects)
# create class "ENMw" with slots:
# "@selected.mdls"  = "data.frame"
# "@MaxEnt.list" (currently "mxnt.mdls")
# "@pred.args" = "character
# "@mxnt.preds" = (list of "RasterBrick" or "RasterStack" objects, containing the predictions)
# maybe # "@mxnt.preds.thrhld" = (list of "RasterBrick" or "RasterStack" objects, containing the thresholded predictions)



