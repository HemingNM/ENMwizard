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
# create batch function ("poly.splt.batch") for  "poly.splt"
# d <- cbind(spp.occ.list$Mleucophrys2$LONG, spp.occ.list$Mleucophrys2$LAT)
#  # hclust.obj <- hclust(dist(d))
#  apclus <- apcluster::apcluster(apcluster::negDistMat(r=2), d)
#  plot(apclus, d)

# 5.f.mscn.
# function "mxnt.p"
# paralelize line 58-64 & 107-118  ##### list of models to average


## 3.f.tuning
# create class "ENMevaluation.list"

## 4.f.tuning.final
# create class "MaxEnt.list" (list of "MaxEnt" objects)
# create class "ENM.w" with slots:
# "@selected.mdls"  = "data.frame"
# "@MaxEnt.list" (currently "mxnt.mdls")
# "@pred.args" = "character
# "@mxnt.preds" = (list of "RasterBrick" or "RasterStack" objects, containing the predictions)
# maybe # "@mxnt.preds.thrhld" = (list of "RasterBrick" or "RasterStack" objects, containing the thresholded predictions)



