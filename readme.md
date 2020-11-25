ENMwizard
======================
### Advanced Tecniques for Ecological Niche Modeling Made Easy

This package provides tools to facilitate the use of advanced techiques related to Ecological Niche Modeling (ENM) and the automation of repetitive tasks (when modeling several species). This package has functions for easier: 1. preparation of occurrence and environmental data; 2. model tunning (thanks to the package ENMeval); 3. model fitting and projection. Tunning and projection can be peformed using a single or multiple cores to speed up processing for multiple species. ENMwizard also implements AICc Model Averaging for MaxEnt models (Gutierrez & Heming, 2018, https://arxiv.org/abs/1807.04346).

-----

# Installation
ENMwizard is downloadable from https://github.com/HemingNM/ENMwizard. You can download it using devtools to install from GitHub.

## Install from GitHub using devtools
Run the following code from your R console:

```r
install.packages("devtools")
devtools::install_github("HemingNM/ENMwizard")

library(ENMwizard)
```

## Citation
Please cite ENMwizard (and other R packages it depends on) by using:

```r
citation("ENMwizard")
citation("spThin")
citation("ENMeval")
citation("raster")
```


-----

# Steps for niche modeling using ENMwizard

## Prepare environmental data

### Load occurrence data

First, lets use occ data available in dismo package.
```r
Bvarieg.occ <- read.table(paste(system.file(package="dismo"),
"/ex/bradypus.csv", sep=""), header=TRUE, sep=",")

head(Bvarieg.occ) # Check first rows

```

Now we make it a named list, where names correspond to species names.
```r
spp.occ.list <- list(Bvarieg = Bvarieg.occ)
```

### Create occ polygon to crop rasters prior to modelling

The occurrence points in the named list are used to create polygons. 
Notice that you can cluster the occ points using several clustering methods. 
See differences and choose one that fits your needs:
```r
occ.polys <- set_calibarea_b(spp.occ.list)
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="AP", q=.01) # less polygons
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="AP", q=.3)
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="AP", q=.8) # more polygons
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="NB", method = "centroid", index = "duda")
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="NB", method = "centroid", index = "sdindex") 

```

### Create buffer

... and the occurrence polygons are buffered using 1.5 degrees.
```r
occ.b <- buffer_b(occ.polys, width = 1.5)
```

### Get and cut enviromental layers
Get climate data for historical (near current) conditions.
In this example, a directory called 'rasters' is created. Then, rasters from historical (near current) are downloaded.
```r
# Create directory to store raster files
dir.create("./rasters")

# Download data for present
library(raster)
predictors <- getData('worldclim', var='bio', res=10, path="rasters")
```

Cut environmental variables for each species (and plot them for visual inspection).
```r
pred.cut <- cut_calibarea_b(occ.b, predictors)

for(i in 1:length(pred.cut)){
  plot(pred.cut[[i]][[1]])
  plot(occ.polys[[i]], border = "red", add = T)
  plot(occ.b[[i]], add = T)
}
```

### Select the least correlated variables
```r
vars <- select_vars_b(pred.cut, cutoff=.75, names.only = T)
# See selected variables for each species
lapply(vars, function(x)x[[1]])
# remove correlated variables from our variable set
pred.cut <- select_vars_b(pred.cut, cutoff=.75, names.only = F)
```


## Prepare occurrence data
### Filter original dataset
Now we want to remove localities that are too close apart. We will do it for all species listed in "spp.occ.list".
```r
thinned.dataset.batch <- thin_b(loc.data.lst = spp.occ.list)
```

### Load occurrence data (filtered localities)
After thinning, we choose one dataset for each species for modelling.
```r
occ.locs <- load_thin_occ(thinned.dataset.batch)
```

## Great! Now we are ready for tunning species' ENMs

-----
## Tunning Maxent's feature classes and regularization multiplier via ENMeval
### Model tuning using ENMeval
Here we will run ENMevaluate_b to call ENMevaluate (from ENMeval package). Here we will test which combination of Feature Classes and Regularization Multipliers give the best results. For this, we will partition our occurrence data using the "block" method.

By providing [at least] two lists, occurrence and environmental data, we will be able to evaluate ENMs for as many species as listed in our occ.locs object. For details see ?ENMeval::ENMevaluate. Notice that you can use multiple cores for this task. This is specially usefull when there are a large number of models and species.
```r
ENMeval.res.lst <- ENMevaluate_b(occ.locs, pred.cut, 
                    RMvalues = c(1, 1.5), fc = c("L", "LQ", "LP"),
                    method="block", algorithm="maxent.jar")
```

-----
## Model fitting (calibration)
After tuning MaxEnt models, we will calibrate them using all occurrence data (i.e. without partition them).

```r

# Run model
mxnt.mdls.preds.lst <- calib_mdl_b(ENMeval.o.l = ENMeval.res.lst, 
                                    a.calib.l = pred.cut, occ.l = occ.locs,
                                    mSel = c("LowAIC", "AUC"))
```

##  Projection
### Prepare projecion area
#### Download environmental data
For projection it is necessary to download raster files with the environmnetal variables of interest. Rasters with historical (near current) climatic conditions was already created. We will download data of climatic conditions for two future (2050 and 2070) scenarios and create one list with all three climate cenarios.

```r
library(raster)
# Get climate data for future conditions (2050)
futAC5085 <- getData('CMIP5', var='bio', res=10, rcp=85, model='AC', year=50,path="rasters")
names(futAC5085) <- names(predictors)

# Get climate data for future conditions (2070)
futAC7085 <- getData('CMIP5', var='bio', res=10, rcp=85, model='AC', year=70,path="rasters")
names(futAC7085) <- names(predictors)

predictors.l <- list(ncurrent = predictors,
                futAC5085 = futAC5085,
                futAC7085 = futAC7085)

```

#### Select area for projection based on the extent of occ points
Now it is time to define the projection area for each species. The projection area can be the same for all species (in this example) of be defined individually. Here, the projection area will be defined as an square area slightly larger than the original occurrence of the species. Then, a two lists with models will be created for a species. In the first list, the projection will be performed using current climatic conditions. In the second list, two cenarios of futurure climate (defined above) are created.

```r
poly.projection <- set_projarea_b(occ.polys, mult = .1, buffer=FALSE)#
plot(poly.projection[[1]], col="gray")
plot(occ.polys[[1]], col="yellow", add=T)

pred.cut.l <- cut_projarea_mscn_b(poly.projection, predictors.l)
plot(poly.projection[[1]], col="gray")
plot(pred.cut.l[[1]][[1]][[1]], add=T)
plot(occ.polys[[1]], add=T)
```

#### ... if the extent to project is the same for all species
When all species are to be projected using the same current and future climates and in the same region, then the following lines can be used to repeat the same lists of cenarios for all species (could be defined differently for each species if wanted).

```r
proj.extent <- extent(c(-109.5, -26.2, -59.5, 18.1))
pred.cut.l <- cut_projarea_rst_mscn_b(proj.extent, predictors.l, occ.polys)
```

### Model projections

Finally, the model(s) can be projected on all climatic cenarios. This is performed by `the proj_mdl_b` function. The function has two arguments: 1) MaxEnt fitted models (see step 4.3 above) and 2) list of rasters representing all cenarios onto which models will be projected.
This function can be run using a single core (default) or multiple cores available in a computer. There two ways of performing parallel processing: by species or by model. If the distribution of few species is being modelled, and models are computationally intensive, then processing by model will provide best results. If there are many species, probably parallel processing by species (split species across the multiple cores of a computer) will be faster.

```r
# For single or multiple species

# using a single core (default)
mxnt.mdls.preds.cf <- proj_mdl_b(mxnt.mdls.preds.lst, a.proj.l = pred.cut.l)

# or using multiple cores
mxnt.mdls.preds.cf <- proj_mdl_b(mxnt.mdls.preds.lst, a.proj.l = pred.cut.l, numCores=2)

# plot projections
par(mfrow=c(1,2), mar=c(1,2,1,2))
plot(mxnt.mdls.preds.cf$Bvarieg$mxnt.preds$ncurrent)
plot(mxnt.mdls.preds.cf$Bvarieg$mxnt.preds$futAC5085)
```

### Apply thresholds on suitability projections
We have the projections for each climatic scenario, now we must select one (or more) threshold criteria and apply on the projections.
```r
# 1. Fixed.cumulative.value.1 (fcv1);
# 2. Fixed.cumulative.value.5 (fcv5);
# 3. Fixed.cumulative.value.10 (fcv10);
# 4. Minimum.training.presence (mtp);
# 5. 10.percentile.training.presence (x10ptp);
# 6. Equal.training.sensitivity.and.specificity (etss);
# 7. Maximum.training.sensitivity.plus.specificity (mtss);
# 8. Balance.training.omission.predicted.area.and.threshold.value (bto);
# 9. Equate.entropy.of.thresholded.and.original.distributions (eetd).

mods.thrshld.lst <- thrshld_b(mxnt.mdls.preds.cf, thrshld.i = c(5,7))
```

## Visualize
### Plot one projection for current climate and another for a future climatic scenario
```r
plot(mods.thrshld.lst$Bvarieg$ncurrent$binary$x10ptp)
plot(mods.thrshld.lst$Bvarieg$futAC5085$binary$x10ptp)
plot_mdl_diff(mxnt.mdls.preds.lst[[1]], mods.thrshld.lst[[1]], sp.nm = "Bvarieg")
plot_mdl_diff_b(mxnt.mdls.preds.cf, mods.thrshld.lst, save=T)
```

### Plot differences between current climate and future climatic scenarios for all thresholds
```r
plot_scn_diff_b(mxnt.mdls.preds.cf, mods.thrshld.lst, 
              ref.scn = "ncurrent", mSel = "LowAIC", save=F)
```


## Compute metrics
### Compute variable contribution and permutation importance
```r
get_cont_permimport_b(mxnt.mdls.preds.cf)
```
### Compute "Fractional Predicted Area" ('n of occupied pixels'/n)
```r
get_fpa_b(mods.thrshld.lst)
```
### Compute species' total suitable area
```r
get_tsa_b(mods.thrshld.lst)
```
