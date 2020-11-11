ENMwizard
======================
### Advanced Tecniques for Ecological Niche Modeling Made Easy

This package provides tools to facilitate the use of advanced techiques related to Ecological Niche Modeling (ENM) and the automation of repetitive tasks (when modeling several species). This package has functions for easier: 1. preparation of occurence and environmental data; 2. model tunning (thanks to the package ENMeval); 3. model fitting and projection. Tunning and projection can be peformed using a single or multiple cores to speed up processing for multiple species. ENMwizard also implements AICc Model Averaging for MaxEnt models (Gutierrez & Heming, 2018, https://arxiv.org/abs/1807.04346).

Please cite ENMwizard (and other R packages it depends on) by using:
```r
citation("ENMwizard")
citation("spThin")
citation("ENMeval")
citation("raster")
```

-----

# Installation
ENMwizard is downloadable from https://github.com/HemingNM/ENMwizard. You can download it using devtools to install from GitHub.

### Installing from GitHub using devtools
Run the following code from your R console:

```r
install.packages("devtools")
devtools::install_github("HemingNM/ENMwizard")

library(ENMwizard)
```


### Install from zip file
You can also download a zip file containing the package and install it from R.

Download from https://github.com/HemingNM/ENMwizard/archive/master.zip and run the following code (where PATH is the path to the zip file)

```r
install.packages("devtools")
library(devtools)
library(ENMwizard)
```



-----


# Using ENMwizard's magic wand

## 1. Prepare environmental data

### 1.1 Load occurence data

First, lets use occ data available in dismo package
```r
Bvarieg.occ <- read.table(paste(system.file(package="dismo"),
"/ex/bradypus.csv", sep=""), header=TRUE, sep=",")

head(Bvarieg.occ)# Check first rows

# Column names must be in capital letters
# colnames(Bvarieg.occ) <- c("SPEC", "LONG", "LAT") # Change column names
```

Now we make it a named list, where names correspond to species names.
```r
spp.occ.list <- list(Bvarieg = Bvarieg.occ)
```

### 1.2 create occ polygon to crop rasters prior to modelling

The occurence points in the named list are used to create polygons. 
Notice that you can cluster the occ points using several clustering methods. 
See differences and choose one that fits your needs:
```r
occ.polys <- set_calibarea_b(spp.occ.list)
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="E")
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="AP", q=.01) # less polygons
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="AP", q=.1)
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="AP", q=.2)
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="AP", q=.3)
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="AP", q=.5)
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="AP", q=.8) # more polygons
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="NB", method = "centroid", index = "duda")
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="NB", method = "centroid", index = "ball")
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="NB", method = "centroid", index = "pseudot2")
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="NB", method = "centroid", index = "trcovw")
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="NB", method = "centroid", index = "kl") 
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="NB", method = "centroid", index = "sdindex") 
occ.polys <- set_calibarea_b(spp.occ.list, k=0, c.m="NB", method = "centroid", index = "sdbw")
```

### 1.2.1 creating buffer

... and the occurrence polygons are buffered using 1.5 degrees.
```r
occ.b <- buffer_b(occ.polys, width = 1.5)
```

### 1.3. Cut enviromental layers with occ.b and save in hardrive.
Specify the path to the environmental variables
it usually is the path on your machine. E.g. "/path/to/variables/WorldClim/2_5min/bio_2-5m_bil"
here we will use variables available in dismo package
```r
path.env <- paste(system.file(package="dismo"), "/ex", sep="")
biovars <- paste0("bio", 1:17)
pattern.env <- 'grd'
```

Get uncut variables
```r
env.uncut <- list.files(path.env, full.names=TRUE)
env.uncut <- env.uncut[grep(paste(paste0(biovars, ".", pattern.env), collapse = "|"), env.uncut)]
env.uncut <- stack(env.uncut)
```

If the variables were already saved in a raster brick, you just need to read them. Run only if previous commant did not work.
```r
env.uncut <- brick(paste(path.env, "bio.grd", sep="/"))
```

Crop environmental variables for each species (and plot them for visual inspection)
```r
occ.b.env <- cut_calibarea_b(occ.b, env.uncut)

for(i in 1:length(occ.b.env)){
  plot(occ.b.env[[i]][[1]])
  plot(occ.polys[[i]], border = "red", add = T)
  plot(occ.b[[i]], add = T)
}
```

Select the least correlated variables
```r
vars <- select_vars_b(occ.b.env, cutoff=.75, names=T)
# See selected variables for each species
lapply(vars, function(x)x[[1]])
# remove correlated variables from our variable set
occ.b.env <- select_vars_b(occ.b.env, cutoff=.75, names=F)
```


## 2. Prepare occurence data
### 2.1 Filtering original dataset
Now we want to remove localities that are too close apart. We will do it for all species listed in "spp.occ.list".
```r
thinned.dataset.batch <- thin_b(loc.data.lst = spp.occ.list)
```

### Great! Now we are ready for tunning species' ENMs


-----


## 3. Tunning Maxent's feature classes and regularization multiplier via ENMeval
### 3.1 Load occurrence data (filtered localities).
After thinning, we choose one dataset for each species for modelling.
```r
occ.locs <- load_thin_occ(thinned.dataset.batch)
```

### 3.2 model tuning using ENMeval
Here we will run ENMevaluateB to call ENMevaluate (from ENMeval package). Here we will test which combination of Feature Classes and Regularization Multipliers give the best results. For this, we will partition our occurence data using the "block" method.

By providing [at least] two lists, occurence and environmental data, we will be able to evaluate ENMs for as many species as listed in our occ.locs object. For details see ?ENMeval::ENMevaluate. Notice that you can use multiple cores for this task. This is specially usefull when there are a large number of models and species.
```r
ENMeval.res.lst <- ENMevaluate_b(occ.locs, occ.b.env, 
                    RMvalues = c(1, 1.5), fc = c("L", "LQ", "LP"),
                    method="block", algorithm="maxent.jar")
```

-----
## 4. Model Fitting (Calibration)
After tuning MaxEnt models, we will calibrate them using all occurence data (i.e. without partition them).

```r

# Run model
mxnt.mdls.preds.lst <- calib_mdl_b(ENMeval.o.l = ENMeval.res.lst, 
                                    a.calib.l = occ.b.env, occ.l = occ.locs,
                                    mSel = c("AvgAIC", "AUC"))# 


# # Comparing single core processing and multiple core processing
# 
# # Parallel processing forking models (best for small sets of species)
# system.time(
# mxnt.mdls.preds.lst <- calib_mdl_b(ENMeval.o.l = ENMeval.res.lst, a.calib.l = occ.b.env, occ.l=occ.locs, wAICsum=0.99, numCores=3, parallelTunning=TRUE)
# )
# 
# # Parallel processing forking species (best for large sets of species)
# system.time(
# mxnt.mdls.preds.lst <- calib_mdl_b(ENMeval.o.l = ENMeval.res.lst, a.calib.l = occ.b.env, occ.l=occ.locs, wAICsum=0.99, numCores=3, parallelTunning=FALSE)
# )

```

##  5. Projection
### 5.1. Downloading environmental data
For projection it is necessary to download raster files with the environmnetal variables of interest. In this example, a directory called 'rasters' is created. Then, rasters from current and future climatic conditions projected for 2050 and 2070 are downloaded and loaded. Finally, two lists are created, one for current conditions and another for the two future cenarios.

```r

# Create directory to store raster files
dir.create("./rasters")

# Download data for present
current <- getData('worldclim', var='bio', res=10, path="rasters")

# Download data for future projection (2050)
futAC5085 <- getData('CMIP5', var='bio', res=10, rcp=85, model='AC', year=50,path="rasters")
names(futAC5085) <- names(current)

# Download data for future projection (2070)
futAC7085 <- getData('CMIP5', var='bio', res=10, rcp=85, model='AC', year=70,path="rasters")
names(futAC7085) <- names(current)

env.proj.l <- list(current = current,
                futAC5085 = futAC5085,
                futAC7085 = futAC7085)

```

### 5.2 Preparing projecion area: save rasters onto which the model will be projected in an object called "pa.env.proj.l"
### 5.2.1 select area for projection based on the extent of occ points
Now it is time to define the projection area for each species. The projection area can be the same for all species (in this example) of be defined individually. Here, the projection area will be defined as an square area slightly larger than the original occurrence of the species. Then, a two lists with models will be created for a species. In the first list, the projection will be performed using current climatic conditions. In the second list, two cenarios of futurure climate (defined above) are created.

```r
poly.projection <- set_projarea_b(occ.polys, mult = .1, buffer=FALSE)#
plot(poly.projection[[1]], col="gray")
plot(occ.polys[[1]], col="yellow", add=T)

pa.env.proj.l <- cut_projarea_mscn_b(poly.projection, env.proj.l)
plot(poly.projection[[1]], col="gray")
plot(pa.env.proj.l[[1]][[1]][[1]], add=T)
plot(occ.polys[[1]], add=T)
```

### 5.2.2 if the extent to project is the same for all species
When all species are to be projected using the same current and future climates and in the same region, then the following lines can be used to repeat the same lists of cenarios for all species (could be defined differently for each species if wanted)

```r
proj.extent <- extent(c(-109.5, -26.2, -59.5, 18.1))
pa.env.proj.l <- cut_projarea_rst_mscn_b(proj.extent, env.proj.l, occ.polys)
```

### 5.3 Model projections

Finally, the model(s) can be projected on all climatic cenarios. This is performed by `the proj_mdl_b` function. The function has two arguments: 1) MaxEnt fitted models (see step 4.3 above) and 2) list of rasters representing all cenarios onto which models will be projected.
This function can be run using a single core (default) or multiple cores available in a computer. There two ways of performing parallel processing: by species or by model. If the distribution of few species is being modelled, and models are computationally intensive, then processing by model will provide best results. If there are many species, probably parallel processing by species (split species across the multiple cores of a computer) will be faster.

```r
# For single or multiple species

# using a single core (default)
mxnt.mdls.preds.cf <- proj_mdl_b(mxnt.mdls.preds.lst, a.proj.l = pa.env.proj.l)

# or using multiple cores
mxnt.mdls.preds.cf <- proj_mdl_b(mxnt.mdls.preds.lst, a.proj.l = pa.env.proj.l, numCores=2)

# plot projections
par(mfrow=c(1,2), mar=c(1,2,1,2))
plot(mxnt.mdls.preds.cf$Bvarieg$mxnt.preds$current)
plot(mxnt.mdls.preds.cf$Bvarieg$mxnt.preds$futAC5085)
```

### 5.4 Applying thresholds on climatic scenarios
We have the projections for each climatic scenario, now we must select one (or more) threshold criteria and 
apply on the projections.
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

mods.thrshld.lst <- thrshld_b(mxnt.mdls.preds.cf, thrshld.i = 5)
```

## 6. Visualizing
### 6.1. Plotting one projection for current climate and another for a future climatic scenario
```r
plot(mods.thrshld.lst$Bvarieg$current$binary$x10ptp)
plot(mods.thrshld.lst$Bvarieg$futAC5085$binary$x10ptp)
plot_mdl_diff(mxnt.mdls.preds.lst[[1]], mods.thrshld.lst[[1]], sp.nm = "Bvarieg")
plot_mdl_diff_b(mxnt.mdls.preds.cf, mods.thrshld.lst, save=T)
```

### 6.2. Plotting differences between current climate and future climatic scenarios for all thresholds
```r
plot_scn_diff_b(mxnt.mdls.preds.cf, mods.thrshld.lst, 
              ref.scn = "current", mSel = "AvgAIC", save=F)
```


## 7. Metrics
### Compute variable contribution and importance
```r
get_cont_permimport_b(mxnt.mdls.preds.cf)
```
### Compute "Fractional predicted area" ('n of occupied pixels'/n) for multiple scenarios
```r
get_fpa_b(mods.thrshld.lst)
```
### Compute species' total suitable area
```r
get_tsa_b(mods.thrshld.lst)
```
