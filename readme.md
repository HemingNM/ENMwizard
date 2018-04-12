ENMwizard
======================
### Advanced Tecniques for Ecological Niche Modeling Made Easy

This package provides tools to facilitate the use of advanced techiques related to Ecological Niche Modeling (ENM) and the automation of repetitive tasks (when modeling several species). This package has functions for easier: 1. preparation of occurence and environmental data; 2. model tunning (thanks to the package ENMeval); 3. model fitting and projection. Tunning and projection can be peformed using a single or multiple cores to speed up processing for multiple species. ENMwizard also implements methods described in Gutiérrez & Heming in prep.

-----

# Installation
ENMwizard is downloadable from https://github.com/HemingNM/ENMwizard. You can download it using devtools to install from GitHub.

### Installing from GitHub using devtools
Run the following code from your R console:

```r
install.packages("devtools")
devtools::install_github("HemingNM/ENMwizard")
devtools::install_github("danlwarren/ENMTools")
devtools::install_github("mlammens/spThin")
library(ENMwizard)
```

You also need to install ENMTools from GitHub.

```r
install_github("danlwarren/ENMTools")
```


### Install from zip file
You can also download a zip file containing the package and install it from R.

Download from https://github.com/HemingNM/ENMwizard/archive/master.zip and run the following code (where PATH is the path to the zip file)

```r
install.packages("devtools")
library(devtools)
install_local("PATH")
library(ENMwizard)
```



-----


# Using ENMwizard's magic wand

## - 1. Prepare environmental data

### - 1.1 Load occurence data

First, lets use occ data available in dismo package
```r
Bvarieg.occ <- read.table(paste(system.file(package="dismo"),
"/ex/bradypus.csv", sep=""), header=TRUE, sep=",")

head(Bvarieg.occ)# Check first rows

# Column names must be in capital letters
colnames(Bvarieg.occ) <- c("SPEC", "LONG", "LAT") # Change column names
```

Now we make it a named list, where names correspond to species names.
```r
spp.occ.list <- list(Bvarieg = Bvarieg.occ)
```

### - 1.2 create occ polygon to crop rasters prior to modelling

The occurence points in the named list are used to create polygons. 
Notice that you can cluster the occ points using several clustering methods. 
See differences and choose one that fits your needs:
```r
occ.polys <- poly.c.batch(spp.occ.list)
occ.polys <- poly.c.batch(spp.occ.list, k=0, c.m="E")
occ.polys <- poly.c.batch(spp.occ.list, k=0, c.m="AP", q=.01) # less polygons
occ.polys <- poly.c.batch(spp.occ.list, k=0, c.m="AP", q=.1)
occ.polys <- poly.c.batch(spp.occ.list, k=0, c.m="AP", q=.2)
occ.polys <- poly.c.batch(spp.occ.list, k=0, c.m="AP", q=.3)
occ.polys <- poly.c.batch(spp.occ.list, k=0, c.m="AP", q=.5)
occ.polys <- poly.c.batch(spp.occ.list, k=0, c.m="AP", q=.8) # more polygons
occ.polys <- poly.c.batch(spp.occ.list, k=0, c.m="NB", method = "centroid", index = "duda")
occ.polys <- poly.c.batch(spp.occ.list, k=0, c.m="NB", method = "centroid", index = "ball")
occ.polys <- poly.c.batch(spp.occ.list, k=0, c.m="NB", method = "centroid", index = "pseudot2")
occ.polys <- poly.c.batch(spp.occ.list, k=0, c.m="NB", method = "centroid", index = "trcovw")
occ.polys <- poly.c.batch(spp.occ.list, k=0, c.m="NB", method = "centroid", index = "kl") 
occ.polys <- poly.c.batch(spp.occ.list, k=0, c.m="NB", method = "centroid", index = "sdindex") 
occ.polys <- poly.c.batch(spp.occ.list, k=0, c.m="NB", method = "centroid", index = "sdbw")
```

### - 1.2.1 creating buffer

... and the occurrence polygons are buffered using 1.5 degrees.
```r
occ.b <- bffr.batch(occ.polys, bffr.width = 1.5)
```

### - 1.3. Cut enviromental layers with occ.b and save in hardrive.
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

Finally, crop environmental variables for each species (and plot them for visual inspection)
```r
occ.b.env <- env.cut(occ.b, env.uncut)

for(i in 1:length(occ.b.env)){
  plot(occ.b.env[[i]][[1]])
  plot(occ.polys[[i]], border = "red", add = T)
  plot(occ.b[[i]], add = T)
}
```


## - 2. Prepare occurence data
### - 2.1 Filtering original dataset
Now we want to remove localities that are too close apart. We will do it for all species listed in "spp.occ.list".
```r
thinned.dataset.batch <- thin.batch(loc.data.lst = spp.occ.list)
```

### Great! Now we are ready for tunning species' ENMs


-----


## - 3. Tunning Maxent's feauture classes and regularization multiplier via ENMeval
### - 3.1 Load occurrence data (filtered localities). So, set working directory as correspond. 
After thinning, we choose one dataset for each species for modelling.
```r
occ.locs <- loadTocc(thinned.dataset.batch)
```

### - 3.3 model tuning using ENMeval
Here we will run ENMevaluate.batch to call ENMevaluate (from ENMeval package). By providing [at least] two lists, occurence and environmental data, we will be able to evaluate ENMs for as many species as listed in our occ.locs object. For details see ?ENMeval::ENMevaluate. Notice that you can use multiple cores for this task. This is specially usefull when there are a large number of models and species.
```r
ENMeval.res.lst <- ENMevaluate.batch(occ.locs, occ.b.env,method="block")
```

-----
## - 4. Model Fitting (Calibration)
### - 4.3 Run top corresponding models and save predictions 
#### - 4.3.1 save maxent best models and predictions for each model

Now, select maxent model calibrations and predictions using the function mxnt.cp.batch. This function can be run using a single core (default) or multiple cores available in a computer. There two ways of performing parallel processing: by species or by model. If the distribution of few species is being modelled, and models are computationally intensive, then processing by model will provide best results. If there are many species, probably parallel processing by species (split species across the multiple cores of a computer) will be faster.

```r

# Run model
mxnt.mdls.preds.lst <- mxnt.c.batch(ENMeval.o.l = ENMeval.res.lst, a.calib.l = occ.b.env, occ.l=occ.locs, wAICsum=0.99)


# Comparing single core processing and multiple core processing

# Parallel processing forking models (best for small sets of species)
system.time(
mxnt.mdls.preds.lst <- mxnt.c.batch(ENMeval.o.l = ENMeval.res.lst, a.calib.l = occ.b.env, occ.l=occ.locs, wAICsum=0.99, numCores=3, parallelTunning=TRUE)
)

# Parallel processing forking species (best for large sets of species)
system.time(
mxnt.mdls.preds.lst <- mxnt.c.batch(ENMeval.o.l = ENMeval.res.lst, a.calib.l = occ.b.env, occ.l=occ.locs, wAICsum=0.99, numCores=3, parallelTunning=FALSE)
)

```


## -  4.1 Projection
### - 4.1. Downloading environmental data

For projection it is necessary to download raster files with the environmnetal variables of interest. In this example, a directory called 'rasters' is created. Then, rasters from current and future climatic conditions projected for 2050 and 2070 are downloaded and loaded. Finally, two lists are created, one for current conditions and another for the two future cenarios.

```r

# Create directory to store raster files
dir.create("./rasters")

# Download data for present
current<-getData('worldclim', var='bio', res=10,path="rasters")

# Download data for future projection (2050)
futAC5085<-getData('CMIP5', var='bio', res=10, rcp=85, model='AC', year=50,path="rasters")

# Download data for future projection (2070)
futAC7085<-getData('CMIP5', var='bio', res=10, rcp=85, model='AC', year=70,path="rasters")


current.l<-list(current=current[[c(1,12,16,17,5,6,7,8)]])
future.l<-list(futAC5085=futAC5085[[c(1,12,16,17,5,6,7,8)]],futAC7085=futAC7085[[c(1,12,16,17,5,6,7,8)]])

names(future.l[[1]])<-names(env.uncut)
names(future.l[[2]])<-names(env.uncut)
```

# 4.1 Preparing projecion area: save rasters onto which the model will be projected in an object called "areas.projection"
# 4.1.1 select area for projection based on the extent of occ points

Now it is time to define the projection area for each species. The projection area can be the same for all species (in this example) of be defined individually. Here, the projection area will be defined as an square area slightly larger than the original occurrence of the species. Then, a two lists with models will be created for a species. In the first list, the projection will be performed using current climatic conditions. In the second list, two cenarios of futurure climate (defined above) are created.

```r
poly.projection <- pred.a.poly.batch(occ.polys, mult = .1, buffer=FALSE)#
plot(poly.projection[[1]])
plot(occ.polys[[1]], col="red", add=T)

pa.current.l <- pred.a.batch.mscn(poly.projection, current.l)
pa.future.l  <- pred.a.batch.mscn(poly.projection, future.l)
```

## 4.1.2 if the extent to project is the same for all species

When all species are to be projected using the same current and future climates and in the same region, then the following lines can be used to repeat the same lists of cenarios for all species (could be defined differently for each species if wanted)

```r
current.all <- lapply(mxnt.mdls.preds.lst,function(x)current.l)
future.all <- lapply(mxnt.mdls.preds.lst,function(x)future.l)
```

### 4.8 projections for present, future, and/or past

Finally, all species can be projected for all cenarios using all models. This is performed by `the mxnt.p.batch.mscn` function. The function has two arguments: 1) Model fit from maxent models within a list (see step 4.3 above) and 2) lists of rasters representing all cenarios to be projected. The last argument must be represented by a list of species. For each species, a list (within the first) is provided with all cenarios to be performed for each species (see last step above).

```r
# For single or multiple species
mxnt.mdls.preds.cf <- mxnt.p.batch.mscn(mxnt.mdls.preds.lst, a.proj.l = pa.current.l)
# using a single core (default)
mxnt.mdls.preds.cf <- mxnt.p.batch.mscn(mxnt.mdls.preds.cf, a.proj.l = pa.future.l)
# or using multiple cores
mxnt.mdls.preds.cf <- mxnt.p.batch.mscn(mxnt.mdls.preds.cf, a.proj.l = pa.future.l, numCores=2)

# plot projections
plot(mxnt.mdls.preds.cf$Bvarieg$mxnt.preds$current)
plot(mxnt.mdls.preds.cf$Bvarieg$mxnt.preds$futAC5085)
```

###  4.8.5 threshold for past and future predictions

We have the projections for each climatic scenario, now we must select one (or more) threshold criteria and 
apply on the projections.

```r
mods.thrshld.lst <- f.thr.batch(mxnt.mdls.preds.cf, thrshld.i = 4:5)
```

#### Plotting one projection for current climate and another for a future climatic scenario
```r
par(mfrow=c(1,2))
plot(mods.thrshld.lst$Bvarieg$current$binary$mtp[[1]])
plot(mods.thrshld.lst$Bvarieg$futAC5085$binary$mtp[[1]])
```

#### Plotting differences between current climate and future climatic scenarios for all thresholds we calculated
```r
f.plot.scn.diff(mxnt.mdls.preds.cf, mods.thrshld.lst)
```
