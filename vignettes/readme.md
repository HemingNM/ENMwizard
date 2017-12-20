ENMwizard
======================
### Advanced Tecniques for Ecological Niche Modeling Made Easy

This package provides tools to facilitate the use of advanced techiques related to Ecological Niche Modeling (ENM) and the automation of repetitive tasks (when modeling several species). This package has functions for easier: 1. preparation of occurence and environmental data; 2. model tunning (thanks to the package ENMeval); 3. model fitting and projection. ENMwizard also implements methods described in Gutiérrez & Heming in prep.

-----

# Installation
ENMwizard is downloadable from https://github.com/HemingNM/ENMwizard. You can download it using devtools to install from GitHub.

### Installing from GitHub using devtools
Run the following code from your R console:


```r
install.packages("devtools")
library(devtools)
install_github("HemingNM/ENMwizard")
library(ENMwizard)
library(raster)
```



### Install from zip file
You can also download a zip file containing the package and install it from R.

Download from https://github.com/HemingNM/ENMwizard/archive/master.zip and run the following code (where PATH is the path to the zip file)

```r
install.packages("devtools")
library(devtools)
install_local("PATH")
library(ENMwizard)
library(raster)
```



-----


# Using ENMwizard's magic wand

## ------- 1. Prepare environmental data

### ------- 1.1 Load occurence data

First, lets use occ data available in dismo package
```r
Bvarieg.occ <- read.table(paste(system.file(package="dismo"),
"/ex/bradypus.csv", sep=""), header=TRUE, sep=",")

head(Bvarieg.occ)# Check first rows

# Column names must be in capital letters
colnames(Bvarieg.occ) <- c("SPEC", "LONG", "LAT") # Change column names
```

Now we make it a named list, where names correspond to species names
```r
spp.occ.list <- list(Bvarieg = Bvarieg.occ)
```

### ------- 1.2 create occ polygon to crop rasters prior to modelling

The occurence points in the named list are used to create polygons ...
```r
occ.polys <- poly.c.batch(spp.occ.list)

```

### ------- 1.2.1 creating buffer

... and the occurrence polygons are buffered 1.5 degrees wider.
```r
occ.b <- bffr.batch(occ.polys, bffr.width = 1.5)
```

### ------- 1.3. Cut enviromental layers with occ.b and save in hardrive.
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


## ------- 2. Prepare occurence data
### ------- 2.1 Filtering original dataset
Now we want to remove localities that are too close apart. We will do it for all species listed in "spp.occ.list".
```r
thinned.dataset.batch <- thin.batch(loc.data.lst = spp.occ.list)
```

### Great! Now we are ready for tunning species' ENMs


-----


## ------- 3. Tunning Maxent's feauture classes and regularization multiplier via ENMeval
### ------- 3.1 Load occurrence data (filtered localities). So, set working directory as correspond. 
After thinning, we choose one dataset for each species for modelling.
```r
occ.locs <- loadTocc(thinned.dataset.batch)
```

### ------- 3.3 model tuning using ENMeval
Here we will run ENMevaluate.batch to call ENMevaluate (from ENMeval package). By providing [at least] two lists, occurence and environmental data, we will be able to evaluate ENMs for as many species as listed in our occ.locs object. For details see ?ENMeval::ENMevaluate. Notice that you can use multiple cores for this task. This is specially usefull when there are a large number of models and species.
```r
ENMeval.res.lst <- ENMevaluate.batch(occ.locs, occ.b.env,method="block")
```

-----
### TODO
# 4. Model Fitting (Calibration) and Projection
# 4.1 Preparing projecion area: save rasters onto which the model will be projected in an object called "areas.projection"
# 4.1.1 select area for projection based on the extent of occ points
```r
area.projection <- pred.a.poly.batch(occ.polys,deg.incr=2, env.uncut, mult = .75, buffer=FALSE)#
plot(area.projection[[1]])
plot(occ.polys[[1]], col="red", add=T)
```

#### 4.3 Run top corresponding models and save predictions 
#### 4.3.1 save maxent best models and predictions for each model
```r
mxnt.mdls.preds.lst <- mxnt.cp.batch(ENMeval.res = ENMeval.res.lst,a.calib.l = occ.b.env, occ.l=occ.locs, wAICsum=0.99,numCores=3)
```


