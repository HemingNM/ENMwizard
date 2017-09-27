ENMwizard
======================
This package has several functions to facilitate the use of Ecological Niche Modelling (ENM) in R

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

## ------- 1. Prepare environmental data

### ------- 1.1 Load occurence data

First, lets use occ data available in dismo package
```r
Bvarieg.occ <- read.table(paste(system.file(package="dismo"), "/ex/bradypus.csv", sep=""), header=TRUE, sep=",")
colnames(Bvarieg.occ) <- c("SPEC", "LONG", "LAT")
```

Now we make it a named list, where names correspond to species names
```r
spp.occ.list <- list(Bvarieg = Bvarieg.occ)
```

### ------- 1.2 create occ polygon to crop rasters prior to modelling

The occurence points in the named list are used to create polygons ...
```r
occ_polys <- f.poly.batch(spp.occ.list, o.path="occ_poly")
```

### ------- 1.2.1 creating buffer

... and the occurrence polygons are buffered 1.5 degrees wider.
```r
occ_b <- f.bffr(occ_polys, bffr.width = 1.5)
```

### ------- 1.3. Cut enviromental layers with occ_b and save in hardrive.
Specify the path to the environmental variables
```r
path.env <- "/path/to/variables/WorldClim/2_5min/bio_2-5m_bil"
biovars <- paste0("bio", 1:17) #c("bio5", "bio8", "bio10", "bio13", "bio16") # "bio18" tem problemas para o cerrado
pattern.env = 'asc'
path.env.out <- "3_envData"
```

Get uncut variables
```r
env_uncut <- list.files(path.env, pattern = pattern.env, full.names=TRUE)
env_uncut <- env_uncut[grepl(paste(paste0(biovars, ".", pattern.env), collapse = "|"), env_uncut)]
env_uncut <- stack(env_uncut) #predictors_uncut
crs(env_uncut) <- crs.set
```

If the variables were already cropped and saved in a raster brick, you just need to read them
```r
env_uncut <- brick(paste(path.env, "bio.grd", sep="/"))
```

Finally, crop environmental variables for each species (and plot them for visual inspection)
```r
occ_b_env <- f.cut.env(occ_b, env_uncut)

for(i in 1:length(occ_b_env)){
  plot(occ_b_env[[i]][[1]])
  plot(occ_polys[[i]], border = "red", add = T)
  plot(occ_b[[i]], add = T)
}
```
