### 4.8 predictions for present, future, and/or past
#### 4.8.b Using multiple cores (parallel processing)
#### Create projection with 2 species (same species repeated)

The same process can be performed using parallel processing. This can dramatially increase processing speed when multiple species are modelled. To exemplify, the polygon with the projection area of a single species is duplicated using the append function. Then, the list of climatic conditions for projection is duplicated. Finally, the maxent model fitted for a single species above is also duplicated.

```r
poly.projection.multi <- append(poly.projection, poly.projection)

pa.current.l.multi <- pred.a.batch.mscn(poly.projection.multi, pa.current.l)

mxnt.mdls.preds.lst.multi <- append(mxnt.mdls.preds.lst, mxnt.mdls.preds.lst)
```

We can run the model for both species sequentially (without parallel processing) or simultaneously (parallel processing) and compare the amount of time required to run the model.

#### Run without parallel processing
```r
system.time(
  mxnt.mdls.preds.cf2 <- mxnt.p.batch.mscn(mxnt.mdls.preds.lst.multi, a.proj.l = pa.current.l.multi)
)
```

#### Run with parallel processing

ps. progress bars are not shown

```r
# Using single core (default)
system.time(
  mxnt.mdls.preds.cf2 <- mxnt.p.batch.mscn(mxnt.mdls.preds.lst.multi, a.proj.l = pa.current.l.multi)
)

# Sending each model to one core (best if single or few species with complex models/heavy rasters)
system.time(
  mxnt.mdls.preds.cf2 <- mxnt.p.batch.mscn(mxnt.mdls.preds.lst.multi, a.proj.l = pa.current.l.multi, numCores=2)
)

# Sending each species to one core (best if multiple species)
system.time(
  mxnt.mdls.preds.cf2 <- mxnt.p.batch.mscn(mxnt.mdls.preds.lst.multi, a.proj.l = pa.current.l.multi, numCores=2,parallelTunning=FALSE)
)

```

