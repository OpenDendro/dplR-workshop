# Simple Signal-Free 

## Introduction

This is an initial draft that demonstrates how to build a simple signal-free chronology using the `ssf` function. Although this function is used to build chronologies, it is sufficiently different in its philosophy that it warrants its own chapter. 

In most regards the `ssf` function follows the procedures laid out by [Melvin and Briffa's](https://crudata.uea.ac.uk/cru/people/melvin/Melvin2008/Melvin2008.pdf) classic 2008 paper in Dendrochronologia [@melvin2008]. This is not a treatise on the pros and cons of the signal-free approach but attempts to demonstrate how the `ssf` function in can be applied to ring-width data. This particular implementation of the signal-free approach is quite limited in its scope and users who want more control and options (such as RCS) should look to the CRUST program detailed in Melvin and Briffa [@melvin2014a,@melvin2014b] available on [GitHub](https://github.com/ClimaticResearchUnit/CRUST).

## Packages
Here are the packages we will use. 


```r
library(dplR)
library(tidyverse)
library(cowplot)
library(ggExtra)
library(PNWColors)
```

## Data

We will demonstrate `ssf` with Bristlecone pine ring widths from Campito Mountain which is included with `dplR`.


```r
data(ca533)
class(ca533)
```

```
## [1] "rwl"        "data.frame"
```

```r
dat <- ca533
```

Note the `class` of the object which is `rwl`. This means, simply, that the data are in a format that `dplR` understands. It has series in columns and years stored as `rownames`. Objects of class `rwl` can be summarized, plotted, used for detrending and so on.

Before we continue we will truncate the the sample depth to a minimum of five series. There is nothing magic about the sample depth but having five series is usually sufficient for having a robust chronology.


```r
sampDepth <- rowSums(!is.na(dat))
datTrunc <- dat[sampDepth > 4,]
```

We will now grab a few variables that will come in handy during this vigentte although they aren't needed for a simple run of `ssf`.


```r
yrs <- time(datTrunc)
medianSegLength <- floor(median(rwl.stats(datTrunc)$year))
sampDepth <- rowSums(!is.na(datTrunc))
normalizedSampleDepth <- sqrt(sampDepth-1)/sqrt(max(sampDepth-1))
```

We will also set up some color palettes to use in plotting.


```r
# get some color palettes from Jake L's PNWColors
seriesColors <- pnw_palette(name="Starfish",n=dim(datTrunc)[2])
divColors <- pnw_palette("Moth",5)
```

Throughout this document, we will make use of both `tidyverse` syntax and do much of the plotting in `ggplot`. E.g., here is a plot of all the ring widths that we will use.


```r
rawRW <- data.frame(yrs = yrs, datTrunc) %>% 
  pivot_longer(!yrs,names_to = "series", values_to = "msmt") 

ggplot(data=rawRW,mapping = aes(x=yrs,y=msmt,color=series)) +
  geom_line(alpha=0.5) +
  scale_color_manual(values = seriesColors) +
  labs(y="mm",x="Years",caption = "Raw Measurements") +
  theme_cowplot() +
  theme(legend.position = "none")
```

```
## Warning: Removed 11868 row(s) containing missing values (geom_path).
```

<img src="04-ssf_files/figure-html/first plot-1.png" width="672" />

NB `ggplot` is issuing a warning above telling us that there are quite a few `NA` values in the data. That is nothing to worry about -- it is merely saying that the series have different start and end dates.

## Making the Simple Signal-Free Chronology

Now we create the simple signal-free chronology. In its most direct use, the `ssf` function outputs a `crn` object that has a simple onboard plotting method (`crn.plot`) which we will show rather than make a `ggplot` of the same data. Don't worry, there is plenty of fancy plotting below.


```r
ssfCrn <- ssf(rwl = datTrunc)
```

```
## Data read. First iteration done.
## Iteration: 2 Median Abs Diff: 0.00415 (12.05% of threshold)
## Iteration: 3 Median Abs Diff: 0.00217 (23.08% of threshold)
## Iteration: 4 Median Abs Diff: 0.00291 (17.19% of threshold)
## Iteration: 5 Median Abs Diff: 0.00133 (37.68% of threshold)
## Iteration: 6 Median Abs Diff: 0.001 (50.02% of threshold)
## Iteration: 7 Median Abs Diff: 0.00076 (66.13% of threshold)
## Iteration: 8 Median Abs Diff: 0.00066 (75.66% of threshold)
## Iteration: 9 Median Abs Diff: 0.00055 (91.4% of threshold)
## Iteration: 10 Median Abs Diff: 0.00047 (105.5% of threshold)
## Simple Signal Free Chronology Complete
## ssf was called with these arguments
## Detrending method: AgeDepSpline
## nyrs: 
## pos.slope: FALSE
## maxIterations: 25
## madThreshold: 5e-04
```

```r
str(ssfCrn) # note class of output object
```

```
## Classes 'crn' and 'data.frame':	1007 obs. of  2 variables:
##  $ sfc       : num  1.219 1.487 0.836 0.664 0.57 ...
##  $ samp.depth: num  5 5 5 5 5 5 5 5 5 5 ...
```

```r
plot(ssfCrn,add.spline=TRUE,nyrs=50,
     crn.line.col=divColors[3],spline.line.col=divColors[1],
     crn.lwd=1.5,spline.lwd=2)
```

<img src="04-ssf_files/figure-html/do ssf-1.png" width="672" />

In the above, the signal-free chronology is shown with a 50-year smoothing spline added for visualization of low-frequency variability. The algorithm converged after 10 iterations.

## Changing the Detrending Method

The help file for `ssf` gives more information on possible arguments. See `?ssf` for details. By default the function uses an age-depended spline (`dplR` function `ads`) for detrending. But we can also use a cubic smoothing spline (`caps`). We will arbitrarily set the stiffness to the median segment length of the input data, which is 627. Whether this is a wise choice or not is debatable.


```r
ssfCrn2 <- ssf(rwl = datTrunc,method="Spline",nyrs=medianSegLength)
```

```
## Data read. First iteration done.
## Iteration: 2 Median Abs Diff: 0.00212 (23.57% of threshold)
## Iteration: 3 Median Abs Diff: 0.00185 (26.97% of threshold)
## Iteration: 4 Median Abs Diff: 0.0017 (29.36% of threshold)
## Iteration: 5 Median Abs Diff: 0.00153 (32.63% of threshold)
## Iteration: 6 Median Abs Diff: 0.00144 (34.65% of threshold)
## Iteration: 7 Median Abs Diff: 0.00138 (36.32% of threshold)
## Iteration: 8 Median Abs Diff: 0.00125 (39.92% of threshold)
## Iteration: 9 Median Abs Diff: 0.00121 (41.25% of threshold)
## Iteration: 10 Median Abs Diff: 0.00114 (43.97% of threshold)
## Iteration: 11 Median Abs Diff: 0.00106 (47.17% of threshold)
## Iteration: 12 Median Abs Diff: 0.00103 (48.5% of threshold)
## Iteration: 13 Median Abs Diff: 0.00098 (50.96% of threshold)
## Iteration: 14 Median Abs Diff: 0.00095 (52.75% of threshold)
## Iteration: 15 Median Abs Diff: 0.00089 (55.93% of threshold)
## Iteration: 16 Median Abs Diff: 0.00086 (58.06% of threshold)
## Iteration: 17 Median Abs Diff: 0.00081 (61.74% of threshold)
## Iteration: 18 Median Abs Diff: 0.00074 (67.72% of threshold)
## Iteration: 19 Median Abs Diff: 0.00067 (74.48% of threshold)
## Iteration: 20 Median Abs Diff: 0.00064 (77.98% of threshold)
## Iteration: 21 Median Abs Diff: 6e-04 (82.93% of threshold)
## Iteration: 22 Median Abs Diff: 0.00057 (88.15% of threshold)
## Iteration: 23 Median Abs Diff: 0.00053 (94.86% of threshold)
## Iteration: 24 Median Abs Diff: 0.00048 (104% of threshold)
## Simple Signal Free Chronology Complete
## ssf was called with these arguments
## Detrending method: Spline
## nyrs: 627
## maxIterations: 25
## madThreshold: 5e-04
```

```r
plot(ssfCrn2,add.spline=TRUE,nyrs=50,
     crn.line.col=divColors[3],spline.line.col=divColors[1],
     crn.lwd=1.5,spline.lwd=2)
```

<img src="04-ssf_files/figure-html/do ssf with caps-1.png" width="672" />

Note that the algorithm takes more iterations to converge but produces a qualitatively similar chronology.

## Walk Through
To demonstrate how the signal free process works we can redo the chronology, this time returning information on the process at each iteration. Note the change in the structure (`str`) of the output object.


```r
ssfCrn <- ssf(rwl = datTrunc,return.info = TRUE,verbose = FALSE)
str(ssfCrn)
```

```
## List of 11
##  $ infoList                :List of 5
##   ..$ method       : chr "AgeDepSpline"
##   ..$ nyrs         : NULL
##   ..$ pos.slope    : logi FALSE
##   ..$ maxIterations: num 25
##   ..$ madThreshold : num 5e-04
##  $ iter0Crn                :Classes 'crn' and 'data.frame':	1007 obs. of  2 variables:
##   ..$ std       : num [1:1007] 1.288 1.557 0.884 0.703 0.6 ...
##   ..$ samp.depth: num [1:1007] 5 5 5 5 5 5 5 5 5 5 ...
##  $ ssfCrn                  :Classes 'crn' and 'data.frame':	1007 obs. of  2 variables:
##   ..$ sfc       : num [1:1007] 1.219 1.487 0.836 0.664 0.57 ...
##   ..$ samp.depth: num [1:1007] 5 5 5 5 5 5 5 5 5 5 ...
##  $ sfRW_Array              : num [1:1007, 1:34, 1:10] NA NA NA NA NA NA NA NA NA NA ...
##  $ sfRWRescaled_Array      : num [1:1007, 1:34, 1:10] NA NA NA NA NA NA NA NA NA NA ...
##  $ sfRWRescaledCurves_Array: num [1:1007, 1:34, 1:10] NA NA NA NA NA NA NA NA NA NA ...
##  $ sfRWI_Array             : num [1:1007, 1:34, 1:10] NA NA NA NA NA NA NA NA NA NA ...
##  $ sfCrn_Mat               : num [1:1007, 1:10] 1.264 1.394 0.849 0.666 0.568 ...
##  $ hfCrn_Mat               : num [1:1007, 1:25] 0.313 0.443 -0.102 -0.285 -0.383 ...
##  $ hfCrnResids_Mat         : num [1:1007, 1:9] -0.02912 0.07009 -0.02075 -0.00941 -0.02138 ...
##  $ MAD_Vec                 : num [1:9] 0.00415 0.00217 0.00291 0.00133 0.001 ...
```

The `ssfCrn` object is now a list with information on all the iterations that algorithm has run through. We can save these to their own objects for easier access.


```r
sfRW_Array <- ssfCrn$sfRW_Array
sfRWRescaled_Array <- ssfCrn$sfRWRescaled_Array
sfRWRescaledCurves_Array <- ssfCrn$sfRWRescaledCurves_Array
sfRWI_Array <- ssfCrn$sfRWI_Array
sfCrn_Mat <- ssfCrn$sfCrn_Mat
hfCrn_Mat <- ssfCrn$hfCrn_Mat
hfCrnResids_Mat <- ssfCrn$hfCrnResids_Mat
MAD_Vec <- ssfCrn$MAD_Vec
```

### Step 1 (Iteration 0)

This is the initial, naive chronology at iteration zero.


```r
plot(ssfCrn$iter0Crn,add.spline=TRUE,nyrs=50,
     crn.line.col=divColors[4],spline.line.col=divColors[2],
     crn.lwd=1.5,spline.lwd=2)
```

<img src="04-ssf_files/figure-html/init crn-1.png" width="672" />

We use this as a starting point and begin the iterations.

### Steps 2 and 3 (Iteration 1)

The algorithm first creates signal-free measurements for each series. At the first iteration these are the measurements divided by the initial chronology.

These signal-free measurements are output in `sfRW_Array` which is an array of years by series, by iteration. The signal-free measurements are then rescaled to have the original mean of each series.


```r
sfRWRescaled <- data.frame(yrs = yrs, sfRWRescaled_Array[,,1]) %>% 
  pivot_longer(!yrs,names_to = "series", values_to = "msmt") 
ggplot(data=sfRWRescaled,mapping = aes(x=yrs,y=msmt,color=series)) +
  geom_line(alpha=0.75) +
  scale_color_manual(values = seriesColors) +
  labs(caption="Rescaled Signal Free Measurements at Iteration 1",
       y="mm", x="Years") +
  theme_cowplot() +
  theme(legend.position = "none")
```

```
## Warning: Removed 11868 row(s) containing missing values (geom_path).
```

<img src="04-ssf_files/figure-html/sf resc iter1-1.png" width="672" />

### Steps 4 and 5 (Iteration 1)

In these steps the algorithm first looks for any places in the rescaled ring width array that have a sample depth of one and replaces the individual signal-free measurements values with the original measurement values. Then the curve fitting is repeated giving the detrending curves for iteration one.


```r
sfRWRescaledCurves <- data.frame(yrs = yrs, sfRWRescaledCurves_Array[,,1]) %>% 
  pivot_longer(!yrs,names_to = "series", values_to = "msmt") 
ggplot(data=sfRWRescaledCurves,mapping = aes(x=yrs,y=msmt,color=series)) +
  geom_line(alpha=0.75,size=1) +
  scale_color_manual(values = seriesColors) +
  labs(caption="Detrending Curves at Iteration 1",
       y="mm",x="Years") +
  theme_cowplot() +
  theme(legend.position = "none")
```

```
## Warning: Removed 11868 row(s) containing missing values (geom_path).
```

<img src="04-ssf_files/figure-html/unnamed-chunk-3-1.png" width="672" />

### Step 6 and 7 (Iteration 1)

With the detrending curves above, we make the first signal-free ring-width indices by dividing the original measurements by the curves above.


```r
sfRWI <- data.frame(yrs = yrs, sfRWI_Array[,,1]) %>% 
  pivot_longer(!yrs,names_to = "series", values_to = "msmt") 
ggplot(data=sfRWI,mapping = aes(x=yrs,y=msmt,color=series)) +
  geom_line(alpha=0.5) +
  geom_hline(yintercept = 1,linetype="dashed") +
  scale_color_manual(values = seriesColors) +
  labs(caption="Signal Free Indices at Iteration 1",
       y="RWI",
       x="Years") +
  theme_cowplot() +
  theme(legend.position = "none")
```

```
## Warning: Removed 11868 row(s) containing missing values (geom_path).
```

<img src="04-ssf_files/figure-html/unnamed-chunk-4-1.png" width="672" />

With these ring width indices we can make the signal-free chronology for the first iteration which we will plot with a 50-year smoothing spline.


```r
sfCrn <- data.frame(yrs = yrs, 
                    msmt = sfCrn_Mat[,1], 
                    msmstSm= caps(sfCrn_Mat[,1],nyrs = 50)) %>%
  pivot_longer(!yrs,names_to = "series", values_to = "msmt") 
ggplot(data=sfCrn,mapping = aes(x=yrs,y=msmt,
                                color=series,
                                size=series,
                                alpha=series)) +
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line() +
  scale_color_manual(values = divColors[c(3,3)]) +
  scale_size_manual(values = c(1,0.5)) +
  scale_alpha_manual(values = c(1,0.5)) +
  labs(caption="Signal Free Chronology at Iteration 1",
       y="RWI",x="Years") +
  theme_cowplot() +
  theme(legend.position = "none")
```

<img src="04-ssf_files/figure-html/unnamed-chunk-5-1.png" width="672" />


### Step 8 (Begin Iterations)

Now we begin iterating through the steps above until the stopping criteria is met or until the maximum number of iterations is reached. Here we show how the chronology changes from iteration one to iteration two and how the high-pass filtering is used to calculate the median absolute difference (MAD) used as the stopping criteria.


```r
sfCrn <- data.frame(yrs = yrs, msmt = sfCrn_Mat[,1:2]) %>%
  rename(`Iteration 1` = 2,`Iteration 2` = 3) %>%
  pivot_longer(!yrs,names_to = "Iteration", values_to = "msmt")
p1 <- ggplot(data=sfCrn,
             mapping = aes(x=yrs,y=msmt,color=Iteration,
                           linetype=Iteration,
                           alpha=Iteration)) +
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line() +
  scale_color_manual(values = divColors[c(4,2)]) +
  scale_alpha_manual(values = c(0.75,0.75)) +
  labs(caption ="Signal Free Chronology at Iteration 1 and 2",
       y="RWI",x="Years") +
  lims(y=c(0,3)) +
  theme_cowplot() +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position=c(.1,.1),legend.title = element_blank())

sfCrnCompare <- data.frame(sfCrn_Mat[,1:2])
p2 <- ggplot(data=sfCrnCompare, mapping = aes(x=X1,y=X2)) + 
  geom_point(alpha=0.5,color=divColors[3]) +
  geom_abline(slope=1,intercept = 0, linetype="dashed") +
  labs(y=element_blank(), x= element_blank()) +
  coord_equal() + theme_cowplot() + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())
p2 <-  ggMarginal(p2, type="boxplot", size=6,
                  xparams = list(fill = divColors[4],
                                 color = "grey30"),
                  yparams = list(fill = divColors[2],
                                 color = "grey30"))

ggdraw(p1) +
  draw_plot(plot = p2, 
            x = .05, y = .5, 
            width = .5, height = .5)
```

<img src="04-ssf_files/figure-html/unnamed-chunk-6-1.png" width="672" />

These chronologies are quite similar but differ at higher RWI values.

The median absolute error between these two chronologies is calculated on the high frequency component which is calculated using a cubic smoothing spline with stiffness set to the median segment length.


```r
highFreqCrn <- data.frame(yrs = yrs, 
                          highFreqIter1 = hfCrn_Mat[,1],
                          highFreqIter2 = hfCrn_Mat[,2]) %>%
  rename(`Iteration 1` = 2,`Iteration 2` = 3) %>%
  pivot_longer(!yrs,names_to = "Iteration", 
               values_to = "msmt")
# plot the high freq chrons
p1 <- ggplot(data=highFreqCrn,
             mapping = aes(x=yrs,y=msmt,color=Iteration,
                           linetype=Iteration,
                           alpha=Iteration)) +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_line() +
  scale_color_manual(values = divColors[c(4,2)]) +
  scale_alpha_manual(values = c(0.75,0.75)) +
  labs(caption="High Freq Chronology at Iteration 1 and 2",
       y="RWI",x="Years") +
  theme_cowplot() +
  theme(legend.position = "top",legend.title = element_blank())

# plot the high freq chron difference (stored in hfCrnResids_Mat)
highFreqDiff <- data.frame(yrs=yrs, 
                           difference =  hfCrnResids_Mat[,1])

p2 <- ggplot(highFreqDiff,aes(x=yrs,y=difference)) +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_line(color=divColors[3]) +
  labs(y="RWI Difference",x="Years",
       caption=paste0("Iteration 1v2 MAD ", round(MAD_Vec[1],4))) +
  lims(y=c(-0.05,0.1)) +
  theme_cowplot() +
  theme(legend.position = "none",
        legend.title = element_blank())

plot_grid(p1, p2,nrow = 2)
```

<img src="04-ssf_files/figure-html/unnamed-chunk-7-1.png" width="672" />

In the lower plot it's clear that the difference in the chronologies is highest at the start and end. We calculate the median absolute difference using the normalized sample depth as:


```r
median(abs(hfCrn_Mat[,2]*normalizedSampleDepth - hfCrn_Mat[,1]*normalizedSampleDepth))
```

```
## [1] 0.004149
```

This is also returned in `MAD_Vec`.

Since this value is greater than the the stopping criteria of 5e-04, the algorithm continues to iterate until the difference is below the threshold or the maximum number of iterations is reached. Here is the difference plot between iterations nine and ten.


```r
highFreqDiff <- data.frame(yrs=yrs, 
                           difference =  hfCrnResids_Mat[,9])

ggplot(highFreqDiff,aes(x=yrs,y=difference)) +
  geom_hline(yintercept = 0,linetype="dashed") +
  geom_line(color=divColors[3]) +
  labs(y="RWI Difference",x="Years",
       caption=paste0("Iteration 9v10 MAD ", round(MAD_Vec[9],4))) +
  lims(y=c(-0.05,0.1)) +
  theme_cowplot() +
  theme(legend.position = "none",
        legend.title = element_blank())
```

<img src="04-ssf_files/figure-html/unnamed-chunk-9-1.png" width="672" />

### Iterations Visualized

We can see the progression of the output chronology relatively clearly by plotting the 50-year smoothing spline of the final chronology for each iteration.


```r
crnCols <- pnw_palette(name="Starfish",n=dim(sfCrn_Mat)[2])

# smooth chrons
sfCrnSm <- data.frame(yrs = yrs, 
                      msmt = apply(sfCrn_Mat,2,caps,nyrs=50)) %>%
  rename_with(.fn = seq, .cols = -1) %>%
  pivot_longer(!yrs,names_to = "Iteration", values_to = "msmt") %>%
  mutate(Iteration = as.numeric(Iteration))

p1 <- ggplot() +
  geom_hline(yintercept = 1,linetype="dashed") +
  geom_line(data=sfCrnSm,
            mapping = aes(x=yrs,y=msmt,color=factor(Iteration)),
            alpha=0.75,size=1) +
  scale_color_manual(values = crnCols) +
  labs(caption="Signal Free Chronology",
       x="Years",y="RWI") +
  theme_cowplot() +
  theme(legend.position="none")

MADdf <- data.frame(Iteration = 2:dim(sfCrn_Mat)[2],
                    MAD = MAD_Vec)

p2 <- ggplot(MADdf, aes(x = Iteration,y = MAD,color=factor(Iteration))) +
  geom_hline(yintercept = 5e-04,color="grey30",linetype="dotted") + 
  geom_line(color="grey30") +
  geom_point(size=2,alpha=0.9) +
  scale_color_manual(values = crnCols) +
  scale_x_continuous(breaks=seq(2,max(MADdf$Iteration),by=2),
                     limits = c(1,max(MADdf$Iteration)+1),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0,max(MADdf$MAD))) +
  theme_cowplot() +
  theme(legend.position = "none",
        plot.background = element_rect(fill = "gray95"))

ggdraw(p1) +
  draw_plot(plot = p2, 
            x = .15, y = .6, 
            width = .4, height = .4)
```

<img src="04-ssf_files/figure-html/unnamed-chunk-10-1.png" width="672" />

We can also use the `gganimate` package to combine some of these elements in a way that allows us to see the progression of the chronology through the fitting process. 




![](sfAnim.gif)<!-- -->

One thing to note with these data is that the upturns and downturns at the ends of the chronology are persistent even though they shrink vastly in the later iterations.

## Final Thoughts
This document shows how the algorithm progresses through the iterations in order to calculate the simple signal-free chronology. We should caution here that while there are excellent theoretical reasons to use the signal-free approach, the method is not a panacea for all applications. Detrending, after all, is a dark art and it is always up to the user to carefully evaluate the data and build chronologies with expectation that the final product is heavily dependent on the site and on the methods used.
