library(WaveletComp)
source('utils.R') # for custom functions

df <- read.csv('data/composite-datasets-binned-250y.csv')
# for time axis in wavelet plots:
# must include a date-format column named 'date', and youngest record at top row
df$date <- as.Date(ISOdate(df$year, 1, 1))  # beginning of year
reorder <- order(df$year, decreasing = FALSE)
df <- df[ reorder, ]

oto <- read.csv('data/Jones-Checkley-2019_cleaned.csv')
oto$date <- as.Date(ISOdate(oto$year, 1, 1))  # beginning of year
reorderOto <- order(oto$year, decreasing = FALSE)
oto <- oto[ reorderOto, ]

scl <- read.csv('data/composite-datasets-binned-2ky.csv')
scl$date <- as.Date(ISOdate(scl$year, 1, 1))
reorderScl <- order(scl$year, decreasing = FALSE)
scl <- scl[ reorderScl, ]

# Univariate --------------------------------------------------------------

# TODO
# Wrap analyze.wavelet and wt.image into function w custom parameters
# Change color scale from power to original units e.g. log deposition rate

## ENSO variables ---------------------------------------------------------

# ENSO index

ensoW <- analyze.wavelet(df, 'ensoi',
                         loess.span = 0, # don't detrend
                         dt = 2, # 1 time unit (observation) = 2 years
                         dj = 1/250,
                         lowerPeriod = 4, 
                         upperPeriod = 64,
                         make.pval = TRUE, 
                         method = 'AR',
                         n.sim = 1000 # for significance tests
                         )
wt.image(ensoW, n.levels = 250,
         show.date = TRUE, date.format = "%Y-%m-%d",
         timelab = '', periodlab = '', main = 'ENSO index', 
         # timelab = 'years', periodlab = "period length (years)",
         legend.params = list(lab = "wavelet power levels")
         )
# color.key = "quantile" # default

# ENSO variability

ensoVarW <- analyze.wavelet(df, 'ensovar',
                         loess.span = 0, # don't detrend
                         dt = 2, 
                         dj = 1/250,
                         lowerPeriod = 4, 
                         upperPeriod = 64,
                         make.pval = TRUE, 
                         method = 'AR',
                         n.sim = 1000
)
wt.image(ensoVarW, n.levels = 250,
         show.date = TRUE, date.format = "%Y-%m-%d",
         timelab = '', periodlab = '', main = 'ENSO var', 
         legend.params = list(lab = "wavelet power levels")
         )

# compare ensovar reconstructed curve with PDO index

ensoVarRec <- reconstruct(ensoVarW, plot.waves = FALSE) 
            # lwd = c(1, 2)
            # legend.coords = "bottomleft"
            # ylim = c(-1.8, 1.8)
df$ensoVarRec <- ensoVarRec$series$ensovar.r

# NB recorded in metadata spreadsheet re: NOAA data:
# PDO has warming trend removed, which is why Barron et al 2013 p 9 suggest 
# it doesn't correlate as well w microfossils as SST
pdo <- read.table('data-raw/ERSST-v5-PDO-index.txt',
                  skip = 1, header = TRUE)
colnames(pdo)[1] <- 'year'
pdo$ann <- rowMeans(pdo[, -1]) # annual mean of monthly measurements

df$pdoAnn <- matchTime(dat = pdo, tmplt = df, yrCol = 'year', xtrctCol = 'ann')
sameTime <- !is.na(df$pdoAnn)
cor.test(df$ensoVarRec[sameTime], df$pdoAnn[sameTime])
 # cor.test(df$ensovar[sameTime], df$pdoAnn[sameTime])

# (may need to extract reconstruction values to plot manually)
# check whether ensovar data were already detrended/how filtered

# estimate all pairwise corelations
df[, -c(1:2, ncol(df))] |> cor() |> round(digits = 2)
# flag correlations with magnitude > 0.2
df[, -c(1:2, ncol(df))] |> cor() |> abs() > 0.2
# correlations from 1855 onwards, to include PDO
df[sameTime, -c(1:2)] |> cor() |> round(digits = 2)

# things to note:
# - use ensovar > ensovar reconstructed
# - DS1 perfectly captures same variances as sili1, and some of dia1
#   so use DS axes instead of dia and sili separately
# - TOC covaries with a lot
# - variation in enso >>> enso, and this independently of Sara's benthic results
#   also variation in enso similar to pdo but not an exact scaling

## TOC --------------------------------------------------------------------

toc <- analyze.wavelet(df, 'TOC',
                       loess.span = 0, # don't detrend
                       dt = 2, 
                       dj = 1/250,
                       lowerPeriod = 4,
                       upperPeriod = 64,
                       make.pval = TRUE, 
                       method = 'AR',
                       n.sim = 1000 # for significance tests. use 1000 for real
)
wt.image(toc, n.levels = 250,
         timelab = 'year', show.date = TRUE, date.format = "%Y-%m-%d",
         periodlab = "period length (years)",
         legend.params = list(lab = "wavelet power levels")
)

## Phytoplankton ----------------------------------------------------------

# DS1 and DS2

ds1 <- analyze.wavelet(df, 'DS1',
                            loess.span = 0,
                            dt = 2,
                            dj = 1/250,
                            lowerPeriod = 4,
                            upperPeriod = 64,
                            make.pval = TRUE, 
                            method = 'AR',
                            n.sim = 1000
)
wt.image(ds1, n.levels = 250,
         show.date = TRUE, date.format = "%Y-%m-%d",
         timelab = '', periodlab = '', main = 'PCA 1', 
         # timelab = 'years', periodlab = "period length (years)",
         legend.params = list(lab = "wavelet power levels")
)

ds2 <- analyze.wavelet(df, 'DS2',
                       loess.span = 0,
                       dt = 2,
                       dj = 1/250,
                       lowerPeriod = 4,
                       upperPeriod = 64,
                       make.pval = TRUE, 
                       method = 'AR',
                       n.sim = 1000
)
wt.image(ds2, n.levels = 250,
         show.date = TRUE, date.format = "%Y-%m-%d",
         timelab = '', periodlab = '', main = 'PCA 2',
         legend.params = list(lab = "wavelet power levels")
)

# Fish scales (2ky) -------------------------------------------------------

# NB: read in Baumgartner data afresh to lengthen the tseries
# (TOC data only goes back half as far)

sard <- analyze.wavelet(scl, 'sardine',
                        loess.span = 0,
                        dt = 10,
                        dj = 1/250,
                        lowerPeriod = 32,
                        upperPeriod = 512,
                        make.pval = TRUE, 
                        method = 'AR',
                        n.sim = 1000
)
wt.image(sard, n.levels = 250,
         show.date = TRUE, date.format = "%Y-%m-%d",
         timelab = '', periodlab = '', main = 'sardine',
         legend.params = list(lab = "wavelet power levels")
)

anch <- analyze.wavelet(scl, 'anchovy',
                        loess.span = 0,
                        dt = 10,
                        dj = 1/250,
                        lowerPeriod = 32,
                        upperPeriod = 512,
                        make.pval = TRUE, 
                        method = 'AR',
                        n.sim = 1000
)
wt.image(anch, n.levels = 250,
         show.date = TRUE, date.format = "%Y-%m-%d",
         timelab = '', periodlab = '', main = 'anchovy',
         legend.params = list(lab = "wavelet power levels")
)

# Fish otoliths (2ky) -----------------------------------------------------

# Myctophidae - laternfish, mesopelagic
# compare to wavelets published in Jones and Checkley
# (agreement is decent; note y axis reversed)
myct <- analyze.wavelet(oto, 'odr_myct',
                        loess.span = 0,
                        dt = 10,
                        dj = 1/250,
                        lowerPeriod = 32,
                        upperPeriod = 512,
                        make.pval = TRUE, 
                        method = 'AR',
                        n.sim = 100
)
wt.image(myct, n.levels = 250,
         show.date = TRUE, date.format = "%Y-%m-%d",
         timelab = '', periodlab = '',
         legend.params = list(lab = "wavelet power levels")
)

# Engraulidae
anch <- analyze.wavelet(oto, 'odr_engr',
                        loess.span = 0,
                        dt = 10,
                        dj = 1/250,
                        lowerPeriod = 32,
                        upperPeriod = 256,
                        make.pval = TRUE, 
                        method = 'AR',
                        n.sim = 1000
)
wt.image(anch, n.levels = 250,
         show.date = TRUE, date.format = "%Y-%m-%d",
         timelab = '', periodlab = '',
         legend.params = list(lab = "wavelet power levels")
)

# Merluccidae

hake <- analyze.wavelet(oto, 'odr_merl',
                        loess.span = 0,
                        dt = 10,
                        dj = 1/250,
                        lowerPeriod = 32,
                        upperPeriod = 256,
                        make.pval = TRUE, 
                        method = 'AR',
                        n.sim = 1000
)
wt.image(hake, n.levels = 250,
         show.date = TRUE, date.format = "%Y-%m-%d",
         timelab = '', periodlab = '',
         legend.params = list(lab = "wavelet power levels")
)

# Bivariate ---------------------------------------------------------------

# TODO export after 1000 simulations (VERY slow)

## 250y -------------------------------------------------------------------

# ensovar vs DS1

ensovar_ds1 <- analyze.coherency(df, my.pair = c('ensovar', 'DS1'),
                                 loess.span = 0,
                                 dt = 2, 
                                 dj = 1/250,
                                 lowerPeriod = 4,
                                 upperPeriod = 64,
                                 make.pval = TRUE, 
                                 method = 'AR',
                                 params = list(AR = list(p = 1)),
                                 n.sim = 50)
wc.image(ensovar_ds1, n.levels = 250,
         show.date = TRUE, date.format = "%Y-%m-%d",
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = '', periodlab = '', main = 'coherence - ensovar and pca 1')

# ensovar vs DS2

ensovar_ds2 <- analyze.coherency(df, my.pair = c('ensovar', 'DS2'),
                                 loess.span = 0,
                                 dt = 2, 
                                 dj = 1/250,
                                 lowerPeriod = 4,
                                 upperPeriod = 64,
                                 make.pval = TRUE, 
                                 method = 'AR',
                                 params = list(AR = list(p = 1)),
                                 n.sim = 50)
wc.image(ensovar_ds2, n.levels = 250,
         show.date = TRUE, date.format = "%Y-%m-%d",
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = '', periodlab = '', main = 'coherence - ensovar and pca 2')

# TOC vs DS1

toc_ds1 <- analyze.coherency(df, my.pair = c('TOC', 'DS1'),
                             loess.span = 0,
                             dt = 2, 
                             dj = 1/250,
                             lowerPeriod = 4,
                             upperPeriod = 64,
                             make.pval = TRUE, 
                             method = 'AR',
                             params = list(AR = list(p = 1)),
                             n.sim = 50)
wc.image(toc_ds1, n.levels = 250,
         show.date = TRUE, date.format = "%Y-%m-%d",
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = '', periodlab = '', main = 'coherence - toc and pca 1')

# TOC vs DS2

toc_ds2 <- analyze.coherency(df, my.pair = c('TOC', 'DS2'),
                             loess.span = 0,
                             dt = 2, 
                             dj = 1/250,
                             lowerPeriod = 4,
                             upperPeriod = 64,
                             make.pval = TRUE, 
                             method = 'AR',
                             params = list(AR = list(p = 1)),
                             n.sim = 50)
wc.image(toc_ds2, n.levels = 250,
         show.date = TRUE, date.format = "%Y-%m-%d",
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = '', periodlab = '', main = 'coherence - toc and pca 2')

## Common era -------------------------------------------------------------

# TOC vs sardines

toc_sard <- analyze.coherency(scl, my.pair = c('TOC', 'sardine'),
                              loess.span = 0,
                              dt = 10, 
                              dj = 1/100,
                              lowerPeriod = 16,
                              upperPeriod = 256,
                              make.pval = TRUE, 
                              method = 'AR',
                              params = list(AR = list(p = 1)),
                              window.type.t = 3,
                              window.type.s = 3, # boxcar
                              n.sim = 100)
wc.image(toc_sard, n.levels = 250,
         show.date = TRUE, date.format = "%Y-%m-%d",
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = '', periodlab = '', main = 'coherence - toc and sardines')

# TOC vs anchovies

toc_anch <- analyze.coherency(scl, my.pair = c('TOC', 'anchovy'),
                              loess.span = 0,
                              dt = 10, 
                              dj = 1/100,
                              lowerPeriod = 16,
                              upperPeriod = 256,
                              make.pval = TRUE, 
                              method = 'AR',
                              params = list(AR = list(p = 1)),
                              window.type.t = 3,
                              window.type.s = 3, # boxcar
                              n.sim = 100)
wc.image(toc_anch, n.levels = 250,
         show.date = TRUE, date.format = "%Y-%m-%d",
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = '', periodlab = '', main = 'coherence - toc and anchovies')
