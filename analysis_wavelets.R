library(WaveletComp)
source('utils.R') # for custom matchTime() function

df <- read.csv('data/composite-datasets-binned-250y.csv')

# Univariate --------------------------------------------------------------

# ENSO index

ensoW <- analyze.wavelet(df, 'ensoi',
                         loess.span = 0, # don't detrend
                         dt = 1, # time unit = 1 obs/unit
                         dj = 1/250,
                         lowerPeriod = 2, # = 4 years
                         upperPeriod = 64,
                         make.pval = TRUE, 
                         n.sim = 100 # for significance tests. use 1000 for real
                         )
wt.image(ensoW, n.levels = 250,
         legend.params = list(lab = "wavelet power levels")
         )
# color.key = "quantile" # default

# ENSO variability

ensoVarW <- analyze.wavelet(df, 'ensovar',
                         loess.span = 0, # don't detrend
                         dt = 1, # time unit = 1 obs/unit
                         dj = 1/250,
                         lowerPeriod = 2, # = 4 years
                         upperPeriod = 64,
                         make.pval = TRUE, 
                         n.sim = 100 # for significance tests. use 1000 for real
)
wt.image(ensoVarW, n.levels = 250,
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

# Bivariate ---------------------------------------------------------------

# ensoi vs DS1
# ensoi vs DS2
# TOC vs DS1
# TOC vs DS2
