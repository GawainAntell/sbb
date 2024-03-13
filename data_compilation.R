library(timeSeries) # library(wql) # linear interpolation within time-series
library(factoextra) # for PCA extraction
library(fUnitRoots) # for timeseries stationarity tests
library(ggplot2)
library(cowplot)
library(RColorBrewer)
source('utils.R') # for custom helper functions
# csv names of raw data files to concatenate:
allFl <- list.files('data/')
allFl[ grep('cleaned', allFl) ]

# Questions about data ----------------------------------------------------

# Barr et al 2013 diatons -
# Several taxa are listed as "spp", e.g. "T. spp"
# But there are other T. species with separate data
# Should T spp be removed as partially redundant?

# Berger et al 2004 TOC -
# Should use non-detrended variable, yes?

# 250y composite df -------------------------------------------------------

comp <- data.frame('year' = 2010:1748)
# template for composite dataframe over the years of the shorter study span

# PCA helper function:
# - return positions of columns that PCA will be performed on -
# - remove age data and any variables with zero variance
# - optionally supply the name of other variables to exclude
okCols <- function(dat, otherExclude = NULL){
  yrCol <- which(colnames(dat) == 'year')
  constCol <- which( apply(dat, 2, sd) == 0 ) # static value; can't use in PCA
  if (! is.null(otherExclude)){
    otherNo <- which(colnames(dat) %in% otherExclude)
    out <- setdiff(1:ncol(dat), c(yrCol, constCol, otherNo)) 
  } else {
    out <- setdiff(1:ncol(dat), c(yrCol, constCol)) 
  }
  out 
}

# ENSO
enso <- read.csv('data/Li-et-al-2011_cleaned.csv')
# matchTime function defined in utils.R
enso2add <- matchTime(dat = enso, tmplt = comp, 
                      xtrctCol = c('ensoi', 'ensovar'))
comp <- cbind(comp, enso2add)

# TOC
toc <- read.csv('data/Berger-et-al-2004_cleaned.csv')
toc2add <- matchTime(dat = toc, tmplt = comp, 
                     xtrctCol = c('TOC', 'TOCdetrend') )
comp <- cbind(comp, toc2add)

# biogenic silica
opal <- read.csv('data/Barron-et-al-2013-opal_cleaned.csv')
opal2add <- matchTime(dat = opal, tmplt = comp, xtrctCol = 'Biogenic_silica')
colnames(opal2add) <- 'opal'
comp <- cbind(comp, opal2add)

# diatom abundances
# NB: there is one non-integer year (Barr13d$year[95] = 1833.5) rounded off
Barr13d <- read.csv('data/Barron-et-al-2013-diatoms_cleaned.csv')
diaCols <- okCols(Barr13d, 'Other_planktic')
pcDia <- prcomp(Barr13d[, diaCols], scale = FALSE) 
# print(pcDia) # show loadings of individual species on each axis
smryDia <- summary(pcDia)
smryDia$importance['Cumulative Proportion', 1:4]
# keep first 2 axes (explain > 66% variation)
#     PC1     PC2     PC3     PC4 
# 0.37505 0.66217 0.79005 0.84585 
yrDia <- get_pca_ind(pcDia)
yrDiaCoords <- yrDia$coord
yrDiaCoords <- cbind(yrDiaCoords, 'year' = round(Barr13d$year))
dia2add <- matchTime(dat = yrDiaCoords, tmplt = comp, 
                     xtrctCol = c('Dim.1', 'Dim.2'))
colnames(dia2add) <- c('dia1', 'dia2')
comp <- cbind(comp, dia2add)

# silicoflagellate abundances
# NB: there is one non-integer year (Barr13s$year[95] = 1833.5) rounded off
Barr13s <- read.csv('data/Barron-et-al-2013-silicoflagellates_cleaned.csv')
siliCols <- okCols(Barr13s)
pcSili <- prcomp(Barr13s[, siliCols], scale = FALSE)
smrySili <- summary(pcSili)
smrySili$importance['Cumulative Proportion', 1:4]
# keep first axis (explains ~60% variation)
#     PC1     PC2     PC3     PC4 
# 0.59709 0.84951 0.89828 0.94291
yrSili <- get_pca_ind(pcSili)
yrSiliCoords <- yrSili$coord
yrSiliCoords <- cbind(yrSiliCoords, 'year' = round(Barr13s$year))
sili2add <- matchTime(dat = yrSiliCoords, tmplt = comp, xtrctCol = 'Dim.1')
colnames(sili2add) <- 'sili1'
comp <- cbind(comp, sili2add)

# option: combine diatoms and silicoflagellates, take first 2 axes:
BarrBoth <- cbind(Barr13d[, diaCols], Barr13s[, siliCols])
pcDS <- prcomp(BarrBoth, scale = FALSE) 
smryDS <- summary(pcDS)
smryDS$importance['Cumulative Proportion', 1:4]
#     PC1     PC2     PC3     PC4 
# 0.43557 0.62612 0.72844 0.81071
yrDS <- get_pca_ind(pcDS)
yrDScoords <- yrDS$coord
yrDScoords <- cbind(yrDScoords, 'year' = round(Barr13s$year))
DS2add <- matchTime(dat = yrDScoords, tmplt = comp, 
                    xtrctCol = c('Dim.1', 'Dim.2'))
colnames(DS2add) <- c('DS1', 'DS2')
comp <- cbind(comp, DS2add)

# 2ky composite df --------------------------------------------------------

# combine fish scale and TOC data over the common era at 10y resolution
# TOC data from Berger (above) extends back over much of the longer interval of 
# Baumgartner, although the fish scales are decadal and only continue to 1970

scl <- read.csv('data/Baumgartner-1992_cleaned.csv')
toc2add2k <- matchTime(dat = toc, tmplt = scl, 
                       xtrctCol = c('TOC', 'TOCdetrend') )
ce <- cbind(scl, toc2add2k)
ce <- ce[ complete.cases(ce), ]
# write.csv(ce, 'data/composite-datasets-binned-2ky.csv', row.names = FALSE)

# Interpolation/downscaling -----------------------------------------------

# downscale (linear interpolation) decadal data to 2yr resolution of other data

sclTimeTmplt <- data.frame('year' = 2010:1740)
# need to include the 1740 data to interpolate through to 1750 annually
sclAnnOrig <- matchTime(dat = scl, tmplt = sclTimeTmplt, 
                        xtrctCol = c('sardine', 'anchovy'))
dateChar <- as.Date(ISOdate(sclTimeTmplt$year, 1, 1)) |> as.character()
sclTs <- timeSeries(sclAnnOrig, charvec = dateChar)
sclAnnInt <- na.omit(sclTs, method = 'ir', interp = 'linear') |> 
  as.data.frame()
sclAnnInt$year <- rownames(sclAnnInt) |> as.Date() |> format('%Y')

scl2add250y <- matchTime(dat = sclAnnInt, tmplt = comp, 
                         xtrctCol = c('sardine', 'anchovy'))
comp <- cbind(comp, scl2add250y)

# Biannual time binning ---------------------------------------------------

top <- 1970 # last yr of Baumgartner fish data; 1987 for TOC; 2007 for plankton
bins <- seq(from = top, to = 1748, by = -2) # defined by BOTTOM year
binnd <- data.frame()
for (b in bins){
  bRows <- which(comp$year %in% c(b, b+1))
  naTally <- apply(comp[bRows, ], 1, function(x){
    is.na(x) |> sum()
  })
  bestRow <- which.min(naTally) |> names() |> as.numeric()
  copy <- cbind('yearBin' = b, comp[bestRow, ])
  binnd <- rbind(binnd, copy)
}
# one manual correction necessary - 
# fish scales sampled in 1970, everything else in 1971
# add 1970 fish values to 1970-1971 time bin where NA
binnd[1, c('sardine', 'anchovy')] <- 
  comp[comp$year == 1970, c('sardine', 'anchovy')]
if (is.na(binnd) |> sum() > 0){ stop('check the data for holes') }

# write.csv(binnd, 'data/composite-datasets-binned-250y.csv', row.names = FALSE)

# Stationarity tests ------------------------------------------------------

for (col in 3:ncol(binnd)){
  v <- binnd[, col]
  tst <- adfTest(v) # testing null that there is NON-stationarity
  p <- signif(tst@test$p.value, 2)
  paste(colnames(binnd)[col], 'adf test p =', p) |> print()
  # inspect autocorr at all lags and magnitude of 1st-order autocorr
    # acf(v) 
    # vAR <- arima(v, order = c(1, 0 , 0))
    # vAR$residuals |> as.numeric() |> acf()
    # vAR$coef
}
# [1] "ensoi adf test p = 0.01"
# [1] "ensovar adf test p = 0.5"
# [1] "TOC adf test p = 0.47"
# [1] "TOCdetrend adf test p = 0.28"
# [1] "opal adf test p = 0.38"
# [1] "dia1 adf test p = 0.01"
# [1] "dia2 adf test p = 0.01"
# [1] "sili1 adf test p = 0.01"
# [1] "DS1 adf test p = 0.01"
# [1] "DS2 adf test p = 0.01"
# [1] "sardine adf test p = 0.01"
# [1] "anchovy adf test p = 0.1"

# 250y time series of fish scales, sans interpolated data
# (confirm that stationarity test values not inflated from pseudoreplicates)
sard250 <- na.omit(scl2add250y$sardine)
adfTest(sard250)
anch250 <- na.omit(scl2add250y$anchovy)
adfTest(anch250)

# Detrending --------------------------------------------------------------

trendVars <- c('ensovar', 'TOC', 'TOCdetrend', 'opal', 'anchovy')
# indicated as non-stationary by tests in previous section

# setup data structures to fill in from detrended data
tblDetrend <- data.frame(matrix(ncol = 6, nrow = length(trendVars)))
colnames(tblDetrend) <- c('Variable', 'LM intercept', 
                          'LM coefficient', 'Coefficient SE', 
                          'ADF test statistic', 'ADF p-val')
rownames(tblDetrend) <- trendVars
detrendVars <- paste0(trendVars, 'Detrend')
binnd[, detrendVars] <- NA

for (v in trendVars){
  fmla <- formula(paste(v, '~ yearBin'))
  timeLm <- lm(fmla, data = binnd)
  # save model coefficients into table for supplemental information
  timeLmCoef <- summary(timeLm)$coefficients
  tblDetrend[v, 'LM intercept'] <- timeLmCoef['(Intercept)', 'Estimate'] |>
    signif(digits = 3)
  tblDetrend[v, c('LM coefficient', 'Coefficient SE')] <- 
    timeLmCoef['yearBin', c('Estimate', 'Std. Error')] |> 
    signif(digits = 2)
  # also recalculate autocorrelation coefficient
  tstVals <- adfTest(binnd[, v]) # testing null that there is NON-stationarity
  p    <- tstVals@test$p.value   |> signif(digits = 2)
  stat <- tstVals@test$statistic |> signif(digits = 2)
  tblDetrend[v, c('ADF test statistic', 'ADF p-val')] <- c(stat, p)
  # add detrended values (residuals) to composite dataframe for future analysis
  newVarNm <- paste0(v, 'Detrend')
  binnd[, newVarNm] <- timeLm$residuals
}

# Tseries plots -----------------------------------------------------------

plotSeries <- function(vNm, trendVars, decadalVars, dat){
  v <- dat[, vNm]
  plotV <- ggplot(dat, aes(x = yearBin, y = v)) +
    theme_bw() +
    scale_x_continuous(name = element_blank(),
                       expand = c(0, 0) # ,limits = c(1748, 1970)
                       ) +
    geom_line() +
    scale_y_continuous(name = vNm)
  
  # overplot data points, on original sampling scale
  if (vNm %in% decadalVars){
    decades <- dat[ , 'yearBin'] %in% seq(1970, 1750, by = -10)
    dat10y <- dat[decades, ]
    v10y <- dat10y[ , vNm]
    plotV <- plotV +
      geom_point(data = dat10y, aes(x = yearBin, y = v10y))
  } else {
    # for data sampled at biannual resolution or finer
    plotV <- plotV +
      geom_point()
  }

  # overplot trend line if applicable
  if (vNm %in% trendVars){
    plotV <- plotV +
      geom_smooth(method = lm, aes(x = yearBin, y = v)) 
    # default CI = 95%; specify otherwise with 'level'
  }
  plotV
}

seriesVars <- c('ensovar', 'TOC', 'opal', 'DS1', 'DS2', 'sardine', 'anchovy')
seriesPlots <- lapply(seriesVars, plotSeries, 
       trendVars = trendVars, 
       decadalVars = c('sardine', 'anchovy'),
       dat = binnd)

plot_grid(seriesPlots[[1]], seriesPlots[[2]], 
          seriesPlots[[3]], seriesPlots[[4]], 
          seriesPlots[[5]], seriesPlots[[6]], 
          seriesPlots[[7]], ncol = 2,
          labels = "AUTO", align = 'hv')
