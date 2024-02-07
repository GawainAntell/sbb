library(factoextra)
list.files('data/')

# Questions about data ----------------------------------------------------

# Barr et al 2013 diatons -
# Several taxa are listed as "spp", e.g. "T. spp"
# But there are other T. species with separate data
# Should T spp be removed as partially redundant?

# Berger et al 2004 TOC -
# Should use non-detrended variable, yes?

# 250y composite df -------------------------------------------------------

comp <- data.frame('year' = 2010:1748)

# Helper function - match rows in data to rows of target composite df:
# - exclude any missing matches and round non-integer years
# - return column(s) of target variables with NAs for non-measured study years
# tmplt = vector/column of years from target composite dataframe
# yrCol = name or position of column in dat AND tmplt with age/year data
# xtrctCol = name(s) of column(s) in dat to add to composite df
matchTime <- function(dat, tmplt = comp, yrCol = 'year', xtrctCol){
  whenDat <- match(dat[, yrCol], tmplt[, yrCol])
  outsideTmplt <- is.na(whenDat) # if data extends beyond target study interval
  newCols <- rep(NA, nrow(tmplt)) |> data.frame()
  reps <- length(xtrctCol) - 1
  if (reps > 0){
    for (i in 1:reps){
      newCols <- cbind(newCols, newCols)
    }
  }
  newCols[ na.omit(whenDat), ] <- dat[ !outsideTmplt, xtrctCol]
  colnames(newCols) <- xtrctCol
  newCols
}

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
enso2add <- matchTime(dat = enso, xtrctCol = c('ensoi', 'ensovar'))
comp <- cbind(comp, enso2add)

# TOC
toc <- read.csv('data/Berger-et-al-2004_cleaned.csv')
toc2add <- matchTime(dat = toc, xtrctCol = 'TOC') # c('TOC', 'TOCdetrend')
comp <- cbind(comp, toc2add)

# biogenic silica
opal <- read.csv('data/Barron-et-al-2013-opal_cleaned.csv')
opal2add <- matchTime(dat = opal, xtrctCol = 'Biogenic_silica')
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
dia2add <- matchTime(dat = yrDiaCoords, xtrctCol = c('Dim.1', 'Dim.2'))
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
sili2add <- matchTime(dat = yrSiliCoords, xtrctCol = 'Dim.1')
colnames(sili2add) <- 'sili1'
comp <- cbind(comp, sili2add)

# Stationarity tests ------------------------------------------------------

# example code copied from EcoRelease repo
# adfTest(globlRich) 
# acf(globlRich) 
# rich_AR <- arima(globlRich, order = c(1, 0, 0)) # fit AR1 model
# rich_e <- as.numeric(rich_AR$residuals)
# acf(rich_e)
# adfTest(rich_e) 
# rich_AR$coef
