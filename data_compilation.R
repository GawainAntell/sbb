library(factoextra)
list.files('data/')
comp <- data.frame('year' = 2010:1748)

# Questions about data ----------------------------------------------------

# Barr et al 2013 diatons -
# Several taxa are listed as "spp", e.g. "T. spp"
# But there are other T. species with separate data
# Should T spp be removed as partially redundant?

# Berger et al 2004 TOC -
# Should use non-detrended variable, yes?

# 250y composite df -------------------------------------------------------

# TODO
# write function that matches rows in data to rows of target composite df
# exclude any missing matches and round non-integer years
# return column(s) of target variables with NAs for non-measured study years

# ENSO
enso <- read.csv('data/Li-et-al-2011_cleaned.csv')
# what row in comp does each Li row's year correspond to?
whenEnso <- match(enso$year, comp$year)
tooDeep <- is.na(whenEnso) # data series extends older than study interval start
newCols <- tail(colnames(enso), -1)
comp[ , newCols] <- NA
comp[na.omit(whenEnso), newCols] <- enso[ !tooDeep, newCols]

# TOC
toc <- read.csv('data/Berger-et-al-2004_cleaned.csv')
# what row in comp does each Berger row's year correspond to?
whenToc <- match(toc$year, comp$year)
tooDeep <- is.na(whenToc) # data series extends older than study interval start
comp$toc <- NA
comp$toc[ na.omit(whenToc) ] <- toc$TOC[ !tooDeep ]

# PCA helper function:
# return positions of columns that PCA will be performed on -
# remove age data and any variables with zero variance
# optionally supply the name of other variables to exclude
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

# diatom abundances
# NB: there is one non-integer year (Barr13d$year[95], 1833.5) rounded off
Barr13d <- read.csv('data/Barron-et-al-2013-diatoms_cleaned.csv')
diaCols <- okCols(Barr13d, 'Other_planktic')
pcDia <- prcomp(Barr13d[, diaCols], scale = FALSE)
summary(pcDia)
# print(pcDia) # show loadings of individual species on each axis
yrDia <- get_pca_ind(pcDia)
# what row in comp does each Barr13d row's year correspond to?
whenDia <- match(round(Barr13d$year), comp$year)
comp$dia2 <- comp$dia1 <- NA
comp$dia1[ whenDia ] <- yrDia$coord[, 1]
comp$dia2[ whenDia ] <- yrDia$coord[, 2]

# silicoflagellate abundances
Barr13s <- read.csv('data/Barron-et-al-2013-silicoflagellates_cleaned.csv')
siliCols <- okCols(Barr13s)
pcSili <- prcomp(Barr13s[, siliCols], scale = FALSE)
summary(pcSili)
yrSili <- get_pca_ind(pcSili)
whenSili <- match(round(Barr13s$year), comp$year)
comp$sili1[ whenDia ] <- yrSili$coord[, 1]

# Stationarity tests ------------------------------------------------------

# example code copied from EcoRelease repo
# adfTest(globlRich) 
# acf(globlRich) 
# rich_AR <- arima(globlRich, order = c(1, 0, 0)) # fit AR1 model
# rich_e <- as.numeric(rich_AR$residuals)
# acf(rich_e)
# adfTest(rich_e) 
# rich_AR$coef
