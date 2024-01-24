library(tabulizer) # requires JDK version of Java, and rJava package
library(fUnitRoots) # for timeseries stationarity tests

# TODO
# standardize data format
# test timeseries stationarity
# ordinate and inspect PCA

# Vet and format data -----------------------------------------------------

## Enviro vars -------------------------------------------------------------

# PP, SST, TOC, 13C, 15N

# TOC - Berger 2004

nmBerg04 <- 'data-raw/Berger-et-al-2004_Baumgartner-1992.pdf'
Berg04raw <- extract_tables(nmBerg04, pages = 21:40) |>
  lapply( function(l) apply(l, c(1,2), as.numeric) )
# all entries treated as character; must convert type before rbind on rows

Berg04 <- do.call(rbind, Berg04raw)
colnames(Berg04) <- c('year', 'TOC', 'TOCdetrend')

# TOC - Wang et al. 2017

# ENSO - Li et al. 2011

# floods/droughts - Sarno 2020

# biogenic opal - Barron 2013

## Primary producers -------------------------------------------------------

# diatoms, silicoflagellates

## Forams ------------------------------------------------------------------

# planktic

## Fish --------------------------------------------------------------------

# anchovy, sardine, mesopelagic fish
# hake - piscivores

# Baumgartner 1992 data included in Berger et al 2004 supplement
Baum92raw <- extract_tables(nmBerg04, pages = 1:5) |> 
  lapply( function(l) apply(l, c(1,2), as.numeric) )
Baum92 <- do.call(rbind, Baum92raw)
colnames(Baum92) <- c('year', 'sardine', 'anchovy')

# Jones and Checkley 2019
# otolith deposition rate for 5 families, and SST and PP proxies from JP Kennett

nmJC19 <- 'data-raw/Jones-Checkley-2019-decadal.csv'
JC19 <- read.csv(nmJC19)
colnames(JC19)[1] <- 'year' # rename from 'year_ad' to match other datasets

# last 4 columns are empty; identify columns with data
findDataCols <- function(v){
  any( !is.na(v) )
} 
keepCols <- apply(JC19, 2,  findDataCols)
JC19 <- JC19[, keepCols]

# missing values currently set as NaN
JC19 <- apply(JC19, 2, function(x) replace(x, is.nan(x), NA))

# Stationarity tests ------------------------------------------------------

# example code copied from EcoRelease repo
  # adfTest(globlRich) 
  # acf(globlRich) 
  # rich_AR <- arima(globlRich, order = c(1, 0, 0)) # fit AR1 model
  # rich_e <- as.numeric(rich_AR$residuals)
  # acf(rich_e)
  # adfTest(rich_e) 
  # rich_AR$coef
