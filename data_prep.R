library(stringr)
library(tabulizer) # read PDF tables; requires JDK vsn of Java, and rJava pkg
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

# diatoms - Barron et al. 2013

Barr13d <- read.csv('data-raw/Barron-et-al-2013-diatoms.csv')
# spaces converted to periods from excel, and periods appear in gen abbreviations
colnames(Barr13d) <- gsub(x = colnames(Barr13d), pattern = '\\.\\.', replacement = '_')  
colnames(Barr13d) <- gsub(x = colnames(Barr13d), pattern = '\\.', replacement = '_')  
# omit trailing space where present
trails <- str_ends(colnames(Barr13d), '_')
colnames(Barr13d)[trails] <- colnames(Barr13d)[trails] |>
  str_sub(end = -2)
# some columns contain metadata, not species abundances (or not relevant taxa)
omitCols <- c('Varve_range', 'Bottom_varve',
              'Benthic', 'Freshwater_planktic', 'Reworked', 'Total_counted')
Barr13d <- Barr13d[ , !colnames(Barr13d) %in% omitCols]
# final rows are completely empty
findBlanks <- function(v){
  any( !is.na(v) )
} 
keepRowsD <- apply(Barr13d, 1,  findBlanks)
Barr13d <- Barr13d[keepRowsD, ]
colnames(Barr13d)[1] <- 'year' # this is the oldest (bottom) year of the sample

# silicoflagellates - Barron et al. 2013

Barr13raw <- read.csv('data-raw/Barron-et-al-2013-silicoflagellates.csv')
# entire dataset is split at 1884, and second half pasted as columns to right :|
dupeCols <- str_ends(colnames(Barr13raw), '1')
topHlf <- Barr13raw[, !dupeCols]
topHlf<- topHlf[ , -ncol(topHlf)] # blank divider column
lwrHlf <- Barr13raw[, dupeCols]
colnames(lwrHlf) <- str_sub(colnames(lwrHlf), end = -3) # strip '.1' suffix
Barr13s <- rbind(topHlf, lwrHlf)
# repeat cleaning steps for Barron et al. 2013 diatom data (above)
colnames(Barr13s) <- gsub(x = colnames(Barr13s), pattern = '\\.\\.', replacement = '_')
colnames(Barr13s) <- gsub(x = colnames(Barr13s), pattern = '\\.', replacement = '_')
trailsS <- str_ends(colnames(Barr13s), '_')
colnames(Barr13s)[trailsS] <- colnames(Barr13s)[trailsS] |>
  str_sub(end = -2)
Barr13s <- Barr13s[ , !colnames(Barr13s) %in% omitCols]
keepRowsS <- apply(Barr13s, 1,  findBlanks) # do this after removing metadata cols
Barr13s <- Barr13s[keepRowsS, ]
colnames(Barr13s)[1] <- 'year' # this is the oldest (bottom) year of the sample

# consider converting count data type from integer to numeric

## Forams ------------------------------------------------------------------

# planktic

## Fish --------------------------------------------------------------------

# anchovy, sardine - Baumgartner 1992
# data included in Berger et al 2004 supplement
Baum92raw <- extract_tables(nmBerg04, pages = 1:5) |> 
  lapply( function(l) apply(l, c(1,2), as.numeric) )
Baum92 <- do.call(rbind, Baum92raw)
colnames(Baum92) <- c('year', 'sardine', 'anchovy')

# mesopelagic fish - Jones and Checkley 2019
# otolith deposition rate for 5 families, and SST and PP proxies from JP Kennett

nmJC19 <- 'data-raw/Jones-Checkley-2019-decadal.csv'
JC19 <- read.csv(nmJC19) 
colnames(JC19)[1] <- 'year' # rename from 'year_ad' to match other datasets

# last 4 columns are empty; identify columns with data
keepCols <- apply(JC19, 2,  findBlanks)
JC19 <- JC19[, keepCols]

# missing values currently set as NaN
JC19 <- apply(JC19, 2, function(x) replace(x, is.nan(x), NA)) |>
  data.frame()

# TODO hake - piscivores
# check Jones and Checkley paper for data source?

# Stationarity tests ------------------------------------------------------

# example code copied from EcoRelease repo
  # adfTest(globlRich) 
  # acf(globlRich) 
  # rich_AR <- arima(globlRich, order = c(1, 0, 0)) # fit AR1 model
  # rich_e <- as.numeric(rich_AR$residuals)
  # acf(rich_e)
  # adfTest(rich_e) 
  # rich_AR$coef
