library(fUnitRoots) # for timeseries stationarity tests

# TODO
# standardize data format
# test timeseries stationarity
# ordinate and inspect PCA

# Vet and format data -----------------------------------------------------

## Enviro vars -------------------------------------------------------------

## PP, SST, TOC, 13C, 15N

## Primary producers -------------------------------------------------------

# diatoms, silicoflagellates

## Forams ------------------------------------------------------------------

# planktic

## Fish --------------------------------------------------------------------

# anchovy, sardine, mesopelagic fish
# hake - piscivores

# Stationarity tests ------------------------------------------------------

# example code copied from EcoRelease repo
  # adfTest(globlRich) 
  # acf(globlRich) 
  # rich_AR <- arima(globlRich, order = c(1, 0, 0)) # fit AR1 model
  # rich_e <- as.numeric(rich_AR$residuals)
  # acf(rich_e)
  # adfTest(rich_e) 
  # rich_AR$coef
