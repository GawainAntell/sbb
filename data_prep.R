library(stringr)
library(tabulizer) # read PDF tables; requires JDK vsn of Java, and rJava pkg
library(ggplot2)
source('utils.R') # for custom helper functions

# Vet and format data -----------------------------------------------------

## Enviro vars -------------------------------------------------------------

# PP, SST - JP Kennett's data included with Jones and Checkley otoliths (below)

# TOC - Berger 2004

nmBerg04 <- 'data-raw/Berger-et-al-2004_Baumgartner-1992.pdf'
Berg04raw <- extract_tables(nmBerg04, pages = 21:40) |>
  lapply( function(l) apply(l, c(1,2), as.numeric) )
# all entries treated as character; must convert type before rbind on rows

Berg04 <- do.call(rbind, Berg04raw) |>
  data.frame()
colnames(Berg04) <- c('year', 'TOC', 'TOCdetrend')
# write.csv(Berg04, 'data/Berger-et-al-2004_cleaned.csv', row.names = FALSE)

# ENSO - Li et al. 2011

Li11 <- read.table('data-raw/Li-et-al-2011.txt', 
                   header = TRUE, comment.char = '#')
colnames(Li11)[1] <- 'year'
# write.csv(Li11, 'data/Li-et-al-2011_cleaned.csv', row.names = FALSE)

# biogenic opal - Barron 2013

Barr13o_splt <- read.csv('data-raw/Barron-et-al-2013-opal.csv')
# entire dataset is split at 1884, and second half pasted as columns to right :|
dupeColsO <- str_ends(colnames(Barr13o_splt), '1')
topHlfO <- Barr13o_splt[, !dupeColsO]
topHlfO<- topHlfO[ , -ncol(topHlfO)] # blank divider column
lwrHlfO <- Barr13o_splt[, dupeColsO]
colnames(lwrHlfO) <- str_sub(colnames(lwrHlfO), end = -3) # strip '.1' suffix
Barr13o <- rbind(topHlfO, lwrHlfO)
# spaces converted to periods from excel, and periods appear in gen abbreviations
colnames(Barr13o) <- gsub(x = colnames(Barr13o), pattern = '\\.\\.', replacement = '_')  
colnames(Barr13o) <- gsub(x = colnames(Barr13o), pattern = '\\.', replacement = '_') 
# omit trailing space in names where present
trails <- str_ends(colnames(Barr13o), '_')
colnames(Barr13o)[trails] <- colnames(Barr13o)[trails] |>
  str_sub(end = -2)
omitCols <- c('Varve_range', 'Bottom_varve')
Barr13o <- Barr13o[ , !colnames(Barr13o) %in% omitCols]
# remove rows at bottom of original spreadsheet - no values
findBlanks <- function(v){
  any( !is.na(v) )
} 
keepRowsO <- apply(Barr13o, 1,  findBlanks) # do this after removing metadata cols
Barr13o <- Barr13o[keepRowsO, ]
colnames(Barr13o)[1] <- 'year' # this is the oldest (bottom) year of the sample
# write.csv(Barr13o, 'data/Barron-et-al-2013-opal_cleaned.csv', row.names = FALSE)

# TOC - Wang et al. 2017

# floods/droughts - Sarno 2020

# 13C, 15N

## Primary producers -------------------------------------------------------

# diatoms - Barron et al. 2013

Barr13d <- read.csv('data-raw/Barron-et-al-2013-diatoms.csv')
# spaces converted to periods from excel, and periods appear in gen abbreviations
colnames(Barr13d) <- gsub(x = colnames(Barr13d), pattern = '\\.\\.', replacement = '_')  
colnames(Barr13d) <- gsub(x = colnames(Barr13d), pattern = '\\.', replacement = '_')  
# omit trailing space in names where present
trailsD <- str_ends(colnames(Barr13d), '_')
colnames(Barr13d)[trailsD] <- colnames(Barr13d)[trailsD] |>
  str_sub(end = -2)
# some columns contain metadata, not species abundances (or not relevant taxa)
omitCols <- c(omitCols, # varve range, bottom varve number
              'Benthic', 'Freshwater_planktic', 'Reworked', 'Total_counted')
Barr13d <- Barr13d[ , !colnames(Barr13d) %in% omitCols]
# final rows are completely empty
keepRowsD <- apply(Barr13d, 1,  findBlanks)
Barr13d <- Barr13d[keepRowsD, ]
colnames(Barr13d)[1] <- 'year' # this is the oldest (bottom) year of the sample
# write.csv(Barr13d, 'data/Barron-et-al-2013-diatoms_cleaned.csv', row.names = FALSE)

# silicoflagellates - Barron et al. 2013

Barr13s_splt <- read.csv('data-raw/Barron-et-al-2013-silicoflagellates.csv')
# entire dataset is split at 1884, and second half pasted as columns to right :|
dupeColsS <- str_ends(colnames(Barr13s_splt), '1')
topHlfS <- Barr13s_splt[, !dupeColsS]
topHlfS <- topHlfS[ , -ncol(topHlfS)] # blank divider column
lwrHlfS <- Barr13s_splt[, dupeColsS]
colnames(lwrHlfS) <- str_sub(colnames(lwrHlfS), end = -3) # strip '.1' suffix
Barr13s <- rbind(topHlfS, lwrHlfS)
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
# write.csv(Barr13s, 'data/Barron-et-al-2013-silicoflagellates_cleaned.csv', row.names = FALSE)

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
# write.csv(Baum92, 'data/Baumgartner-1992_cleaned.csv', row.names = FALSE)

# mesopelagic fish - Jones and Checkley 2019
# otolith deposition rate for 5 families, and SST and PP proxies from JP Kennett

jc19 <- read.csv('data-raw/Jones-Checkley-2019-decadal.csv') 
colnames(jc19)[1] <- 'year' # rename from 'year_ad' to match other datasets

# last 4 columns are empty; identify columns with data
keepCols <- apply(jc19, 2,  findBlanks)
jc19 <- jc19[, keepCols]

# missing values currently set as NaN
jc19 <- apply(jc19, 2, function(x) replace(x, is.nan(x), NA)) |>
  data.frame()
# write.csv(jc19, 'data/Jones-Checkley-2019_cleaned.csv', row.names = FALSE)

# Data distribution plot --------------------------------------------------

timescale <- 0:2010 # 1748:2010
plotDat <- data.frame('year' = timescale)
datNms <- c('enso', 'toc', 'SST, PP', 'opal', 
            'diatoms', 'silicos', 'scales', 'otoliths')
for (dat in list(Li11, Berg04, jc19, Barr13o, Barr13d, Barr13s, Baum92, jc19)){
  hasDat <- plotDat[ ,'year'] %in% dat[ ,'year'] |>
    as.numeric()
  plotDat <- cbind(plotDat, hasDat)
}
colnames(plotDat)[-1] <- datNms

plotDatLng <- reshape(plotDat, direction = 'long',
                       v.names = 'sampled', varying = datNms,
                       timevar = 'dataset', times = datNms,
                       idvar = 'year')
plotDatLng$dataset <- factor(plotDatLng$dataset, 
                             levels = datNms,
                             ordered = TRUE)


plotPrint <- ggplot(plotDatLng[ plotDatLng$sampled==1, ], 
       aes(x = year, y = dataset)) +
  theme_bw() +
  scale_x_continuous(name = 'year AD', 
                     limits = range(timescale),
                     expand = expansion(add = c(0, 20)),
                     breaks = seq(0, 2000, by = 250)) +
  geom_point(size = 0.3)

plotNm <- nameOutput('tseries-rangethrough-durations', 'pdf')
pdf(plotNm, width = 6, height = 4)
print(plotPrint)
dev.off()
