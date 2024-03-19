# This script houses helper functions called in multiple places
# throughout the data processing and analysis workflow

# Clean a character vector (e.g. column names) containing periods -
# introduced by conversion when reading Excel special characters into R
cleanNames <- function(char){
  # substitute special characters and spaces converted to periods from Excel
  char <- gsub(x = char, pattern = '\\.\\.', replacement = '_')  
  char <- gsub(x = char, pattern = '\\.', replacement = '_')  
  # omit trailing underscore(s) in character strings where present
  trails <- str_ends(char, '_')
  char[trails] <- char[trails] |>
    str_sub(end = -2)
  # deal with case of double trailing underscores
  trails2 <- str_ends(char, '_')
  char[trails2] <- char[trails2] |>
    str_sub(end = -2) 
  return(char)
}

# Match rows in data to rows of target composite df, based on age value:
# - exclude any missing matches and round non-integer years
# - return column(s) of target variables with NAs for non-measured study years
# tmplt = vector/column of years from target composite dataframe
# yrCol = name or position of column in dat AND tmplt with age/year data
# xtrctCol = name(s) of column(s) in dat to add to composite df
matchTime <- function(dat, tmplt, yrCol = 'year', xtrctCol){
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

# Generate filepath name for charts to save externally
nameOutput <- function(string, ext = 'pdf'){
  day <- as.Date(date(), format="%a %b %d %H:%M:%S %Y")
  paste0('figs/', string, '_', day, '.', ext)
}
