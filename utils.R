# This script houses helper functions called in multiple places
# throughout the data processing and analysis workflow

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
