# return number of decimals 
# https://stackoverflow.com/questions/5173692/how-to-return-number-of-decimal-places-in-r/5173906
.decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}


substrRight <- function(x, n){
  substr(x, nchar(x) - n + 1, nchar(x))
}