#Force an amount of significant digits.
#Outputs numbers as strings
#Useful for tables
spec_dec <- function(x, k) format(round(x, k), nsmall=k)