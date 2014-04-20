prog.bar <- function(x, y)
    if((x %% floor(y/100)) == 0) if((x %% (floor(y/100)*10)) == 0) cat("|") else cat(".")
