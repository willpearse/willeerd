prog.bar <- function(x, y)
    tryCatch(if((x %% floor(y/100)) == 0) if((x %% (floor(y/100)*10)) == 0) cat("|") else cat("."), error=function(z) cat("."))
