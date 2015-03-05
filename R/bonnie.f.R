#' \code{bonnie.f} Calculate nested F-ratios - "correctly" (for Bonnie Waring, anyway)
#' 
#' @param model standard linear model to be calculated
#' @param high Highest hierarchical-level in model
#' @param low Lowest hierarchical-level in model
#' @details Divides the F-value of the highest level by that of the interaction between it and the lowest-value. No p-value, you can do that yourself!
#' @return Numeric
#' @author Will Pearse
#' @export
bonnie.f <- function(model, high, low){
    high <- "top"
    low <- "bottom"
    y <- rnorm(100)
    top <- sample(letters[1:3], 100, TRUE)
    bottom <- sample(letters[4:6], 100, TRUE)
    model <- aov(y ~ top/bottom)
    coefs <- summary(model)[[1]]
    names <- gsub(" ", "", rownames(coefs))
    top.fac <- which(names == high)
    bottom.fac <- which(grepl(paste(high,low,sep=":"), names))
    return(coefs[top.fac,2] / coefs[bottom.fac,2])
}
