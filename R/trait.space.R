#' \code{trait.space} (maybe) make a trait space of species co-existence
#' 
#' @param comm Community matrix (communities in rows, species in column); NAs would probably be a bad idea
#' @param abundance Work with species abundances? Default FALSE
#' @param plot Produce a plot of the space at the end of the run? Default TRUE
#' @param coexist Provide a co-existence matrix to use for fitting (see below). Default to NULL to calculate one as part of the run. The matrix will be returned at the end of the run so you can work with it later. It is *definitely* worth calculating this once, as I've almost maliciously not optimised this code.
#' @param mask matrix of TRUE and FALSE, where only elements of the coexistence matrix where this is TRUE will be used in the calculations.
#' @param ... additional arguments to be passed to the optimisation routine. See the notes for advice on this.
#' @details Attempts to place all species in a 2-D space where distance between species is directly related to their proportional coexistence. I.e., if species are >=1 unit distant from each other, they never coexit, if they're on top of one-another, they always co-exist. See notes.
#' @note It is *strongly* advised that you check the convergence of this routine, and play around with the algorithm used to assess convergence. This came about as a result of a NutNet (http://www.nutnet.umn.edu/) meeting; please don't run away and use this without bringing me on-board, or at the very least informing the organisers of NutNet. Please don't become "that guy" and ignore this...
#' @return The data that were plotted last, the PCA and clade.matrix to speed later plots, and the colours used.
#' @author Will Pearse
#' @export
#' @importFrom fields rdist
trait.space <- function(comm, abundance=FALSE, plot=TRUE, coexist=NULL, mask=NULL, ...){
    #Internal functions
    prop.coexist <- function(x){
        #Life is too short to do matrix math
        output <- matrix(0, ncol=ncol(x), nrow=ncol(x))
        for(i in seq(ncol(x)))
            for(j in seq(ncol(x)))
                output[i,j] <- sum(x[,i] > 0 & x[,j] > 0) / sum(x[,i] > 0)
        return(1 - output)
    }
    prop.coexist.abundance <- function(x){
        #Life is too short to do matrix math
        output <- matrix(0, ncol=ncol(x), nrow=ncol(x))
        for(i in seq(ncol(x)))
            for(j in seq(ncol(x)))
                output[i,j] <- sum(x[x[,j]>0,i]) / sum(x[,i])
        return(1 - output)
    }
    lik.coexist <- function(x, mat, mask=NULL){
        logit <- function(x) 1 / (1 + exp(-x))
        points <- cbind(x[seq(ncol(mat))], x[seq(from=ncol(mat)+1, length.out=ncol(mat))])
        dist <- rdist(points)
        if(is.null(mask))
            dist <- logit(abs(dist - mat)) else dist <- logit(abs(dist[mask] - mat[mask]))
        return(sum(dist, na.rm=TRUE))
    }
    
    #Make a coexistence matrix if necessary
    if(is.null(coexist)){
        if(abundance) coexist <- prop.coexist(comm) else
            coexist <- prop.coexist.abundance(comm)
    }

    #Fit parameters and neaten results
    results <- optim(par=rnorm(ncol(comm)*2), lik.coexist, mat=coexist, mask=mask, ...)
    neat.results <- cbind(results$par[seq(ncol(comm))], results$par[seq(from=ncol(comm)+1, length.out=ncol(comm))])
    rownames(neat.results) <- colnames(comm)
    
    #Plot if necessary
    if(plot == TRUE){
        plot(neat.results[,1] ~ neat.results[,2], type="n")
        text(neat.results[,1] ~ neat.results[,2], labels=rownames(neat.results))
    }

    #Return
    return(list(results=neat.results, coexist=coexist))
}
