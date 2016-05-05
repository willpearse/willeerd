#' Dropping missing data from data.frames
#'
#' It's often hard to figure out what data to drop, particularly in
#' comparative datasets where we have lots of trait-data for species,
#' but some species are missing a lot of data and some traits are
#' missing for a lot of species. These functions help you choose which
#' species and traits to drop, while trying to minimise the total
#' information lost. \code{auto.drop} uses a greedy algorithm to pick
#' species and traits to drop, and returns a \code{data.frame} of the
#' species and traits to be dropped (along with the order in which
#' they were picked). \code{run.drop} actually performs the drop. As
#' the examples show, these are not difficult functions to use, and
#' the code for them is quite transparent.
#' @param data \code{data.frame} to be checked or dropped. Output will
#' make more sense (and likely be safer) if it's got rownames.
#' @param stop.density stop dropping when the information density is
#' greater than or equal to this number. A density of 1 means
#' everything in the modified \code{data} input has no \code{NA}s.
#' @param drop.history output from \code{auto.drop} to be run through
#' @param na.omit whether to calculate density on na.omit-ted data,
#' and whether to run \code{\link{na.omit}} on datasets before
#' returning. If you want to do some sort of analysis on your data,
#' you almost certainly want this option.
#' @return \code{auto.drop}: a data.frame with the (greedy) order in
#' which species and data.frame columns were chosen to be
#' dropped. \code{run.drop}: your cleaned \code{data.frame}.
#' @examples
#' data <- data.frame(x=rnorm(1000), y=rnorm(1000,3), z=rnorm(1000,-3))
#' data$x[data$x < -2] <- NA
#' data$z[data$z < -2] <- NA
#' new.data <- run.drop(data, auto.drop(data))
#' #...the 'z' column was giving us trouble, and so has been dropped,
#' #   retaining most of the rows in our dataset.
#' @rdname data.drop
#' @name data.drop
#' @export
auto.drop <- function(data, stop.density=0.75, na.omit=TRUE){
    #Internal functions
    rank.density <- function(data, dim=1){
        if(!xor(dim == 1, dim == 2))
            stop("'dim' must be 1 (species) or 2 (traits)")
        if(dim == 1) other.dim <- 2 else other.dim <- 1
        improvement <- setNames(
            apply(data, dim, function(x) sum(is.na(x))) / dim(data)[other.dim],
            dimnames(data)[[dim]]
        )
        return(sort(improvement, decreasing=TRUE))
    }

    data.density <- function(x, na.omit=TRUE){
        if(na.omit){
            max.info <- prod(dim(x))
            x <- na.omit(x)
            return(sum(!is.na(as.vector(x))) / max.info)
        }
        return(sum(!is.na(as.vector(x))) / prod(dim(x)))
    }
    
    #Setup
    history <- rep(NA, prod(dim(data))+1)
    history <- data.frame(species=history, trait=history, density=history, species.gain=history, trait.gain=history)
    history$density[1] <- data.density(data, na.omit=TRUE)

    #Do work
    for(i in seq(2:nrow(history))){
        spp <- rank.density(data, 1)
        traits <- rank.density(data, 2)
        history$trait.gain[i] <- traits[1]
        history$species.gain[i] <- spp[1]
        if(spp[1] > traits[1]){
            history$species[i] <- names(spp[1])
            data <- data[-which(rownames(data)==names(spp[1])),]
        } else {
            history$trait[i] <- names(traits[1])
            data <- data[,-which(names(data)==names(traits[1]))]
        }
        history$density[i] <- data.density(data, na.omit=TRUE)
        if(history$density[i] >= stop.density)
            break
    }

    #Cleanup and return
    history <- history[!is.na(history$density),]
    return(history)
}
#' @rdname pez.metrics
#' @name pez.metrics
#' @export
run.drop <- function(data, drop.history, na.omit=TRUE){
    traits <- Filter(Negate(is.na), drop.history$trait)
    spp <- Filter(Negate(is.na), drop.history$species)
    data <- data[!rownames(data) %in% spp, !names(data) %in% traits]
    if(na.omit)
        data <- na.omit(data)
    return(data)
}
