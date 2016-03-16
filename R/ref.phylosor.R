#' \code{ref.phylosor} Calculate phylosor using a single assembalge as a reference
#' 
#' @param samp community matrix
#' @param reference community/assemblage to be used as a reference
#' @param tree ape::phylo phylogeny
#' @details Exactly as phylosor, except calculates the distance of everything in samp from the reference. A "good" way of dealing with large datasets, if you can't fit everything in memory at the same time.
#' @return A vector of the distances of the communities from the reference
#' @author Will Pearse, based on picante::phylosor
#' @examples \dontrun{
#' require(picante)
#' data(phylocom)
#' phylosor(phylocom$sample, phylocom$phylo)
#' ref.phylosor(phylocom$sample, phylocom$sample[6,], phylocom$phylo)
#' }
#' @importFrom picante pd
#' @importFrom ape is.rooted
#' @export
ref.phylosor <- function (samp, reference, tree){
    #Argument checking
    if(is.null(tree$edge.length))
        stop("Tree has no branch lengths, cannot compute PD")
    if(!is.rooted(tree))
        stop("Rooted phylogeny required for phylosor calculation")
    if(length(reference) != ncol(samp) | !is.null(dim(reference)))
        stop("Reference must be a vector with dimensions compatible with sample")
    samp <- as.matrix(samp)
    s <- nrow(samp)
    phylodist <- numeric(s)
    names(phylodist) <- rownames(samp)
    samp_comb <- matrix(NA, nrow=s, ncol=ncol(samp), dimnames=list(rownames(samp), colnames(samp)))
    for(i in seq(s))
        samp_comb[i,] <- samp[i,] + reference
    
    pdsamp <- pd(samp, tree)
    pdsamp_comb <- pd(samp_comb, tree)
    t <- matrix(c(reference, reference), nrow=2, byrow=TRUE)
    rownames(t) <- c("A", "B")
    colnames(t) <- names(reference)
    pdl <- pd(t, tree)$PD[1]
    for(i in seq(s)){
        pdk <- pdsamp[i, "PD"]
        pdcomb <- pdsamp_comb[i, "PD"]
        pdsharedlk <- pdl + pdk - pdcomb
        phylodist[i] = 2 * pdsharedlk/(pdl + pdk)
    }
    return(phylodist)
}
