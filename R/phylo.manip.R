#' \code{label.tree.tips} Label tips of phylogeny according to node labels
#' @param tree ape::phylo object; must have node.label slot (otherwise there're no node labels to use!)
#' @param clade.mat If not NULL (the default), the \code{clade.mat} slot of a caper::clade.matrix function call on the phylogeny to be labelled. You may already have this lying around, if so it speeds the function somewhat.
#' @details Returns a named vector describing what node.label the tips in your phylogeny subtends from. Each tip will only be labelled by one node - the youngest. This is essentially a way of figuring out if a species is in a particular clade (e.g., Fagaceae) from \code{node.label}s in a phylogeny.
#' @return Named vector describing tree.
#' @author Will Pearse
#' @importFrom caper clade.matrix
#' @export
label.tree.tips <- function(tree, clade.mat=NULL){
    #Assertions and argument checking
    if(is.null(tree$node.label)) stop("Phylogeny must have internally-labelled nodes")
    if(!inherits(tree, "phylo")) stop("Phylogeny must be an ape::phylo object")
    if(is.null(clade.mat)) clade.mat <- clade.matrix(tree)$clade.matrix

    #Label tips according to internal node labels
    named.nodes <- which(tree$node.label != "")
    names(named.nodes) <- tree$node.label[named.nodes]
    ages.nodes <- sapply(named.nodes, function(x) sum(tree$edge.length[clade.mat[x+length(tree$tip.label),] == 1]))
    named.nodes <- named.nodes[order(ages.nodes, decreasing=TRUE)]
    groups <- rep("background", length(tree$tip.label))
    for(i in seq_along(named.nodes))
        groups[clade.mat[named.nodes[i]+length(tree$tip.label),] == 1] <- names(named.nodes[i])
    names(groups) <- tree$tip.label   

    return(groups)
}

#' \code{label.tree.nodes} Label nodes of phylogeny according to tip labels
#' @param tree ape::phylo object; must have node.label slot (otherwise there're no node labels to use!)
#' @param groups The groups to be painted onto the nodes
#' @param clade.mat If not NULL (the default), the \code{clade.mat} slot of a caper::clade.matrix function call on the phylogeny to be labelled. You may already have this lying around, if so it speeds the function somewhat.
#' @details Returns a named vector describing what the node.label of a phylogeny should be to match the tips. It's essentially a ppor man's ancestral state reconstruction, such that you can paint values for species onto the clade that contains them all on a phylogeny. For instance, if you know thirty species are in Fagaceae, it will paint those labels onto a clade subtending from them.
#' @note A friend of \code{label.tree.tips}
#' @return Vector of labels; could be put straight into an ape::phylo object in the node.label slot.
#' @author Will Pearse
#' @importFrom caper clade.matrix
#' @export
label.tree.nodes <- function(tree, groups, clade.mat=NULL){
    #Assertions and argument checking
    if(is.null(tree$node.label)) stop("Phylogeny must have internally-labelled nodes")
    if(!inherits(tree, "phylo")) stop("Phylogeny must be an ape::phylo object")
    if(is.null(clade.mat)) clade.mat <- clade.matrix(tree)$clade.matrix

    #Go along and label
    node.groups <- numeric(nrow(clade.mat) - ncol(clade.mat))
    for(i in seq_along(node.groups)){
        names <- unique(groups[clade.mat[i+ncol(clade.mat),]==1])
        if(length(names) > 1) node.groups[i] <- NA else node.groups[i] <- names
    }
    
    return(node.groups)
}

#' \code{drip.node.labels} 'Drip' node.labels through a phylogeny
#' @param tree ape::phylo object; must have node.label slot (otherwise there're no node labels to use!)
#' @param clade.mat If not NULL (the default), the \code{clade.mat} slot of a caper::clade.matrix function call on the phylogeny to be labelled. You may already have this lying around, if so it speeds the function somewhat.
#' @details Returns a vector that names all nodes according to the label of their ancestral nodes. So if the Fagaceae node is labelled, all nodes within that clade will be labelled Fagaceae also. A new node label, such a Quercus, will be used instead if one is encountered.
#' @note Convenient if you do a lot of OU modelling...
#' @return Vector of labels; could be put straight into an ape::phylo object in the node.label slot.
#' @author Will Pearse
#' @importFrom caper clade.matrix
#' @importFrom phytools nodeheight
#' @export
drip.node.labels <- function(tree, clade.mat=NULL){
    #Assertions and argument checking
    if(is.null(tree$node.label)) stop("Phylogeny must have internally-labelled nodes")
    if(!inherits(tree, "phylo")) stop("Phylogeny must be an ape::phylo object")
    if(is.null(clade.mat)) clade.mat <- clade.matrix(tree)$clade.matrix

    #Order node labels according to age
    named.nodes <- which(tree$node.label != "")
    names(named.nodes) <- tree$node.label[named.nodes]
    ages.nodes <- sapply(named.nodes, function(x) nodeheight(tree, which(tree$node.label==x)))
    named.nodes <- named.nodes[order(ages.nodes)]

    #Loop through nodes and label down the tree (slow, I know...)
    node.labels <- character(tree$Nnode)
    for(i in seq_along(named.nodes)){
        next.edges <- numeric(0)
        curr.edges <- tree$edge[tree$edge[,1]==named.nodes[i]+length(tree$tip.label),2]
        curr.edges <- curr.edges[curr.edges > length(tree$tip.label)]
        while(length(curr.edges) > 1){
            for(j in seq_along(curr.edges)){
                node.labels[curr.edges[j]] <- names(named.nodes[i])
                t <- tree$edge[tree$edge[,1]==curr.edges[j],2]
                next.edges <- append(next.edges, t[t > length(tree$tip.label)])
            }
            curr.edges <- next.edges
            next.edges <- numeric(0)
        }
    }

    return(node.labels[-seq_along(tree$tip.label)])
}
