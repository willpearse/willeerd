#' \code{bind.ultrametric.by.replacement} Bind one phylogeny into another, resulting in an ultrametric phylogeny
#' 
#' @param backbone the backbone phylogeny into which the donor is to be bound
#' @param donor the phylogeny to bound into the backbone phylogeny
#' @param replacing.tip.label the species in the donor phylogeny that's being replaced by the donor phylogeny
#' @param donor.length how deep the donor phylogeny should be cut into the backbone phylogeny. If NA (default), then the bladj algorithm is followed (or, in plain English, it's put half-way along the branch)
#' @details Binds a phylogeny (donor) into a bigger phylogeny ('backbone'); useful if you're building a phylogeny a la Phylomatic.
#' @return The bound phylogeny
#' @author Will Pearse
#' @importFrom ape bind.tree
bind.ultrametric.by.replacement <- function(backbone, donor, replacing.tip.label, donor.length=NA){	
    bind.point <- which(backbone$tip.label == replacing.tip.label)
    backbone <- bind.tree(backbone, donor, where=bind.point)
    which.tip <- which(backbone$tip.label == donor$tip.label[1])
    which.node <- backbone$edge[which(backbone$edge[,2] == which.tip),1]
    which.edge <- which(backbone$edge[,2] == which.node)
    tip.length <- backbone$edge.length[which.edge]
    if(is.na(donor.length)){
        backbone$edge.length[which.edge] <- tip.length/2
    } else {
        backbone$edge.length[which.edge] <- tip.length - donor.length/2
    }
    return(backbone)	
}	

#' \code{make.polytomy} Make a polytomy (optionally with edge lengths)
#' 
#' @param species tip.labels for polytomy
#' @param tip.length edge length for polytomy
#' @details A light wrapper around some ape functions to make a polytomy. If you need it you'll know.
#' @return ape::phylo polytomy
#' @author Will Pearse
#' @importFrom ape as.phylo.formula
make.polytomy <- function(species, tip.length=NA){	
    d.f <- data.frame(spp=factor(species))	
    polytomy <- as.phylo.formula(~spp, data=d.f)	
    if(!is.na(tip.length)) polytomy$edge.length <- rep(tip.length, length(species))	
    return(polytomy)	
}	

#' \code{find.unique.branch.length} Find the terminal edge length leading to a species in a phylogeny
#' 
#' @param tree the phylogeny to be searched
#' @param tip the tip whose edge length must be found
#' @details If you find yourself manually manipulating phylogenies a lot, this is for you
#' @return The edge length leading to that tip
#' @author Will Pearse
find.unique.branch.length <- function(tree, tip){	
    which.tip <- which(tree$tip.label == tip)
    which.edge <- which(tree$edge[,2] == which.tip)
    tip.length <- tree$edge.length[which.edge]
    return(tip.length)	
}

#' \code{make.clean.taxon.lookup} Take input species names, quickly scrub them using The Plant List, and return a lookup table
#' 
#' @param species Vector of species to be cleaned
#' @details Binds a phylogeny (donor) into a bigger phylogeny ('backbone'); useful if you're building a phylogeny a la Phylomatic.
#' @return data.frame with the original species names ('original'), and the cleaned names if in TPL or NULL if not found in TPL ('clean'). Can be used as a lookup table in other functions, or in your own stuff.
#' @author Will Pearse
#' @importFrom Taxonstand TPLck
#' @importFrom plyr ldply
#' @export
make.clean.taxon.lookup <- function(species){
    #Basic name scrubbing
    original <- species
    species <- gsub("_", " ", species)
    species <- gsub("  ", " ", species)
    
    #Run through TPL
    search <- ldply(species, function(x) TPLck(x))

    #Make lookup (including missing) and return
    lookup <- with(search, data.frame(original=original, clean=ifelse(Plant.Name.Index, paste(New.Genus, New.Species), NULL)))
    return(lookup)
}

#' \code{congeneric.merge} Bind missing species into a phylogeny based on taxonomy
#' 
#' @param lookup Either a data.frame from 'make.clean.taxon.lookup' or a vector of species names to be bound into the tree if missing from it
#' @param tree The phylogeny to have those species inserted into it
#' @param cite Set to FALSE to have it stop bleating at you to cite phyloGenerator.
#' @param split split by which tip.labels will be cut to find congeners
#' @details Binds missing species into a phylogeny by replacing that
#' clade , which this function (and all its internal calls)
#' @return The bound phylogeny
#' @author Will Pearse
#' @importFrom ape drop.tip
congeneric.merge <- function(lookup, tree, split="_", cite=TRUE){
    if(cite)
        cat("\nCite phyloGenerator when using this: DOI-10.1111/2041-210X.12055")
    if(!is.data.frame(lookup))
        lookup <- data.frame(clean=lookup, stringsAsFactors=FALSE)
    before <- sum(lookup$clean %in% tree$tip.label)
    for(i in seq(nrow(lookup))){
        prog.bar(i, nrow(lookup))
        if(!is.null(lookup$clean[i]) & !lookup$clean[i] %in% tree$tip.label){
            genus <- strsplit(lookup$clean[i], split, fixed=TRUE)[[1]][1]
            matches <- unique(grep(genus, tree$tip.label, value=TRUE))
            if(length(matches) > 0){
                tree <- drop.tip(tree, matches[-1])
                tip.length <- find.unique.branch.length(tree, matches[1])
                polytomy <- make.polytomy(unique(c(matches, lookup$clean[i])), (tip.length/2))
                tree <- bind.ultrametric.by.replacement(tree, polytomy, matches[1], tip.length)
            }
        }
    }
    cat("\nNumber of species in tree before:", before)
    cat("\nNumber of species in tree now:   ", sum(lookup$clean %in% tree$tip.label), "\n")
    return(tree)
}

#' \code{make.composite.with.polytomies} Bind groups of species into a phylogeny as polytomies
#' 
#' @param tree the backbone phylogeny into which the species are to be placed
#' @param genera vector, of same length as species, that delimits the groups that the species are to be placed together with
#' @param species the species that will be bound into the phylogeny
#' @param max.genus.age the maximum age of a group (genus) to be bound into the phylogeny; if NA (the default) there is no maximum age, and it just follows the bladj algorithm
#' @details Binds species into the phylogeny using the values in 'genera' to group the species into polytomies. Each entry in 'genus' should be present in the tree, since the polytomy will be bound in replacing this genus.
#' This allows you, essentially, to reproduce much of the functionality of Phylomatic. It's intended to be used to replace a representative of a genus in a phylogeny with lots of members of that genus. See the example for usage.
#' @return The bound phylogeny
#' @author Will Pearse
#' @examples \dontrun{
#' require(ape)
#' tree <- read.tree(text='((Areplacement:10,GenusB:10):10, GenusC:10);')
#' genera <- c("Areplacement","Areplacement","Areplacement","GenusB","GenusB","GenusC")
#' species <- c("A robur","A ilex","A crud","B homo","B sapiens","C us")
#' tree <- make.composite.with.polytomies(tree, genera, species)
#' }
make.composite.with.polytomies <- function(tree, genera, species, max.genus.age=NA){
    genera <- as.character(genera)
    species <- as.character(species)
    for(genus in unique(genera)){
        species.to.bind <- species[genera == genus]
        if(length(species.to.bind) == 1){
            tree$tip.label[tree$tip.label == genus] <- species.to.bind
        } else {
            tip.length <- find.unique.branch.length(tree, genus)
            edge.warning <- NA
            if(!is.na(max.genus.age)){
                if(max.genus.age*2 < tip.length){
                    tip.length <- min(tip.length, max.genus.age*2)
                    edge.warning <- tip.length
                }
            }
            polytomy <- make.polytomy(species.to.bind, (tip.length/2))
            tree <- bind.ultrametric.by.replacement(tree, polytomy, genus, edge.warning)
        }
    }
    return(tree)
}

#' \code{randomly.resolve} Moves species around randomly within a genus
#' 
#' @param tree phylogeny containing the species to be shunted around
#' @param to.resolve species to be shunted around
#' @param split split by which tip.labels will be cut to find congeners
#' @details Moves species around within a genus, binding them back in
#' as a sister to a species in the same genus. Does nothing if a
#' specified species is the only member of a genus in the tree.  This
#' isn't the perfect way to deal with species for which you have no
#' DNA data, but it's definitely better than just sticking them in at
#' the base of a genus (which makes no sense to me). Do this many,
#' many times! If you're desperate to move things around at a depth
#' different from the genus, either modify the code (...) or use
#' Family_Genus_species names in your phylogeny.
#' @return The randomised phylogeny
#' @author Will Pearse
#' @importFrom ape drop.tip
#' @export
randomly.resolve <- function(tree, to.resolve, split="_"){
    for(species in to.resolve){
        genus <- strsplit(species, split)[[1]][1]
        if(sum(grepl(genus, tree$tip.label)) > 1){
            tree <- drop.tip(tree, species)
            sister <- sample(grep(genus, tree$tip.label, value=TRUE), 1)
            edge.length <- tree$edge.length[which(tree$edge[,2] == which(tree$tip.label == sister))]
            poly <- make.polytomy(c(sister, species), edge.length/2)
            tree <- bind.ultrametric.by.replacement(tree, poly, sister)
        }
    }
    return(tree)
}
