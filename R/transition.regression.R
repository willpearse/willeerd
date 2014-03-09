#' \code{transition.calc} regression of aged transitions among character states
#' 
#' @param tree a phylogeny (of class phylo) you wish to plot
#' @param gif name of GIF you would like to create
#' @param slices the time slices you would like to plot out; NULL (default) means one every unit time
#' @param colours colours for each speciation event (thus must be same length as 'slices'); NULL (default) uses 'rainbow'
#' @param pca PCA (of class prcomp) of phylogenetic dissimilarity matrix; NULL calculates one, I recommend you use the output from a previous run to speed things up
#' @param clade.mat clade matrix (of class clade.matrix) of phylogeny; NULL calculates one, I recommend you use the output from a previous run to speed things up
#' @param delay the delay between each slice's frame in the output GID; default 0.2 seconds
#' @details Probably best to just plot it out and see what happens. There are much smater ways of plotting out what species goes where, but this is what I've done... As with everything I have written, this is very much unchecked! Beware!!!
#' @return The data that were plotted last, the PCA and clade.matrix to speed later plots, and the colours used.
#' @author Will Pearse
#' @examples \dontrun{
#' tree <- rtree(50)
#' continuous <- rnorm(50); names(continuous) <- tree$tip.label
#' discrete <- factor(sample(letters[1:6], 50, replace=TRUE)); names(discrete) <- tree$tip.label
#' silly <- transition.calc(tree, continuous, discrete, simmap.nsim=10)
#' plot(silly)
#' }
#' @import phytools
#' @import caper
#' @export
transition.calc <- function(tree, continuous, discrete, simmap.model="ER", simmap.nsim=1000, simmap.pi="estimated", anc.ML.maxit=100000){
  #Assertions and argument checking
  if(!inherits(tree, "phylo")) stop("Error: '", deparse(substitute(simmap)), "' must be of class 'phylo'")
  if(!is.factor(discrete)) stop("Error: '", deparse(substitute(discrete)), "' must be a factor; preferably a discrete character!")
  if(is.null(names(discrete))) stop("Error: '", deparse(substitute(discrete)), "' must be named")
  if(!is.numeric(continuous)) stop("Error: '", deparse(substitute(continuous)), "' must be a numeric; preferably a continuous character!")
  if(is.null(names(continuous))) stop("Error: '", deparse(substitute(continuous)), "' must be named")
  if(!identical(sort(tree$tip.label), sort(names(discrete)))) stop("Error: mismatch between'", deparse(substitute(discrete)), "' and phylogeny")
  if(!identical(sort(tree$tip.label), sort(names(continuous)))) stop("Error: mismatch between'", deparse(substitute(continuous)), "' and phylogeny")
  
  #Make simmap
  simmap <- make.simmap(tree, discrete, model=simmap.model, nsim=simmap.nsim, pi=simmap.pi)
  t.ltt <- ltt(simmap[[1]], plot=FALSE, gamma=FALSE)
  
  #Find transitions (can be multiple per branch); pre-allocation could make too long a data.frame (CHECK!)
  states <- levels(discrete)
  tmp <- rep(NA, length(simmap)*length(simmap[[1]]$maps))
  transitions <- data.frame(from=tmp, to=tmp, node=tmp)
  x <- 1
  for(i in seq_along(simmap)){
      for(j in seq(from=length(simmap[[i]]$tip.label)+1, to=length(simmap[[i]]$maps))){
          if(length(simmap[[i]]$maps[[j]]) > 1){
              for(k in seq(from=1, length.out=length(simmap[[i]]$maps[[j]])-1)){
                  transitions$from[x] <- states[which(states==names(simmap[[i]]$maps[[j]])[k])]
                  transitions$to[x] <- states[which(states==names(simmap[[i]]$maps[[j]])[k+1])]
                  transitions$node[x] <- j
                  x <- x + 1
                  if(x >= nrow(transitions))
                      transitions <- rbind(transitions, data.frame(from=tmp, to=tmp, node=tmp))
              }
          } else {
              transitions$from[x] <- states[which(states==names(simmap[[i]]$maps[[j]])[1])]
              transitions$to[x] <- states[which(states==names(simmap[[i]]$maps[[j]])[1])]
              transitions$node[x] <- j
              x <- x + 1
              if(x == nrow(transitions))
                  transitions <- rbind(transitions, data.frame(from=tmp, to=tmp, node=tmp))
          }
      }
  }
  transitions <- transitions[!is.na(transitions$from),]
  #Age the transitions
  transitions$age <- max(t.ltt$times) - t.ltt$times[match(transitions$node, names(t.ltt$times))]
  transitions$first <- c(TRUE, rep(FALSE, nrow(transitions)-1))
  
  #Reconstruct continuous state
  anc.continuous <- anc.ML(tree, continuous, maxit=anc.ML.maxit)
  #Get the modal reconstructed nodal value and plot against that
  transitions$transition <- with(transitions, paste(from, to, sep="_"))
  #Prepare output and return
  output <- list(transitions=transitions, cont.sim=anc.continuous)
  class(output) <- "transition.calc"
  return(output)
}

#This has got a naughty which()[1] that couldn cause trouble...
plot.transition.calc <- function(x, ...){
    counts <- with(x$transitions, table(transition, node))
    modal.trans <- rownames(counts)[unlist(apply(counts, 2, function(x) which(max(x)==x)[1]))]
    plot(x$cont.sim$ace[-1] ~ factor(modal.trans))
}
