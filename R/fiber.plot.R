#' \code{fiber.plot} (fibrously) plots a phylogeny
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
#' data(perissodactyla)
#' output <- fiber.plot(perissodactyla.tree, "test.gif")
#' fiber.plot(perissodactyla.tree, "test.gif", clade.mat=output$clade.mat, pca=output$pca)
#' #...that was much quicker (...on bigger phylogenies...), because we used the PCA and clade.matrix from last time
#' }
#' @importFrom caper clade.matrix
#' @importFrom ape branching.times
#' @importFrom animation saveGIF
#' @export
fiber.plot <- function(tree, gif, slices=NULL, colours=NULL, pca=NULL, clade.mat=NULL, delay=0.2){
  #Assertions and argument checking
  if(!inherits(tree, "phylo"))
    stop("'", deparse(substitute(tree)), "' must be of class 'phylo'")
  if(!is.null(slices) & !is.numeric(slices))
      stop("'", deparse(substitute(slices)), "' must be a numeric!")
  if(!is.null(colours) & length(colours) != nrow(tree$edge))
    stop("'", deparse(substitute(colours)), "' is not long enough to colour the phylogeny")
  if(!is.null(pca) & !inherits(pca, "prcomp"))
    stop("'", deparse(substitute(pca)), "' must be of class 'prcomp'")
  if(!is.character(gif))
    stop("'", deparse(substitute(gif)), "' needs to be a filename!")
  #if(!is.null(clade.mat) & !inherits(clade.mat, "clade.matrix"))
    #stop("'", deparse(substitute(clade.mat)), "' must be of class 'clade.matrix'")
  
  #Setup
  timing <- branching.times(tree)
  if(is.null(slices))
    slices <- seq(from=max(timing), to=0)
  if(is.null(colours))
    colours <- rainbow(length(slices))
    #colours <- rainbow(nrow(tree$edge))
  #if(is.null(pca))
    #pca <- prcomp(cophenetic(tree), scale=TRUE, center=TRUE)
  if(is.null(clade.mat))
    clade.mat <- clade.matrix(tree)$clade.matrix
  #spp.val <- pca$x[,1]
  dimension <- floor(sqrt(length(spp.val))) + 1
  #data <- data.frame(species=names(spp.val), x=rep(seq(dimension), dimension)[seq_along(spp.val)], y=rep(seq(dimension), each=dimension)[seq_along(spp.val)])

  #n.groups <- 10; n.pc <- 1:5
  #clust <- hclust(dist(pca$x[,n.pc]))
  #groups <- cutree(clust, k=n.groups)
  #spp.val <- spp.val[order(groups, pca$x[,1])]
  #clade.mat <- clade.mat[,order(groups, pca$x[,1])]
  spp.val <- seq(5020)
  #Loop over all slices and print
  saveGIF({
    curr.colours <- rep(1, length(tree$tip.label))
    x <- 2
    curr.colours <- rep(NA, dimension^2)
    curr.colours[seq_along(spp.val)] <- 1
    image(matrix(curr.colours, nrow=dimension), col=colours, main=slices[1], zlim=c(1, length(colours)), bty="n", xaxt="n", yaxt="n", ylab="", xlab="")
    #contour(matrix(curr.colours, nrow=dimension), main="", bty="n", xaxt="n", yaxt="n", ylab="", xlab="", labels="", add=TRUE)
    for(i in seq_along(slices)[-1]){
      to.update <- names(timing)[timing > slices[i] & timing <= slices[i-1]]
      for(j in seq_along(to.update)){
        curr.colours[which(clade.mat[to.update[j],] == 1)] <- i
        x <- x + 1
      }
      image(matrix(curr.colours, nrow=dimension), col=colours, main=slices[i], zlim=c(1, length(colours)), bty="n", xaxt="n", yaxt="n", ylab="", xlab="")
      #contour(matrix(curr.colours, nrow=dimension), col=colours, main="", zlim=c(1, length(colours)), bty="n", xaxt="n", yaxt="n", ylab="", xlab="", labels="", add=TRUE)
    }
  }, interval = delay, movie.name = gif, ani.width = 600, ani.height = 600)
  
  #Invisibly return
  invisible(list(plot=matrix(curr.colours, nrow=dimension), colours=curr.colours, clade.mat=clade.mat, pca=pca))
}
