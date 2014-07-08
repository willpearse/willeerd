\name{transition.calc}
\alias{transition.calc}
\title{\code{transition.calc} regression of aged transitions among character states}
\usage{
transition.calc(tree, continuous, discrete, simmap.model = "ER",
  simmap.nsim = 1000, simmap.pi = "estimated", anc.ML.maxit = 1e+05)
}
\arguments{
  \item{tree}{a phylogeny (of class phylo) you wish to
  plot}

  \item{gif}{name of GIF you would like to create}

  \item{slices}{the time slices you would like to plot out;
  NULL (default) means one every unit time}

  \item{colours}{colours for each speciation event (thus
  must be same length as 'slices'); NULL (default) uses
  'rainbow'}

  \item{pca}{PCA (of class prcomp) of phylogenetic
  dissimilarity matrix; NULL calculates one, I recommend
  you use the output from a previous run to speed things
  up}

  \item{clade.mat}{clade matrix (of class clade.matrix) of
  phylogeny; NULL calculates one, I recommend you use the
  output from a previous run to speed things up}

  \item{delay}{the delay between each slice's frame in the
  output GID; default 0.2 seconds}
}
\value{
The data that were plotted last, the PCA and clade.matrix
to speed later plots, and the colours used.
}
\description{
\code{transition.calc} regression of aged transitions among
character states
}
\details{
Probably best to just plot it out and see what happens.
There are much smater ways of plotting out what species
goes where, but this is what I've done... As with
everything I have written, this is very much unchecked!
Beware!!!
}
\examples{
\dontrun{
tree <- rtree(50)
continuous <- rnorm(50); names(continuous) <- tree$tip.label
discrete <- factor(sample(letters[1:6], 50, replace=TRUE)); names(discrete) <- tree$tip.label
silly <- transition.calc(tree, continuous, discrete, simmap.nsim=10)
plot(silly)
}
}
\author{
Will Pearse
}

