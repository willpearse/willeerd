\name{factorise.tree}
\alias{factorise.tree}
\title{\code{factorise.tree} 'Factorise' a tree by removing out species, making it easier to plot/manipulate}
\usage{
factorise.tree(tree, scale.factor = 0.5)
}
\arguments{
  \item{tree}{ape::phylo phylogeny to be 'factorised'}

  \item{scale.factor}{multiplier for the number of species
  within a terminal polytomy. E.g., 0.1 means each terminal
  polytomy will be ~10% its current size}
}
\value{
List where first element is the factorised phylogeny, the
second the tips that were dropped from each node (on the
original tree; I can't guarantee the returning tree's
structure)u
}
\description{
\code{factorise.tree} 'Factorise' a tree by removing out
species, making it easier to plot/manipulate
}
\details{
Thins out additional species, making a phylogeny smaller by
reducing the size of each terminal polytomy by scale.factor
}
\examples{
\dontrun{
tree <- read.tree(text="((A,B,C,D,E),F);")
t <- factorise.tree(tree, 0.5)
plot(t$tree)
}
}
\author{
Will Pearse
}

