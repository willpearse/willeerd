\name{congeneric.merge}
\alias{congeneric.merge}
\title{\code{congeneric.merge} Bind missing species into a phylogeny based on taxonomy}
\usage{
congeneric.merge(lookup, tree, split = "_", cite = TRUE)
}
\arguments{
  \item{lookup}{Either a data.frame from
  'make.clean.taxon.lookup' or a vector of species names to
  be bound into the tree if missing from it}

  \item{tree}{The phylogeny to have those species inserted
  into it}

  \item{cite}{Set to FALSE to have it stop bleating at you
  to cite phyloGenerator.}
}
\value{
The bound phylogeny
}
\description{
\code{congeneric.merge} Bind missing species into a
phylogeny based on taxonomy
}
\details{
Binds missing species into a phylogeny by replacing that
clade , which this function (and all its internal calls)
}
\author{
Will Pearse
}

