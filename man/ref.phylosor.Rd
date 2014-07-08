\name{ref.phylosor}
\alias{ref.phylosor}
\title{\code{ref.phylosor} Calculate phylosor using a single assembalge as a reference}
\usage{
ref.phylosor(samp, reference, tree)
}
\arguments{
  \item{samp}{community matrix}

  \item{reference}{community/assemblage to be used as a
  reference}

  \item{tree}{ape::phylo phylogeny}
}
\value{
A vector of the distances of the communities from the
reference
}
\description{
\code{ref.phylosor} Calculate phylosor using a single
assembalge as a reference
}
\details{
Exactly as phylosor, except calculates the distance of
everything in samp from the reference. A "good" way of
dealing with large datasets, if you can't fit everything in
memory at the same time.
}
\examples{
\dontrun{
require(picante)
data(phylocom)
phylosor(phylocom$sample, phylocom$phylo)
ref.phylosor(phylocom$sample, phylocom$sample[6,], phylocom$phylo)
}
}
\author{
Will Pearse, based on picante::phylosor
}

