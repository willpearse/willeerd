% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/clean.taxonomy.R
\name{make.clean.taxon.lookup}
\alias{make.clean.taxon.lookup}
\title{\code{make.clean.taxon.lookup} Take input species names, quickly scrub them using The Plant List, and return a lookup table}
\usage{
make.clean.taxon.lookup(species)
}
\arguments{
\item{species}{Vector of species to be cleaned}
}
\value{
data.frame with the original species names ('original'), and the cleaned names if in TPL or NULL if not found in TPL ('clean'). Can be used as a lookup table in other functions, or in your own stuff.
}
\description{
\code{make.clean.taxon.lookup} Take input species names, quickly scrub them using The Plant List, and return a lookup table
}
\details{
Binds a phylogeny (donor) into a bigger phylogeny ('backbone'); useful if you're building a phylogeny a la Phylomatic.
}
\author{
Will Pearse
}

