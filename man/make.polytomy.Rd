% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/clean.taxonomy.R
\name{make.polytomy}
\alias{make.polytomy}
\title{\code{make.polytomy} Make a polytomy (optionally with edge lengths)}
\usage{
make.polytomy(species, tip.length = NA)
}
\arguments{
\item{species}{tip.labels for polytomy}

\item{tip.length}{edge length for polytomy}
}
\value{
ape::phylo polytomy
}
\description{
\code{make.polytomy} Make a polytomy (optionally with edge lengths)
}
\details{
A light wrapper around some ape functions to make a polytomy. If you need it you'll know.
}
\author{
Will Pearse
}

