% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/phylo.plots.R
\name{tipring}
\alias{tipring}
\title{\code{tipring} Add tiplabels to a circular phylogeny}
\usage{
tipring(tips, col, radial.adj = 1, ...)
}
\arguments{
\item{tips}{Numbers of tips in phylogeny's \code{tip.labels} to be plotted}

\item{col}{Colour (can be of length > 1) to plot tips}

\item{radial.adj}{adjustment factor (multiplier) for spacing of tip-labels}

\item{...}{additional arguments passed to \code{\link[graphics]{lines}}}
}
\value{
...exactly as nodelabels
}
\description{
\code{tipring} Add tiplabels to a circular phylogeny
}
\details{
A way of getting evenly spaced node labels on a radial phylogeny. Just run the examples, it's quite simple. This is a very minorly-editted version of nodelabels; I can't take much credit! Don't cite this, cite ape!
}
\examples{
\dontrun{
tree <- stree(128, type="balanced")
plot(tree, type="radial", show.tip.label=FALSE)
tipring(pch=20)
tipring(pch=20, radial.adj=1.05, col="red")
tipring(pch=20, radial.adj=1.1, col="blue")
}
}
\author{
Will Pearse
}

