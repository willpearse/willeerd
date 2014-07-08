\name{cartoon.plot}
\alias{cartoon.plot}
\title{\code{cartoon.plot} Plot phylogenies with cartoon-ised polytomies}
\usage{
cartoon.plot(tree, tip.groups = vector("list", 0), clade.col = NULL,
  br.clade.col = NULL, auto.polies = FALSE, ...)
}
\arguments{
  \item{tree}{ape::phylo object}

  \item{tip.groups}{list where each element are the tip
  *numbers* that should be in a cartoon polytomy (see
  auto.polies below)}

  \item{clade.col}{a vector of length tip.groups with
  colours for each polytomy OR a single colour for all
  polytomies OR NULL (the default) to make all polytomies
  rainbow-coloured}

  \item{br.clade.col}{the colours of the branches inside
  the cartoon. If NULL (the default, which you probably
  want), they don't show up in the plot.}

  \item{auto.polies}{generate tip.groups by making a
  cartoon polytomy for all terminal polytomies. Setting
  this to TRUE (default is FALSE) will over-ride any
  tip.groups information. Note the function's return values
  is you want to set up a weird interaction of polytomies}

  \item{...}{additional arguments to pass to
  ape::plot.phylo}
}
\value{
A list containing the tip.groups plotted, the clade colours
they were plotted under, the edges that were obscured by
the cartoon printing, and a list of the number of the node
that each polytomy started under
}
\description{
\code{cartoon.plot} Plot phylogenies with cartoon-ised
polytomies
}
\details{
Makes a simple 'cartoon' phylogeny where polytomies are
joined together in a big triangle. Could cause problems if
you're outputting to PDF, sorry!
}
\examples{
\dontrun{
require(ape)
tree <- read.tree(text="(((((A,B,C,D,E),(F,G,H,I,J)),H),K),L);")
cartoon.plot(tree, auto.polies=TRUE)
cartoon.plot(tree, list(1:5, 6:10), colours="grey30")
}
}
\author{
Will Pearse
}

