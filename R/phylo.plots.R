#' \code{cartoon.plot} Plot phylogenies with cartoon-ised polytomies
#' 
#' @param tree ape::phylo object
#' @param tip.groups list where each element are the tip *numbers* that should be in a cartoon polytomy (see auto.polies below)
#' @param clade.col a vector of length tip.groups with colours for each polytomy OR a single colour for all polytomies OR NULL (the default) to make all polytomies rainbow-coloured
#' @param br.clade.col the colours of the branches inside the cartoon. If NULL (the default, which you probably want), they don't show up in the plot.
#' @param auto.polies generate tip.groups by making a cartoon polytomy for all terminal polytomies. Setting this to TRUE (default is FALSE) will over-ride any tip.groups information. Note the function's return values is you want to set up a weird interaction of polytomies
#' @param ... additional arguments to pass to ape::plot.phylo
#' @details Makes a simple 'cartoon' phylogeny where polytomies are joined together in a big triangle. Could cause problems if you're outputting to PDF, sorry!
#' @return A list containing the tip.groups plotted, the clade colours they were plotted under, the edges that were obscured by the cartoon printing, and a list of the number of the node that each polytomy started under
#' @author Will Pearse
#' @examples
#' \dontrun{
#' require(ape)
#' tree <- read.tree(text="(((((A,B,C,D,E),(F,G,H,I,J)),H),K),L);")
#' cartoon.plot(tree, auto.polies=TRUE)
#' cartoon.plot(tree, list(1:5, 6:10), clade.col="grey30")
#' cartoon.plot(tree, list(1:5, 6:10), clade.col=c("blue", "red"))
#' }
#' @import ape
#' @import caper
#' @export
cartoon.plot <- function(tree, tip.groups=vector("list", 0), clade.col=NULL, br.clade.col=NULL, auto.polies=FALSE, ...){
    if(auto.polies == TRUE){
        node.table <- table(tree$edge[,1])
        polytomies <- as.numeric(names(node.table)[node.table>2])
        tip.groups <- vector("list", length(polytomies))
        for(i in seq_along(tip.groups)){
            t <- tree$edge[tree$edge[,1]==polytomies[i],2]
            #...can't handle non-terminal polytomies (including 'Adams consensus-style' polytomies, which need to be detected and removed at present...)
            if(all(t <= length(tree$tip.label)))
                tip.groups[[i]] <- t
        }
        lengths <- sapply(tip.groups, length)
        tip.groups <- tip.groups[lengths>0]
    }
    if(is.null(clade.col))
        clade.col <- rainbow(length(tip.groups))
    if(length(clade.col) < length(tip.groups))
        clade.col <- rep(clade.col[1], length(tip.groups))
    to.be.joined <- rep(FALSE, nrow(tree$edge))
    nodes <- vector("list", length(tip.groups))
    to.be.joined <- tree$edge[,2] %in% unlist(tip.groups)
    clade.mems <- lapply(seq(from=length(tree$tip.label)+1,to=max(tree$edge[,2])), clade.members, tree)
    for(i in seq_along(tip.groups)){
        nodes[[i]] <- which(sapply(clade.mems, function(x) identical(x,tip.groups[[i]])))+length(tree$tip.label)
    }
    to.be.joined <- to.be.joined | tree$edge[,1] %in% unlist(nodes)
    if(is.null(br.clade.col))
        plot(tree, edge.col=ifelse(to.be.joined, "white", "black"), ...) else plot(tree, plot=FALSE, ...)
    pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    for(i in seq_along(tip.groups)){
        range <- tip.groups[[i]][c(1,length(tip.groups[[i]]))]
        polygon(pp$xx[c(min(nodes[[i]]),range)], pp$yy[c(min(nodes[[i]]),range)], col=clade.col[i], border=clade.col[i])
    }
    if(!is.null(br.clade.col)){
        par(new=TRUE)
        plot(tree, edge.col=ifelse(to.be.joined, br.clade.col, "black"), ...)
    }
    invisible(list(tip.groups, clade.col, to.be.joined, nodes))
}

#' \code{ringlabels} Label particular tip(s) with text around the edge of a circular phylogeny
#' 
#' @param tip.groups list where each element are the tip *numbers* that should be labelled
#' @param text list when each element is the text to be plotted
#' @param radial.adj a multiplier for how far out each tip label should be
#' @param ... additional arguments for plotrix::arctext
#' @details Add text to the outside of a circular phylogeny. Useful if you've made a cartoon phylogeny and need to label clades.
#' @return The centers of each piece of text (in radians)
#' @author Will Pearse
#' @examples \dontrun{
#' tree <- read.tree(text="(((((A,B,C,D,E),(F,G,H,I,J)),H),K),L);")
#' ringlabels(tip.groups=list(1:5, 6:10) text=list("this is yet", "another test"))
#' tree <- read.tree(text="(((((A,B,C,D,E),(F,G,H,I,J)),H),K),L);")
#' }
#' @import ape
#' @import plotrix
#' @export
ringlabels <- function(tip.groups, text, radial.adj=1.05, ...){
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    raw <- seq(0, 2 * pi * (1 - 1/lastPP$Ntip) - 2 * pi * 1/360, length.out=lastPP$Ntip)
    edges <- lastPP$edge[,2]
    edges <- order(edges[edges <= lastPP$Ntip])
    radians <- numeric(length(tip.groups))
    for(i in seq_along(tip.groups)){
        radians[i] <- median(raw[edges %in% tip.groups[[i]]])
        arctext(x=text[[i]], radius=max(lastPP$xx)*radial.adj, middle=radians[i], ...)
    }
    invisible(radians[i])
}

tipring <- function(tips, col, radial.adj=1, ...){
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (missing(tips)) 
        tips <- seq(lastPP$Ntip)
    edges <- lastPP$edge[,2]
    edges <- order(edges[edges <= lastPP$Ntip])
    raw <- seq(0, 2 * pi * (1 - 1/lastPP$Ntip) - 2 * pi * 1/360, length.out=lastPP$Ntip)
    XX <- cos(raw) * max(lastPP$xx) * radial.adj
    YY <- sin(raw) * max(lastPP$xx) * radial.adj
    prev <- length(XX)
    coords <- matrix(NA, nrow=length(XX)+1, ncol=2)
    for(i in seq_along(XX)){
        x.adj <- (XX[i]-XX[prev])/2
        y.adj <- (YY[i]-YY[prev])/2
        coords[i,] <- c(XX[i]+x.adj, YY[i]+y.adj)
        prev <- i
    }
    coords[i+1,] <- coords[1,]
    for(i in seq_along(XX))
        if(i %in% tips)
            lines(c(coords[edges[i],1],coords[edges[i]+1,1]), c(coords[edges[i],2],coords[edges[i]+1,2]), col=col, ...)
    invisible(coords)
}

#' \code{willeerd.tiplabels} Plot tip labels with radial spacing
#' 
#' @param radial.adj A multiplier for how far out each tip label should be
#' @param ...everything else is exactly as tiplabels
#' @details A way of getting evenly spaced tip labels on a radial phylogeny. Just run the examples, it's quite simple. This is a very minorly-editted version of tiplabels; I can't take much credit! Don't cite this, cite ape!
#' @return ...exactly as tiplabels
#' @author Will Pearse
#' @examples \dontrun{
#' tree <- stree(128, type="balanced")
#' plot(tree, type="radial", show.tip.label=FALSE)
#' willeerd.tiplabels(tip=seq(128), pch=20)
#' willeerd.tiplabels(tip=seq(128), pch=20, radial.adj=1.05, col="red")
#' willeerd.tiplabels(tip=seq(128), pch=20, radial.adj=1.1, col="blue")
#' }
#' @import ape
#' @export
willeerd.tiplabels <- function (text, tip, adj = c(0.5, 0.5), radial.adj=1, frame = "rect", pch = NULL, thermo = NULL, pie = NULL, piecol = NULL, col = "black", bg = "yellow", horiz = FALSE, width = NULL, height = NULL, ...) 
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (missing(tip)) 
        tip <- 1:lastPP$Ntip
    lastPP$xx <- lastPP$xx * radial.adj
    lastPP$yy <- lastPP$yy * radial.adj
    XX <- lastPP$xx[tip]
    YY <- lastPP$yy[tip]
    BOTHlabels(text, tip, XX, YY, adj, frame, pch, thermo, pie, 
        piecol, col, bg, horiz, width, height, ...)
}

#' \code{willeerd.nodelabels} Plot tip labels with radial spacing
#' 
#' @param radial.adj A multiplier for how far out each node label should be
#' @param ...everything else is exactly as nodelabels
#' @details A way of getting evenly spaced node labels on a radial phylogeny. Just run the examples, it's quite simple. This is a very minorly-editted version of nodelabels; I can't take much credit! Don't cite this, cite ape!
#' @return ...exactly as nodelabels
#' @author Will Pearse
#' @examples \dontrun{
#' tree <- stree(128, type="balanced")
#' plot(tree, type="radial", show.tip.label=FALSE)
#' willeerd.nodelabels(pch=20)
#' willeerd.nodelabels(pch=20, radial.adj=1.05, col="red")
#' willeerd.nodelabels(pch=20, radial.adj=1.1, col="blue")
#' }
#' @import ape
#' @export
willeerd.nodelabels <- function (text, node, adj = c(0.5, 0.5), radial.adj=1, frame = "rect", pch = NULL, thermo = NULL, pie = NULL, piecol = NULL, col = "black", bg = "lightblue", horiz = FALSE, width = NULL, height = NULL, ...) 
{
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    if (missing(node)) 
        node <- (lastPP$Ntip + 1):length(lastPP$xx)
    lastPP$xx <- lastPP$xx * radial.adj
    lastPP$yy <- lastPP$yy * radial.adj
    XX <- lastPP$xx[node]
    YY <- lastPP$yy[node]
    willeerd.BOTHlabels(text, node, XX, YY, adj, frame, pch, thermo, pie, 
        piecol, col, bg, horiz, width, height, ...)
}

#' \code{willeerd.plot.phylo} Plot a rooted phylogeny (with more control over the root)
#' @details Almost *identical* to plot.phylo, but with root edge width control. Just try the example. Don't cite this, cite ape!
#' @return ...exactly as plot.phylo
#' @author Will Pearse
#' @examples \dontrun{
#' tree <- rtree(20)
#' par(mfrow=c(1,2))
#' plot(tree, edge.width=6, root.edge=TRUE)
#' willeerd.plot.phylo(tree, edge.width=6, root.edge=TRUE)
#' }
#' @import ape
#' @export
willeerd.plot.phylo <- function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
    show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
    edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
    adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
    label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
    direction = "rightwards", lab4ut = "horizontal", tip.color = "black", 
    plot = TRUE, rotate.tree = 0, open.angle = 0, ...) 
{
    Ntip <- length(x$tip.label)
    if (Ntip < 2) {
        warning("found less than 2 tips in the tree")
        return(NULL)
    }
    if (any(tabulate(x$edge[, 1]) == 1)) 
        stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles()")
    .nodeHeight <- function(Ntip, Nnode, edge, Nedge, yy) .C("node_height", 
        as.integer(Ntip), as.integer(Nnode), as.integer(edge[, 
            1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(yy), 
        DUP = FALSE, PACKAGE = "ape")[[6]]
    .nodeDepth <- function(Ntip, Nnode, edge, Nedge) .C("node_depth", 
        as.integer(Ntip), as.integer(Nnode), as.integer(edge[, 
            1]), as.integer(edge[, 2]), as.integer(Nedge), double(Ntip + 
            Nnode), DUP = FALSE, PACKAGE = "ape")[[6]]
    .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, 
        edge.length) .C("node_depth_edgelength", as.integer(Ntip), 
        as.integer(Nnode), as.integer(edge[, 1]), as.integer(edge[, 
            2]), as.integer(Nedge), as.double(edge.length), double(Ntip + 
            Nnode), DUP = FALSE, PACKAGE = "ape")[[7]]
    Nedge <- dim(x$edge)[1]
    Nnode <- x$Nnode
    ROOT <- Ntip + 1
    type <- match.arg(type, c("phylogram", "cladogram", "fan", 
        "unrooted", "radial"))
    direction <- match.arg(direction, c("rightwards", "leftwards", 
        "upwards", "downwards"))
    if (is.null(x$edge.length)) 
        use.edge.length <- FALSE
    if (type %in% c("unrooted", "radial") || !use.edge.length || 
        is.null(x$root.edge) || !x$root.edge) 
        root.edge <- FALSE
    if (type == "fan" && root.edge) {
        warning("drawing root edge with type = 'fan' is not yet supported")
        root.edge <- FALSE
    }
    phyloORclado <- type %in% c("phylogram", "cladogram")
    horizontal <- direction %in% c("rightwards", "leftwards")
    xe <- x$edge
    if (phyloORclado) {
        phyOrder <- attr(x, "order")
        if (is.null(phyOrder) || phyOrder != "cladewise") {
            x <- reorder(x)
            if (!identical(x$edge, xe)) {
                ereorder <- match(x$edge[, 2], xe[, 2])
                if (length(edge.color) > 1) {
                  edge.color <- rep(edge.color, length.out = Nedge)
                  edge.color <- edge.color[ereorder]
                }
                if (length(edge.width) > 1) {
                  edge.width <- rep(edge.width, length.out = Nedge)
                  edge.width <- edge.width[ereorder]
                }
                if (length(edge.lty) > 1) {
                  edge.lty <- rep(edge.lty, length.out = Nedge)
                  edge.lty <- edge.lty[ereorder]
                }
            }
        }
        yy <- numeric(Ntip + Nnode)
        TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
        yy[TIPS] <- 1:Ntip
    }
    z <- reorder(x, order = "pruningwise")
    if (phyloORclado) {
        if (is.null(node.pos)) {
            node.pos <- 1
            if (type == "cladogram" && !use.edge.length) 
                node.pos <- 2
        }
        if (node.pos == 1) 
            yy <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
        else {
            ans <- .C("node_height_clado", as.integer(Ntip), 
                as.integer(Nnode), as.integer(z$edge[, 1]), as.integer(z$edge[, 
                  2]), as.integer(Nedge), double(Ntip + Nnode), 
                as.double(yy), DUP = FALSE, PACKAGE = "ape")
            xx <- ans[[6]] - 1
            yy <- ans[[7]]
        }
        if (!use.edge.length) {
            if (node.pos != 2) 
                xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge) - 
                  1
            xx <- max(xx) - xx
        }
        else {
            xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
                z$edge.length)
        }
    }
    else {
        twopi <- 2 * pi
        rotate.tree <- twopi * rotate.tree/360
        switch(type, fan = {
            TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
            xx <- seq(0, twopi * (1 - 1/Ntip) - twopi * open.angle/360, 
                length.out = Ntip)
            theta <- double(Ntip)
            theta[TIPS] <- xx
            theta <- c(theta, numeric(Nnode))
            theta <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, 
                theta)
            if (use.edge.length) {
                r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, 
                  Nedge, z$edge.length)
            } else {
                r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
                r <- 1/r
            }
            theta <- theta + rotate.tree
            xx <- r * cos(theta)
            yy <- r * sin(theta)
        }, unrooted = {
            nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
            XY <- if (use.edge.length) unrooted.xy(Ntip, Nnode, 
                z$edge, z$edge.length, nb.sp, rotate.tree) else unrooted.xy(Ntip, 
                Nnode, z$edge, rep(1, Nedge), nb.sp, rotate.tree)
            xx <- XY$M[, 1] - min(XY$M[, 1])
            yy <- XY$M[, 2] - min(XY$M[, 2])
        }, radial = {
            X <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
            X[X == 1] <- 0
            X <- 1 - X/Ntip
            yy <- c((1:Ntip) * twopi/Ntip, rep(0, Nnode))
            Y <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
            xx <- X * cos(Y + rotate.tree)
            yy <- X * sin(Y + rotate.tree)
        })
    }
    if (phyloORclado) {
        if (!horizontal) {
            tmp <- yy
            yy <- xx
            xx <- tmp - min(tmp) + 1
        }
        if (root.edge) {
            if (direction == "rightwards") 
                xx <- xx + x$root.edge
            if (direction == "upwards") 
                yy <- yy + x$root.edge
        }
    }
    if (no.margin) 
        par(mai = rep(0, 4))
    if (is.null(x.lim)) {
        if (phyloORclado) {
            if (horizontal) {
                x.lim <- c(0, NA)
                pin1 <- par("pin")[1]
                strWi <- strwidth(x$tip.label, "inches")
                xx.tips <- xx[1:Ntip] * 1.04
                alp <- try(uniroot(function(a) max(a * xx.tips + 
                  strWi) - pin1, c(0, 1e+06))$root, silent = TRUE)
                if (is.character(alp)) 
                  tmp <- max(xx.tips) * 1.5
                else {
                  tmp <- if (show.tip.label) 
                    max(xx.tips + strWi/alp)
                  else max(xx.tips)
                }
                x.lim[2] <- tmp
            }
            else x.lim <- c(1, Ntip)
        }
        else switch(type, fan = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                x.lim <- c(min(xx) - offset, max(xx) + offset)
            } else x.lim <- c(min(xx), max(xx))
        }, unrooted = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                x.lim <- c(0 - offset, max(xx) + offset)
            } else x.lim <- c(0, max(xx))
        }, radial = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.03 * cex)
                x.lim <- c(-1 - offset, 1 + offset)
            } else x.lim <- c(-1, 1)
        })
    }
    else if (length(x.lim) == 1) {
        x.lim <- c(0, x.lim)
        if (phyloORclado && !horizontal) 
            x.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label) 
            x.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
                cex)
        if (type == "radial") 
            x.lim[1] <- if (show.tip.label) 
                -1 - max(nchar(x$tip.label) * 0.03 * cex)
            else -1
    }
    if (phyloORclado && direction == "leftwards") 
        xx <- x.lim[2] - xx
    if (is.null(y.lim)) {
        if (phyloORclado) {
            if (horizontal) 
                y.lim <- c(1, Ntip)
            else {
                y.lim <- c(0, NA)
                pin2 <- par("pin")[2]
                strWi <- strwidth(x$tip.label, "inches")
                yy.tips <- yy[1:Ntip] * 1.04
                alp <- try(uniroot(function(a) max(a * yy.tips + 
                  strWi) - pin2, c(0, 1e+06))$root, silent = TRUE)
                if (is.character(alp)) 
                  tmp <- max(yy.tips) * 1.5
                else {
                  tmp <- if (show.tip.label) 
                    max(yy.tips + strWi/alp)
                  else max(yy.tips)
                }
                y.lim[2] <- tmp
            }
        }
        else switch(type, fan = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                y.lim <- c(min(yy) - offset, max(yy) + offset)
            } else y.lim <- c(min(yy), max(yy))
        }, unrooted = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
                y.lim <- c(0 - offset, max(yy) + offset)
            } else y.lim <- c(0, max(yy))
        }, radial = {
            if (show.tip.label) {
                offset <- max(nchar(x$tip.label) * 0.03 * cex)
                y.lim <- c(-1 - offset, 1 + offset)
            } else y.lim <- c(-1, 1)
        })
    }
    else if (length(y.lim) == 1) {
        y.lim <- c(0, y.lim)
        if (phyloORclado && horizontal) 
            y.lim[1] <- 1
        if (type %in% c("fan", "unrooted") && show.tip.label) 
            y.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
                cex)
        if (type == "radial") 
            y.lim[1] <- if (show.tip.label) 
                -1 - max(nchar(x$tip.label) * 0.018 * max(yy) * 
                  cex)
            else -1
    }
    if (phyloORclado && direction == "downwards") 
        yy <- max(yy) - yy
    if (phyloORclado && root.edge) {
        if (direction == "leftwards") 
            x.lim[2] <- x.lim[2] + x$root.edge
        if (direction == "downwards") 
            y.lim[2] <- y.lim[2] + x$root.edge
    }
    asp <- if (type %in% c("fan", "radial", "unrooted")) 
        1
    else NA
    plot(0, type = "n", xlim = x.lim, ylim = y.lim, ann = FALSE, 
        axes = FALSE, asp = asp, ...)
    if (plot) {
        if (is.null(adj)) 
            adj <- if (phyloORclado && direction == "leftwards") 
                1
            else 0
        if (phyloORclado && show.tip.label) {
            MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
            loy <- 0
            if (direction == "rightwards") {
                lox <- label.offset + MAXSTRING * 1.05 * adj
            }
            if (direction == "leftwards") {
                lox <- -label.offset - MAXSTRING * 1.05 * (1 - 
                  adj)
            }
            if (!horizontal) {
                psr <- par("usr")
                MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - 
                  psr[1])
                loy <- label.offset + MAXSTRING * 1.05 * adj
                lox <- 0
                srt <- 90 + srt
                if (direction == "downwards") {
                  loy <- -loy
                  srt <- 180 + srt
                }
            }
        }
        if (type == "phylogram") {
            phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal, 
                edge.color, edge.width, edge.lty)
        }
        else {
            if (type == "fan") {
                ereorder <- match(z$edge[, 2], x$edge[, 2])
                if (length(edge.color) > 1) {
                  edge.color <- rep(edge.color, length.out = Nedge)
                  edge.color <- edge.color[ereorder]
                }
                if (length(edge.width) > 1) {
                  edge.width <- rep(edge.width, length.out = Nedge)
                  edge.width <- edge.width[ereorder]
                }
                if (length(edge.lty) > 1) {
                  edge.lty <- rep(edge.lty, length.out = Nedge)
                  edge.lty <- edge.lty[ereorder]
                }
                circular.plot(z$edge, Ntip, Nnode, xx, yy, theta, 
                  r, edge.color, edge.width, edge.lty)
            }
            else cladogram.plot(x$edge, xx, yy, edge.color, edge.width, 
                edge.lty)
        }
        if (root.edge) 
            switch(direction, rightwards = segments(0, yy[ROOT], 
                x$root.edge, yy[ROOT], lwd = edge.width), leftwards = segments(xx[ROOT], 
                yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT], lwd = edge.width), 
                upwards = segments(xx[ROOT], 0, xx[ROOT], x$root.edge, 
                  lwd = edge.width), downwards = segments(xx[ROOT], 
                  yy[ROOT], xx[ROOT], yy[ROOT] + x$root.edge, 
                  lwd = edge.width))
        if (show.tip.label) {
            if (is.expression(x$tip.label)) 
                underscore <- TRUE
            if (!underscore) 
                x$tip.label <- gsub("_", " ", x$tip.label)
            if (phyloORclado) 
                text(xx[1:Ntip] + lox, yy[1:Ntip] + loy, x$tip.label, 
                  adj = adj, font = font, srt = srt, cex = cex, 
                  col = tip.color)
            if (type == "unrooted") {
                if (lab4ut == "horizontal") {
                  y.adj <- x.adj <- numeric(Ntip)
                  sel <- abs(XY$axe) > 0.75 * pi
                  x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
                    1.05
                  sel <- abs(XY$axe) > pi/4 & abs(XY$axe) < 0.75 * 
                    pi
                  x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
                    (2 * abs(XY$axe)[sel]/pi - 0.5)
                  sel <- XY$axe > pi/4 & XY$axe < 0.75 * pi
                  y.adj[sel] <- strheight(x$tip.label)[sel]/2
                  sel <- XY$axe < -pi/4 & XY$axe > -0.75 * pi
                  y.adj[sel] <- -strheight(x$tip.label)[sel] * 
                    0.75
                  text(xx[1:Ntip] + x.adj * cex, yy[1:Ntip] + 
                    y.adj * cex, x$tip.label, adj = c(adj, 0), 
                    font = font, srt = srt, cex = cex, col = tip.color)
                }
                else {
                  adj <- abs(XY$axe) > pi/2
                  srt <- 180 * XY$axe/pi
                  srt[adj] <- srt[adj] - 180
                  adj <- as.numeric(adj)
                  xx.tips <- xx[1:Ntip]
                  yy.tips <- yy[1:Ntip]
                  if (label.offset) {
                    xx.tips <- xx.tips + label.offset * cos(XY$axe)
                    yy.tips <- yy.tips + label.offset * sin(XY$axe)
                  }
                  font <- rep(font, length.out = Ntip)
                  tip.color <- rep(tip.color, length.out = Ntip)
                  cex <- rep(cex, length.out = Ntip)
                  for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                    cex = cex[i], x$tip.label[i], adj = adj[i], 
                    font = font[i], srt = srt[i], col = tip.color[i])
                }
            }
            if (type %in% c("fan", "radial")) {
                xx.tips <- xx[1:Ntip]
                yy.tips <- yy[1:Ntip]
                angle <- atan2(yy.tips, xx.tips)
                if (label.offset) {
                  xx.tips <- xx.tips + label.offset * cos(angle)
                  yy.tips <- yy.tips + label.offset * sin(angle)
                }
                s <- xx.tips < 0
                angle <- angle * 180/pi
                angle[s] <- angle[s] + 180
                adj <- as.numeric(s)
                font <- rep(font, length.out = Ntip)
                tip.color <- rep(tip.color, length.out = Ntip)
                cex <- rep(cex, length.out = Ntip)
                for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                  x$tip.label[i], font = font[i], cex = cex[i], 
                  srt = angle[i], adj = adj[i], col = tip.color[i])
            }
        }
        if (show.node.label) 
            text(xx[ROOT:length(xx)] + label.offset, yy[ROOT:length(yy)], 
                x$node.label, adj = adj, font = font, srt = srt, 
                cex = cex)
    }
    L <- list(type = type, use.edge.length = use.edge.length, 
        node.pos = node.pos, show.tip.label = show.tip.label, 
        show.node.label = show.node.label, font = font, cex = cex, 
        adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
        x.lim = x.lim, y.lim = y.lim, direction = direction, 
        tip.color = tip.color, Ntip = Ntip, Nnode = Nnode)
    assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)), 
        envir = .PlotPhyloEnv)
    invisible(L)
}

#' \code{willeerd.BOTHlabels} Internal function for plotting tiplabels and nodelabels
#' @details Almost identical to the similar function in ape! This is an internal function
#' @author Will Pearse
#' @import ape
willeerd.BOTHlabels <- function (text, sel, XX, YY, adj, frame, pch, thermo, pie, piecol, 
    col, bg, horiz, width, height, ...) 
{
    if (missing(text)) 
        text <- NULL
    if (length(adj) == 1) 
        adj <- c(adj, 0.5)
    if (is.null(text) && is.null(pch) && is.null(thermo) && is.null(pie)) 
        text <- as.character(sel)
    frame <- match.arg(frame, c("rect", "circle", "none"))
    args <- list(...)
    CEX <- if ("cex" %in% names(args)) 
        args$cex
    else par("cex")
    if (frame != "none" && !is.null(text)) {
        if (frame == "rect") {
            width <- strwidth(text, units = "inches", cex = CEX)
            height <- strheight(text, units = "inches", cex = CEX)
            if ("srt" %in% names(args)) {
                args$srt <- args$srt%%360
                if (args$srt == 90 || args$srt == 270) {
                  tmp <- width
                  width <- height
                  height <- tmp
                }
                else if (args$srt != 0) 
                  warning("only right angle rotation of frame is supported;\n         try  `frame = \"n\"' instead.\n")
            }
            width <- xinch(width)
            height <- yinch(height)
            xl <- XX - width * adj[1] - xinch(0.03)
            xr <- xl + width + xinch(0.03)
            yb <- YY - height * adj[2] - yinch(0.02)
            yt <- yb + height + yinch(0.05)
            rect(xl, yb, xr, yt, col = bg, border=bg)
        }
        if (frame == "circle") {
            radii <- 0.8 * apply(cbind(strheight(text, units = "inches", 
                cex = CEX), strwidth(text, units = "inches", 
                cex = CEX)), 1, max)
            symbols(XX, YY, circles = radii, inches = max(radii), 
                add = TRUE, bg = bg)
        }
    }
    if (!is.null(thermo)) {
        parusr <- par("usr")
        if (is.null(width)) {
            width <- CEX * (parusr[2] - parusr[1])
            width <- if (horiz) 
                width/15
            else width/40
        }
        if (is.null(height)) {
            height <- CEX * (parusr[4] - parusr[3])
            height <- if (horiz) 
                height/40
            else height/15
        }
        if (is.vector(thermo)) 
            thermo <- cbind(thermo, 1 - thermo)
        thermo <- if (horiz) 
            width * thermo
        else height * thermo
        if (is.null(piecol)) 
            piecol <- rainbow(ncol(thermo))
        xl <- XX - width/2 + adj[1] - 0.5
        xr <- xl + width
        yb <- YY - height/2 + adj[2] - 0.5
        yt <- yb + height
        if (horiz) {
            rect(xl, yb, xl + thermo[, 1], yt, border = NA, col = piecol[1])
            for (i in 2:ncol(thermo)) rect(xl + rowSums(thermo[, 
                1:(i - 1), drop = FALSE]), yb, xl + rowSums(thermo[, 
                1:i]), yt, border = NA, col = piecol[i])
        }
        else {
            rect(xl, yb, xr, yb + thermo[, 1], border = NA, col = piecol[1])
            for (i in 2:ncol(thermo)) rect(xl, yb + rowSums(thermo[, 
                1:(i - 1), drop = FALSE]), xr, yb + rowSums(thermo[, 
                1:i]), border = NA, col = piecol[i])
        }
        s <- apply(thermo, 1, function(xx) any(is.na(xx)))
        xl[s] <- xr[s] <- NA
        rect(xl, yb, xr, yt, border = "black")
        if (!horiz) {
            segments(xl, YY, xl - width/5, YY)
            segments(xr, YY, xr + width/5, YY)
        }
    }
    if (!is.null(pie)) {
        if (is.vector(pie)) 
            pie <- cbind(pie, 1 - pie)
        xrad <- CEX * diff(par("usr")[1:2])/50
        xrad <- rep(xrad, length(sel))
        XX <- XX + adj[1] - 0.5
        YY <- YY + adj[2] - 0.5
        for (i in seq_along(sel)) {
            if (any(is.na(pie[i, ]))) 
                next
            floating.pie.asp(XX[i], YY[i], pie[i, ], radius = xrad[i], 
                col = piecol)
        }
    }
    if (!is.null(text)) 
        text(XX, YY, text, adj = adj, col = col, ...)
    if (!is.null(pch)) 
        points(XX + adj[1] - 0.5, YY + adj[2] - 0.5, pch = pch, 
            col = col, bg = bg, ...)
}

#' \code{factorise.tree} 'Factorise' a tree by removing out species, making it easier to plot/manipulate
#' 
#' @param tree ape::phylo phylogeny to be 'factorised'
#' @param scale.factor multiplier for the number of species within a terminal polytomy. E.g., 0.1 means each terminal polytomy will be ~10% its current size
#' @details Thins out additional species, making a phylogeny smaller by reducing the size of each terminal polytomy by scale.factor
#' @return List where first element is the factorised phylogeny, the second the tips that were dropped from each node (on the original tree; I can't guarantee the returning tree's structure)
#' @author Will Pearse
#' @examples \dontrun{
#' tree <- read.tree(text="((A,B,C,D,E),F);")
#' t <- factorise.tree(tree, 0.5)
#' plot(t$tree)
#' }
#' @import ape
#' @export
factorise.tree <- function(tree, scale.factor=0.5){
    #Get the terminal nodes
    nodes <- table(tree$edge[tree$edge[,2]<=length(tree$tip.label),1])
    nodes <- nodes[nodes > 2]

    #Reduce the diversity of those nodes and collect tips to drop
    nodes <- round(nodes * scale.factor)
    x <- 1
    to.drop <- numeric(length(tree$tip.label))
    for(i in seq_along(nodes)){
        prog.bar(i,length(nodes))
        t <- tree$edge[tree$edge[,1]==names(nodes[i]),2]
        #Not all "terminal" polytomies have only tips descending (outgroups!)
        t <- t[t <= length(tree$tip.label)]
        t <- t[-1:-nodes[i]]
        to.drop[seq(from=x, length.out=length(t))] <- t
        x <- x+length(t)
    }
    
    #Drop those tips and return
    tree <- drop.tip(tree, unique(to.drop))
    return(list(tree=tree, dropped=to.drop))
}
t2 <- factorise.tree(tpl, 0.001)
