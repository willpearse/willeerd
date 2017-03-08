#' Predict missing traits assuming Brownian motion
#' 
#' @param tree a phylogeny
#' @param trait a vector of traits, with NAs for missing values
#' @importFrom phytools phylosig
#' @importFrom ape vcv
#' @importFrom geiger rescale.phylo
#' @examples
#' real <- pred <- numeric(200)
#' for(i in seq_along(real)){
#'     tree <- sim.bdtree()
#'     trait <- setNames(sim.char(tree, .01)[,,1], tree$tip.label)
#'     real[i] <- trait[50]; trait[50] <- NA
#'     pred[i] <- predict.trait(trait, tree)
#' }
#' plot(pred ~ real)
#' @export
predict.trait <- function(trait, tree){
    if(is.null(attr(trait,"names")))
        stop("Requires named trait to work")
    got.data <- tree$tip.label[!is.na(trait)]
    to.predict <- tree$tip.label[is.na(trait)]

    signal <- phylosig(tree, trait, method="lambda")
    tree <- rescale(tree, "lambda", signal$lambda)
    vcv <- vcv(tree)

    predictions <- setNames(numeric(length(to.predict)), to.predict)
    for(i in seq_along(predictions))
        predictions[i] <- weighted.mean(trait[got.data], vcv[!rownames(vcv) %in% to.predict[i], to.predict[i]])
    return(predictions)
}


