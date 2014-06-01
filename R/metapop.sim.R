#' \code{sim.meta.comm} simulate a (sort of) meta-community
#' 
#' @param size the length and width of the meta-community in grid cells
#' @param n.spp number of species
#' @param timesteps number of time-steps (each discrete)
#' @param p.migrate probability that a group of species in each grid cell will migrate to another grid cell each timestep (i.e., 10 cells occuped by species A --> 10*p.migrate chance of migration)
#' @param env.lam Lambda value for Poisson distribution used to distribute environmental quality; essentially the carrying capacity (for each species separately) for that cell
#' @param abund.lam Lambda value for Poisson distribution used to distribute initial abundances and abundance after migration
#' @param stoch.lam Lambda value for Poisson distribution of noise added to the next step abundance calculation. With equal chance, this is taken as either a positive or a negative number (see details if you're confused as to why this is Poisson!)
#' @details Simulates species moving through a metacommunity. At each time-step each cell's next abundance for each species is env.quality - current.abundance + stochastic, and a species gets as many chances to migrate in each time-step as it has cells (the same cell could migrate multiple times). I use a Poisson for everything because I don't want half-species (these are individuals, because I eventually want to make a phylogeny of these), and keeping everything in Poisson makes it easier to compare the relative rates of everything.
#' @return List with the species abundances (as a 3D array) and the environmental quality (carrying capacities)
#' @author Will Pearse
#' @examples \dontrun{
#' tree <- sim.bd.tree(0.1, 0, 10)
#' plot(tree)
#' }
#' @export
sim.meta.comm <- function(size=10, n.spp=8, timesteps=10, p.migrate=0.05, env.lam=10, abund.lam=5, stoch.lam=1){
    #Setup environment and abundances
    env <- matrix(rpois(size^2, env.lam), nrow=size, ncol=size)
    abundance <- array(rpois(size^2*n.spp, abund.lam), dim=c(size,size,n.spp))
    
    #Loop over for timesteps
    for(i in seq(timesteps)){
        #Loop over species
        for(j in seq(n.spp)){
            #Calculate new abundance in each cell
            present <- abundance[,,j]>0
            stoch <- rpois(1,stoch.lam) * sample(c(1,-1),1)
            abundance[,,j][present] <- env[present] - abundance[,,j][present] + stoch
            #Migration
            # - randomly choose number of migration events, then choose cells for it to happen in (for simulation ease)
            # - remember the species may have died out...
            present <- abundance[,,j]>0
            if(sum(present) > 0){
                n.migrate <- sum(rbinom(sum(present),1,p.migrate))
                cells <- which(present, arr.ind=TRUE)
                for(k in seq(n.migrate)){
                    index <- cells[sample(seq(nrow(cells)), 1),]
                    index[1] <- index[1] + sample(c(-1,0,1),1)
                    index[2] <- index[2] + sample(c(-1,0,1),1)
                    #Not a wrapped world; you can fall off the edge!
                    if(all(index<=size))
                        abundance[index] <- env[index] - rpois(1,abund.lam)
                }
            }
        }
        #Clean up negative species etc. (probably not necessary...)
        abundance[abundance < 0] <- 0
    }
    
    #Return
    return(list(species=abundance, environment=env))
}

