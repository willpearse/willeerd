#' Load US NCDC data (climate, etc.)
#'
#' @param file location (on hard drive or URL) of NCDC file to be loaded
#' @param tempout name of temporary file to be written out (and then
#' deleted) as part of read process. If this file already exists,
#' it'll get overwritten and then deleted!
#' @examples
#' \dontrun{
#' #The location/name of this file changes as it's updated, so below is only a guide...
#' palmer <- read.ncdc("http://www1.ncdc.noaa.gov/pub/data/cirs/climdiv/climdiv-pdsidv-v1.0.0-20160105")
#' }
#' @export
read.ncdc <- function(file, tempout="delete.me.txt"){
    #Process and load raw file
    data <- read.delim(file, header=FALSE, as.is=TRUE)[,1]
    data <- gsub("[ ]+", " ", data)
    write.table(data, tempout, row.names=FALSE, col.names=FALSE, quote=FALSE)
    data <- read.table(tempout, header=FALSE, as.is=TRUE, colClasses="character")
    unlink(tempout)
    
    #Match up meta-data
    states <- setNames(c("Alabama","Arizona","Arkansas","California","Colorado","Connecticut","Delaware","Florida","Georgia","Idaho","Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana","Maine","Maryland","Massachusetts","Michigan","Minnesota","Mississippi","Missouri","Montana","Nebraska","Nevada","New Hampshire","New Jersey","New Mexico","New York","North Carolina","North Dakota","Ohio","Oklahoma","Oregon","Pennsylvania","Rhode Island","South Carolina","South Dakota","Tennessee","Texas","Utah","Vermont","Virginia","Washington","West Virginia","Wisconsin","Wyoming"), formatC(1:48, 1, flag="0"))[substr(data[,1], 1, 2)]
    division <- substr(data[,1], 3, 4)
    types <- setNames(c("Precipitation","Average Temperature","PDSI","PHDI","ZNDX","PMDI","Heating Degree Days","Cooling Degree Days","Maximum Temperature","Minimum Temperature","1-month Standardized Precipitation Index","2-month Standardized Precipitation Index","3-month Standardized Precipitation Index","6-month Standardized Precipitation Index","9-month Standardized Precipitation Index","12-month Standardized Precipitation Index","24-month Standardized Precipitation Index"),c("01","02","05","06","07","08","25","26","27","28","71","72","73","74","75","76","77"))[substr(data[,1], 5, 6)]
    year <- substr(data[,1], 7, 8)
    
    #Neaten and return
    data <- data.frame(states, division, types, year, data[,-1])
    names(data)[-1:-4] <- c("january","february","march","april","may","june","july","august","september","october","november","december")
    return(data)
}
