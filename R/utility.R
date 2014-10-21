prog.bar <- function(x, y)
    tryCatch(if((x %% floor(y/100)) == 0) if((x %% (floor(y/100)*10)) == 0) cat("|") else cat("."), error=function(z) cat("."))

#' \code{lean} Local Environmental Adaptation Null model
#' 
#' @param union Allow for juvenile adaptation
#' @details This is extremely unlikely to be of any interest to anyone
#' other than those attending the Silwood Park Christmas 2010 meeting,
#' where the importance of LEAN was disucssed in relation to student
#' teacher. Seriously, don't even bother running it.
#' @return Nothing; plot-type
#' @author Will Pearse
#' @examples \dontrun{
#' lean()
#' lean(FALSE)
#' }
#' @export
lean <- function(union=FALSE)
    if(union==TRUE) browseURL("http://youtu.be/XWqpXhWa3ug") else browseURL("http://youtu.be/2H4fk97IKnA")
