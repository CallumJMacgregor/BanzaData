network <- function(x) {
  p <- c("bipartite")                                  # list of necessary packages
  packages <- p[!(p %in% installed.packages()[,"Package"])]      # checks each is installed
  
  if(length(packages)) {
    stop("Install required packages - ",packages)                # returns an error saying which packages are missing
  
    } else {
    lapply(p, require, character.only = TRUE)                    # loads up packages if all are installed
    }
  
  if (nrow(x)>2) {                   # checks if more than two insect species were sampled
    x <- x[c(1,3:75)]                # remove the Sample column (every entry is identical within each dframe)
    rownames(x) <- x[,1]             # set the row names as the first column (Family_Species)
    x <- x[,-1]                      # remove the first column, leaving Family_Species as the row names only
    x <- t(x)                        # transpose the matrix so pollinators are in columns
    data <- data.frame(networklevel(x, index = "ALLBUTDD"))  # produces network metrics for the matrix
    data[is.na(data)] <- 0                                 # makes NAs into 0s
    secex <- second.extinct(x, participant = "both", method = "random")   # simulates secondary extinctions
    robust <- data.frame(robustness(secex))               # calculates robustness from secex simulations
    return(data)                                          # outputs the data
    return(robust)
    
  } else {                                               # if only one or two insect species sampled...
    
    warning("Not enough data to create network")         # returns a warning, and...
    print("Fail")                                        # prints "Fail" as a character to the output list
    
    }
}