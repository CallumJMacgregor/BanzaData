network <- function(x) {
  p <- c("bipartite")                                  # list of necessary packages
  packages <- p[!(p %in% installed.packages()[,"Package"])]      # checks each is installed
  
  if(length(packages)) {
    stop("Install required packages - ",packages)                # returns an error saying which packages are missing
    
  } else {
    lapply(p, require, character.only = TRUE)                    # loads up packages if all are installed
  }
  
  if (nrow(x)>2) {                   # checks if more than two insect species were sampled
    x <- x[c(1,3:length(x))]                # remove the Sample column (every entry is identical within each dframe)
    rownames(x) <- x[,1]             # set the row names as the first column (Family_Species)
    x <- x[,-1]                      # remove the first column, leaving Family_Species as the row names only
    x <- t(x)                        # transpose the matrix so pollinators are in columns
    data <- data.frame(networklevel(x, index = "ALLBUTDD"))  # produces network metrics for the matrix
    data[is.na(data)] <- 0                                 # makes NAs into 0s
    data$metric <- row.names(data)                         # creates a column containing the metric names
    secex <- second.extinct(x, participant = "both", method = "random")   # simulates secondary extinctions
    robust <- data.frame(robustness(secex))                # calculates robustness from secex simulations
    robust$metric <- row.names(robust)                     # creates a column for the metric names
    robust$metric <- "robustness"                          # sets the metric name as 'robustness'
    row.names(robust) <- robust$metric                     # makes 'robustness' the row name
    names(robust) <- names(data)                           # makes the column names of the two dframes identical
    results <- rbind(data,robust)                          # tacks the robustness value on the end of the other metrics
    rownames(results) <- results$metric                    # resets the rownames with the full list of metric names
    results <- results[c(0:1)]                             # trims off the metric names column, preserving the rownames
    return(results)                                        # outputs the data
    
  } else {                                               # if only one or two insect species sampled...
    
    warning("Not enough data to create network")         # returns a warning, and...
    print("Fail")                                        # prints "Fail" as a character to the output list
    
  }
}




prepare <- function(x) {
  x <- x[c(1,3:length(x))]                # remove the Sample column (every entry is identical within each dframe)
  rownames(x) <- x[,1]             # set the row names as the first column (Family_Species)
  x <- x[,-1]                      # remove the first column, leaving Family_Species as the row names only
  x <- t(x)                        # transpose the matrix so pollinators are in columns
  return(x)                        # outputs the data
}