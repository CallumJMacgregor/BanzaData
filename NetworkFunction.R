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




specieslev <- function(x) {
  p <- c("bipartite")                                  # list of necessary packages
  packages <- p[!(p %in% installed.packages()[,"Package"])]      # checks each is installed
  
  if(length(packages)) {
    stop("Install required packages - ",packages)                # returns an error saying which packages are missing
    
  } else {
    lapply(p, require, character.only = TRUE)                    # loads up packages if all are installed
  }
  
  if (nrow(x)>2) {                   # checks if more than two insect species were sampled
    sample <- x[1,2]
    x <- x[c(1,3:length(x))]                # remove the Sample column (every entry is identical within each dframe)
    rownames(x) <- x[,1]             # set the row names as the first column (Family_Species)
    x <- x[,-1]                      # remove the first column, leaving Family_Species as the row names only
    x <- t(x)                        # transpose the matrix so pollinators are in columns
    data <- specieslevel(x, index = "ALLBUTD")  # produces species-level metrics for the matrix, as a list of two dframes
    higher.data <- data.frame(data[1])  # take only higher-level dframe
    lower.data <- data.frame(data[2])   # take only lower-level dframe
    
    namesh <- colnames(higher.data)    # extract higher-level metric names
    namesh <- gsub('higher.level.','',namesh)   # generalise higher-level metric names
    
    namesl <- colnames(lower.data)    # extract lower-level metric names
    namesl <- gsub('lower.level.','',namesl)   # generalise lower-level metric names
    
    
    colnames(higher.data) <- namesh    # replace metric names with generalised versions
    colnames(lower.data) <- namesl
  
    lower.data <- lower.data[,colnames(higher.data)] # reorder lower-level metrics to be same order as higher-level
      
    higher.data$level <- "insect"        # create a variable indicating whether a species is an insect or a plant
    lower.data$level <- "plant"
  
    higher.data$Species <- rownames(higher.data)   # creates a column for species
    
    lower.data$Species <- rownames(lower.data)     # creates a column for species
    lower.data$Species <- gsub(paste(sample[1],'',sep="."),'',lower.data$Species)  # removes information about season in species column
        
    combined <- rbind(higher.data,lower.data)  # merge the insects and plants together
    
    combined$sample <- sample[1]   # create a variable indicating what season the metric refers to
    
    
    return(combined)                                        # outputs the data
    
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