
test <- data.frame(dframes[1])
test <- test[c(1,3:75)]                # remove the Sample column (every entry is identical within each dframe)
rownames(test) <- test[,1]             # set the row names as the first column (Family_Species)
test <- test[,-1]                      # remove the first column, leaving Family_Species as the row names only
test <- t(test)                        # transpose the matrix so pollinators are in columns
data <- data.frame(networklevel(test, index = "binary"))  # produces network metrics for the matrix
data[is.na(data)] <- 0                                 # makes NAs into 0s
data$metric <- row.names(data)
secex <- second.extinct(test, participant = "both", method = "random")   # simulates secondary extinctions
robust <- data.frame(robustness(secex))               # calculates robustness from secex simulations
robust$metric <- row.names(robust)
robust$metric <- "robustness"
row.names(robust) <- robust$metric
names(robust) <- names(data)
results <- rbind(data,robust)
rownames(results) <- results$metric
results <- results[c(0:1)]

summary(results)


tmetrics <- metrics
