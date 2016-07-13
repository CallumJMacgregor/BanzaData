sacplot <- function(x) {
  p <- c("vegan")
  packages <- p[!(p %in% installed.packages()[,"Package"])]
  if(length(packages)) {
    stop("Install required packages - ",packages)
  } else {
      lapply(p, require, character.only = TRUE)
  }
  x <- x[,-c(1:9)]
  sac <- specaccum(x)
  plot(sac, ci.type="polygon",ci.col="yellow")
}


samplebased <- function(x) {
  p <- c("vegan")
  packages <- p[!(p %in% installed.packages()[,"Package"])]
  if(length(packages)) {
    stop("Install required packages - ",packages)
  } else {
    lapply(p, require, character.only = TRUE)
  }
  x <- x[,-c(1:9)]
  SR <- estimateR(x)
  return(SR)
}


sitebased <- function(x) {
  p <- c("vegan")
  packages <- p[!(p %in% installed.packages()[,"Package"])]
  if(length(packages)) {
    stop("Install required packages - ",packages)
  } else {
    lapply(p, require, character.only = TRUE)
  }
  x <- x[,-c(1:9)]
  SR <- t(specpool(x))
  return(SR)
}
