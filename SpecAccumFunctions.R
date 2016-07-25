sacplot <- function(x,cols) {
  p <- c("vegan")
  packages <- p[!(p %in% installed.packages()[,"Package"])]
  if(length(packages)) {
    stop("Install required packages - ",packages)
  } else {
      lapply(p, require, character.only = TRUE)
  }
  x <- x[,-c(1:cols)]
  sac <- specaccum(x)
  plot(sac, ci.type="polygon",ci.col="yellow")
}


samplebased <- function(x,cols) {
  p <- c("vegan")
  packages <- p[!(p %in% installed.packages()[,"Package"])]
  if(length(packages)) {
    stop("Install required packages - ",packages)
  } else {
    lapply(p, require, character.only = TRUE)
  }
  x <- x[,-c(1:cols)]
  SR <- estimateR(x)
  return(SR)
}


sitebased <- function(x,cols) {
  p <- c("vegan")
  packages <- p[!(p %in% installed.packages()[,"Package"])]
  if(length(packages)) {
    stop("Install required packages - ",packages)
  } else {
    lapply(p, require, character.only = TRUE)
  }
  x <- x[,-c(1:cols)]
  SR <- t(specpool(x))
  return(SR)
}
