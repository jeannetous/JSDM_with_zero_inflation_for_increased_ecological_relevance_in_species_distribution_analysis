sigmoid <- function(x){ 1 / (1 + exp(-x))}

logit <- function(p){ log(p/(1 - p))}
