#'@rdname normalize_weight
#'@title normalize_weight
#'@description takes log weights and return normalized weights
#'@export
normalize_weight <- function(logw){
  w <- exp(logw - max(logw))
  w <- w / sum(w)
  return(w)
}
