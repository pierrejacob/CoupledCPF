#'@export

hilbert_order <- function(x){
  return(hilbert_order_(t(x)) + 1)
}