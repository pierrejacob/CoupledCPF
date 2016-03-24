#'@rdname setmytheme
#'@title Set My Theme
#'@description set theme of plots
#'@export
setmytheme <- function(){
  library(ggplot2)
  library(gridExtra)
  library(reshape)
  theme_set(theme_bw())
  theme_update(axis.text.x = element_text(size = 20),
               axis.text.y = element_text(size = 20),
               axis.title.x = element_text(size = 25),
               axis.title.y = element_text(size = 25, angle= 90),
               legend.text = element_text(size = 20),
               legend.title = element_text(size = 20),
               title = element_text(size = 30),
               strip.text = element_text(size = 25),
               strip.background = element_rect(fill="white"),
               legend.position = "bottom")
}