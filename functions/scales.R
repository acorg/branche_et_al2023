
shrinkrange <- function(x, amount) {
  
  xn <- c(amount / 2, 1 - amount / 2)
  xn * diff(range(x)) + x[1]
  
}

scale_x_titer <- function(
    threshold = "<20",
    logthreshold = 0,
    axisname = "Titer"
) {
  
  scale_x_continuous(
    name = axisname,
    breaks = function(x) {
      logthreshold:ceiling(max(x))
    },
    labels = function(x) {
      output <- 2^x*10
      output[x == logthreshold] <- threshold
      output
    },
    minor_breaks = NULL
  )
  
}

scientific_10 <- function(x) {
  #parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
  expression(paste("10"^x))
}

scale_y_titer <- function(
    threshold = "<20",
    log10scale = TRUE,
    logthreshold = 0,
    axisname = "Titer",
    ymin = NULL,
    ymax = NULL,
    ...
) {
  
  if(!log10scale){
     scale_y_continuous(
    name = axisname,
    breaks = function(x) {
      if (is.null(ymin)) ymin <- ceiling(min(x))
      if (is.null(ymax)) ymax <- floor(max(x))
      ymin:ymax
    },
    labels = function(x) {
      output <- 2^x*10
      output[x == logthreshold] <- threshold
      output
    },
    minor_breaks = NULL,
    ...
  )
  } else {
     scale_y_continuous(
    name = axisname,
    breaks = c(logthreshold, 1:4),
    labels = function(x) {
      output <- parse(text = paste0("10^",x+1))
      output[x == logthreshold] <- threshold
      output
    },
    minor_breaks = NULL,
    ...
  )
  }
 
  
}

titerplot_theme <- function(){
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.background = element_rect(
      fill = NA,
      colour = "grey40"
    ),
    panel.grid.major = element_line(
      linetype = "solid",
      colour = rgb(0,0,0,0.05)
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      # face = "bold",
      size = 10
    )
  )
  
}