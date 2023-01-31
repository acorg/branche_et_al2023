padding <- 0.25

screenshot_html_landscape_to_png <- function(landscape, save_name ){

  to_save <- file.path(paste0(save_name, ".html"))
  saveWidget(landscape, to_save, selfcontained = FALSE)

}

plot_single_landscape_panel <- function(landscape, label, label_size = 10, label_x_pos = 2, label_y_pos = 9,
                                        sr_group_label = "", sr_group_y_pos = 0, sr_group_size = 3, show_border = FALSE,
                                        delete_html = TRUE, save_name = "temp"){
  
  
  to_save <- file.path(paste0(save_name, ".html"))
  png_save <- gsub(".html", ".png", to_save)
  saveWidget(landscape, to_save, selfcontained = FALSE)
  webshot(url=to_save,file = png_save)
  temp_plot <- readPNG(png_save)
  
  qplot(c(1:10),c(1:10), geom="blank") +
    annotation_custom(rasterGrob(temp_plot, height = unit(0.7, "npc")), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    annotate(geom="text", x=label_x_pos, y=label_y_pos, label=label,size= label_size, hjust = 0) + 
    annotate(geom="text", x=label_x_pos, y=sr_group_y_pos, label=sr_group_label,size= sr_group_size, hjust = 0) +
    theme_void() -> p

if(show_border) {
    p + theme(panel.border = element_rect(color = "grey50",
                                      fill = NA,
                                      size = 0.5))-> p
}
  
if(delete_html){
  if (file.exists(to_save)) {
    #Delete file if it exists
    file.remove(to_save)
  }
  if (file.exists(png_save)) {
    #Delete file if it exists
    file.remove(png_save)
  }
}
  
 return(p) 
}


# sams landscape functions to add landscape from lndscp fits list
get_titertable <- function(data, group) {
  
  data %>% 
    select(
      variant,
      subjid,
      titer
    ) %>%
    mutate(
      titer = replace(titer, is.na(titer), "*")
    ) %>%
    pivot_wider(
      id_cols = subjid,
      names_from = variant,
      values_from = titer
    ) %>% 
    as.matrix() -> titermatrix
  
  attr(titermatrix, "arm") <- group$arm_code
  attr(titermatrix, "visit") <- group$visit_code
  rownames(titermatrix) <- titermatrix[,"subjid"]
  titermatrix <- titermatrix[,-1]
 
  return(titermatrix)
  
}

# Functions to remove buttons
addObject3js <- function(
    data3js,
    object,
    number_of_ids = 1
){
  
  # Generate an object ID
  if(is.null(data3js$lastID)){ data3js$lastID <- 0 }
  object$ID <- max(data3js$lastID) + seq_len(number_of_ids)
  
  # If object is interactive and highlighted add a reference to itself to
  # it's highlight group by default
  if(!is.null(object$properties$interactive)){
    object$group <- object$ID
  }
  
  # Add the object to the plot data
  data3js$plot[[length(data3js$plot)+1]] <- object
  
  # Update the ID of the last object added
  data3js$lastID <- object$ID
  
  # Return the new data
  data3js
  
}

remove_buttons <- function(data3js){
  
  new_data3js = data3js
  
  new_data3js = data3js
  
  new_data3js[['lastID']] = 0
  new_data3js[['plot']] = list()
  
  N = data3js[['lastID']] 
  
  
  
  
  for (i in 1:N)
  {
    obj = data3js[['plot']][[i]]
    
    
    
    if ('toggle' %in% names(obj[['properties']])){
      obj[['properties']][['toggle']] <- NULL
    }
    
    new_data3js = addObject3js(new_data3js,obj)
    
    
  }
  
  
  
  return (new_data3js)
  
}


base_plot_data3js <- function(map, lndscp_fits, highlighted_ags, lims, ag_plot_names, alternative_ba5 = FALSE, opti_nr = 1,
                              add_border = TRUE, add_axis = TRUE){
  
  if(alternative_ba5){
    # set points and coordinates of highlighted ags
    x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1], agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 1])
    y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2], agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 2])
    z_coords <- rep(0.02, length(highlighted_ags))
    ag_point_size <- c(rep(14, length(highlighted_ags)), 12) / 5
    ag_col <- c(agOutline(map)[agNames(map) %in% highlighted_ags], agOutline(map)[agNames(map) =="BA.4/BA.5"])
    ag_fill <- c(agFill(map)[agNames(map) %in% highlighted_ags], agFill(map)[agNames(map) =="BA.4/BA.5"])
    labels <- c(ag_plot_names[agNames(map) %in% highlighted_ags], "BA.4/BA.5(2)")
    
  } else {
    
    x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1])
    x_coords["BA.4/BA.5"]<- agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 1]
    y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2])
    y_coords["BA.4/BA.5"]<- agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 2]
    z_coords <- rep(0.02, length(highlighted_ags)-1)
    ag_point_size <- c(rep(14, length(highlighted_ags)-1)) / 5
    ag_col <- c(agOutline(map)[agNames(map) %in% highlighted_ags])
    ag_fill <- c(agFill(map)[agNames(map) %in% highlighted_ags])
    labels <- c(ag_plot_names[agNames(map) %in% highlighted_ags])
    
  }
  
  border_col <- "grey50"
  
  z_lims <- c(0,12)
  axis_at <- seq(z_lims[1], z_lims[2],2)
  # Setup plot
  data3js <- ablandscapes:::lndscp3d_setup(
    xlim = lims$xlim,
    ylim = lims$ylim,
    zlim = z_lims,
    aspect.z = 0.5,
    options = list(
      lwd.grid =  0.05,
      sidegrid.lwd = 1,
      sidegrid.col = border_col,
      sidegrid.at = list("z" = axis_at),
      zaxt = "log"
    ),
    show.axis = FALSE
  )
  
  if(add_axis){

    axis_labels <- 2^axis_at*10
    
    data3js <- r3js::axis3js(
      data3js,
      side = "z",
      at = axis_at,
      labels = axis_labels,
    # labeloffset = 0.11,
      cornerside = "f",
      size = 30,
      alignment = "right"
    )
  }

  # Add basemap
  data3js <- lndscp3d_map(
    data3js = data3js,
    fit = lndscp_fits[[1]],
    xlim = lims$xlim,
    ylim = lims$ylim,
    zlim = c(0, 10),
    show.map.sera = FALSE,
    options = list(
      opacity.basemap = 0.3
    )
  )
  
  data3js <- r3js::points3js(
    data3js,
    x          = x_coords,
    y          = y_coords,
    z          = z_coords,
    size       = ag_point_size,
    col        = ag_col,
    fill       = ag_fill,
    lwd        = 0.5,
    opacity    = 1,
    highlight  = list(col = "red"),
    label      = labels,
    toggle     = "Basepoints",
    depthWrite = FALSE,
    shape      = "circle filled"
  )
  
  if(add_border){
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[1]), y = c(lims$ylim[1], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    data3js <- lines3js(data3js, x = c(lims$xlim[2],lims$xlim[2]), y = c(lims$ylim[1], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    
    # y border
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[2]), y = c(lims$ylim[1], lims$ylim[1]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[2]), y = c(lims$ylim[2], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)

    data3js <- r3js::box3js(
      data3js,
      col   = border_col
    )
    
  }
  
  return(data3js)
}


plot_landscapes_from_list <- function(data3js, titertables_groups, lndscp_fits,map, gmt_data, highlighted_ags,
                                      ag_plot_names, alternative_ba5 = FALSE, opti_nr = 1, hide_buttons = TRUE){
  
  if(alternative_ba5){
    x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1], agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 1])
    y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2], agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 2])
    z_coords <- rep(0.02, length(highlighted_ags))
    ag_point_size <- c(rep(14, length(highlighted_ags)), 12) / 5
   # text_x <- c(c(x_coords[1:4] + ag_point_size[1:4]*0.15),c(x_coords[5:6] - ag_point_size[5:6]*0.25))
  #  text_y <- c(y_coords[1:4], c(y_coords[5:6] - ag_point_size[5:6]*0.2))
    text_x <- c(x_coords[1:6] - ag_point_size[1:6]*0.2)
    text_y <- c(y_coords[1:6] - ag_point_size[1:6]*0.2)
    text_plot <- c(ag_plot_names[agNames(map)[agNames(map) %in% highlighted_ags]], "BA.4/BA.5(2)")
    
  } else {
    x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1])
    x_coords["BA.4/BA.5"]<- agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 1]
    y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2])
    y_coords["BA.4/BA.5"]<- agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 2]
    z_coords <- rep(0.02, length(highlighted_ags)-1)
    ag_point_size <- c(rep(14, length(highlighted_ags)-1)) / 5
  #  text_x <- c(c(x_coords[1:3] + ag_point_size[1:3]*0.15),c(x_coords[4:5] - ag_point_size[4:5]*0.2))
  #  text_y <- c(y_coords[1:3], c(y_coords[4:5] - ag_point_size[4:5]*0.2))
    text_x <- c(x_coords - ag_point_size*0.3)
    text_y <- c(y_coords - ag_point_size*0.2)
    text_plot <- c(ag_plot_names[agNames(map)[agNames(map) %in% highlighted_ags]])
    
  }
  
    min_offset <- -0.1
    max_offset <- 0.1
    offset <- seq(from = min_offset, to = max_offset, by = (max_offset - min_offset)/length(lndscp_fits))
    
  for (i in seq_along(lndscp_fits)) {
    
   # message(i)
    arm <- titertables_groups$arm_code[i]
    visit <- titertables_groups$visit_code[i]
    v_manuf <- titertables_groups$v_manuf_code[i]
    lndscp_fit <- lndscp_fits[[i]]
    
    coords <- cbind(x_coords, y_coords)
    
    # Add titers
    gmts <- filter(gmt_data, visit_code == visit, arm_code == arm, v_manuf_code == v_manuf)
    
    if (is.na(filter(gmts, variant == "BA.4/BA.5")$gmt)){
      gmts <- gmts[!is.na(gmts$gmt),]
      coords <- coords[rownames(coords) != "BA.4/BA.5",]
    } # next
    
    gmts <- gmts[match(rownames(coords), gmts$variant),]
    
    for (j in seq_len(nrow(coords))) {
      
      data3js <- r3js::lines3js(
        data3js,
        x = rep(coords[j, 1], 2),
        y = rep(coords[j, 2], 2),
        z = c(0, gmts$gmt[j]),
        col = "grey50",
        highlight = list(col = "red"),
        interactive = FALSE,
     #   toggle = sprintf("Titers, %s, %s", arm, visit),
        geometry = TRUE,
        lwd = 0.4 #was 1
      )

     # print(coords[j,])
     # print(offset[j])
     # print(gmts$gmt[j])
      
      data3js <- r3js::points3js(
        data3js,
        x         = coords[j, 1] + offset[i],
        y         = coords[j, 2],
        z         = gmts$gmt[j],
      #  size      = 0.7, #was 2
        size      = 0.9, #was 2
      #  col       = "grey50",
        col  = arm_cols[arm],
        highlight = list(col = "red"),
   #     label     = gmts$variant[j],
        #   toggle    = sprintf("Titers, %s, %s", arm, visit),
        opacity   = 1 # was 1
     
      )
      
     # text_x <- c(c(agCoords(map)[agNames(map) %in% highlighted_ags[1:4], 1] + agSize(map)[agNames(map) %in% highlighted_ags[1:4]]*0.025), c(agCoords(map)[agNames(map) == "BA.4/BA.5", 1]- agSize(map)[agNames(map) == "BA.4/BA.5"]*0.07))
    #  text_y <- c(c(agCoords(map)[agNames(map) %in% highlighted_ags[1:4], 2]), c(agCoords(map)[agNames(map) == "BA.4/BA.5", 2] - agSize(map)[agNames(map) == "BA.4/BA.5"]*0.07))
      
      # set points and coordinates of highlighted ags

      
      data3js <- r3js::text3js(
        data3js,
        x          = text_x,
        y          = text_y,
        z          = z_coords,
        text       = text_plot,
        toggle     = "Labels",
        size       = c(rep(22*0.02, length(text_x))), #agSize(map)[agNames(map) %in% highlighted_ags]*0.02,
        alignment  = "right"
      )
      
    }
    
    # Add landscapes
    data3js <- lndscp3d_surface(
      data3js = data3js,
      object = lndscp_fit,
      # zlim = c(0, 10),
      crop2chull = FALSE,
      # crop2base = TRUE,
      toggle = sprintf("Landscape, %s, %s", arm, visit),
      grid_spacing = 0.5,
      padding = padding,
      options = list(
        col.surface = arm_cols[arm],
       # opacity.surface = 0.5
        opacity.surface = 1
      )
    )
    
  }
  
  if(hide_buttons){
    data3js <- remove_buttons(data3js)
  }
 
  
  
  return(data3js)
}




