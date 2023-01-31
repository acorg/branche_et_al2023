inf_code_shapes <- c("inf" = 24, "non_inf" = 21)

ag_order <- c("D614G" = 1, "B.1.617.2" = 3, "B.1.351" = 5, "BA.1" = 7, "BA.4/5" = 9, "BA.2.12.1" = 11)

ag_plot_names <- c("D614G" = "D614G", "B.1.617.2" = "Delta", "B.1.351" = "Beta", "BA.1" = "BA.1", "BA.2.12.1" = "BA.2.12.1", "BA.4/5" = "BA.4/BA.5")

arm_code_order_no_vacc <- c("P", "P+O","O", "D+O", "B+O", "B+O, B+O",
                            "P+B", "B")
arm_code_order <- c("P:M","P+O:M", "O:M", "D+O:M", "B+O:M", "B+O, B+O:M",
                    "P:Pf","P+O:Pf", "O:Pf", "P+B:Pf", "B+O:Pf", "B:Pf",
                    "P:S", "P+B:S", "B:S",
                    arm_code_order_no_vacc)

arm_plot_names <- c("Prototype", "Prototype + Omicron BA.1","Omicron BA.1", "Delta + Omicron BA.1", "Beta + Omicron BA.1", "2x(Beta + Omicron BA.1)", #Moderna
                    "Prototype", "Prototype + Omicron BA.1", "Omicron BA.1", "Prototype + Beta", "Beta + Omicron BA.1", "Beta", #Pfizer
                    "Prototype", "Prototype + Beta", "Beta", #Sanofi
                    "Prototype", "Prototype + Omicron BA.1","Omicron BA.1", "Delta + Omicron BA.1", "Beta + Omicron BA.1", "2x(Beta + Omicron BA.1)", "Prototype + Beta", "Beta",
                    "D29-D1", "D91-D1", "D91-D29", "D1", "D29", "D91", "D15", "D57", "Moderna", "Pfizer", "Sanofi") #no vacc manuf
names(arm_plot_names) <- c(arm_code_order, "D29  D1", "D91  D1", "D91  D29", "D1", "D29", "D91", "D15", "D57", "M", "Pf", "S")

inf_code_order <- c("inf", "non_inf")
age_code_order <- c("<65", ">65", "combined")
visit_code_order <- c("D1","D15", "D29", "D57","D85","D91")

line_vals <- c("<65" = "solid", ">65" = "solid", "combined" = "solid", "NA" = "solid",
"inf" = "solid", "non_inf" = "solid", "V3" = "solid", "V4" = "solid", "V1" = "solid",
"M" = "solid", "Pf" = "dashed", "S" = "dotted",
"D1" = "solid", "D29"= "solid", "D15"= "solid", "D91"= "solid", "D85"= "solid", "D57"= "solid")

point_spacing <- 0.2

`+.uneval` <- function(a,b) {
  `class<-`(modifyList(a,b), "uneval")
}

x_by_visit_diff <- function(data){
  
  visits <- unlist(lapply(data$visit_diff, function(x){
    d <- strsplit(x, " - ")[[1]][1]
    d <- gsub("D", "")
    as.numeric(d)
  }))
  
  data$x <- visits

  return(data)
}

x_position_spacing <- function(data, position_by = "arm_code"){
  
  if(position_by == "visit_diff"){
    return(x_by_visit_diff(data))
  }
  unique_labels <- unique(data %>% pull(position_by))
  
  unique_pos <- length(unique_labels)
  total_len <- point_spacing*(unique_pos-1)
  
  seq_spaces <- seq(from = -total_len/2, to = total_len/2, by = point_spacing)
  names(seq_spaces) <- unique_labels
  
  data$x <- unlist(lapply(1:nrow(data), function(x){
    target <- as.character(data[[x, position_by]])
    ag_order[data$ag_name[x]] + seq_spaces[target]
  }))
  
  
  return(data)
  
}
titerlineplot_dodge <- function(data, sr_group_colors, titer_thresh = 13, antigens = c("D614G","B.1.351", "B.1.617.2", "BA.1"),
                          facet_n_row = 1, sr_group_order = NULL, gmt_facetter = "visit_code", nrow_gmt = 1, 
                          color_by = "sr_group", shape_by = "inf_code", line_by = "age_code", x_position_by = "age_code",
                          cols_to_keep = c("arm_code", "age_code", "inf_code"), to_long = T, 
                          show_gmt_label = F, show_group_count = F, show_mean_line = F, mean_line_color = "black",
                          gmt_grid_row = "v_manuf_code", gmt_grid_col = "arm_code", dodge_group = "visit_code", pre_titer_adj = FALSE) {
  
  ag_order <- ag_order[antigens]
  # format titer to longer
  if(to_long){
    data %>%
      pivot_longer(cols = antigens, names_to = "ag_name", values_to = "titer") %>%
      mutate(logtiter = log2(titer/10),
             logtiter = ifelse(titer < titer_thresh, log2((titer_thresh/20)), logtiter)) -> data_long
  } else{
    data_long <- data 
  }
  
  # create plot colors
  
  sr_groups <- unique(as.character(data_long$sr_group))
  
  if(length(color_by) == 1){
    if(color_by == "lab_code"){
    plot_colors <- c("Monogram" = "#ed94a7", "Montefiori" = "#6C992E", "Adjusted" = "#bb94f1",
      "Monogram-PNEUT50" = "#ed94a7", "Monogram-PNEUT80" = "#ed94a7", "Montefiori-PNEUT50" = "#6C992E", "Montefiori-PNEUT80" = "#9bd748")
    } else {

      plot_colors <- unlist(lapply(sr_groups, function(x) {
        blend_sr_group_colors(x, sr_group_colors, split_by = color_by)
      }))
  
      sr_group_colors[sr_groups, "Color"] <- plot_colors
  
      plot_colors <- sr_group_colors$Color
  
      names(plot_colors)<- rownames(sr_group_colors)

    }
  } else {
      plot_colors <- unlist(lapply(sr_groups, function(x) {
        blend_sr_group_colors(x, sr_group_colors, split_by = color_by)
      }))
  
      sr_group_colors[sr_groups, "Color"] <- plot_colors
  
      plot_colors <- sr_group_colors$Color
  
      names(plot_colors)<- rownames(sr_group_colors)

  }

  
  
  # plot_colors <-  plot_colors[as.character(unique(data_long[, color_by]))]
  
  # set sr_order
  if(!is.null(sr_group_order)) {
    data_long$sr_group <- factor(as.character(data_long$sr_group), levels =sr_group_order)
  }
  
  
  #  max_y <- max(ceiling(data_long$logtiter), na.rm = T) +0.5
  max_y <- 15
  
  
  sr_group_gmt_plotdata <- sr_group_gmt_calc(data_long, titer_thresh, cols_to_keep = cols_to_keep, pre_adj = pre_titer_adj) %>%
    filter(!all_below_thresh) %>%
    calculate_fold_change_from_ag() %>%
    mutate(y = max_y,
           label = paste0(round(as.numeric(titer), 0), "\n", fold_change))

  # sr_group_gmt_plotdata <- sr_group_to_column(sr_group_gmt_plotdata)
  
  
  # add other factor levels for arm code
  data_long$arm_code <- factor(data_long$arm_code, levels = arm_code_order)
  
  sr_group_gmt_plotdata$arm_code <- factor(sr_group_gmt_plotdata$arm_code, levels = arm_code_order)
  
  # order by arm code
  data_long <- data_long[order(data_long$arm_code),]
  sr_group_gmt_plotdata <- sr_group_gmt_plotdata[order(sr_group_gmt_plotdata$arm_code),]
  
  data_long <- x_position_spacing(data_long, position_by = x_position_by)
  sr_group_gmt_plotdata <- x_position_spacing(sr_group_gmt_plotdata, position_by = x_position_by)
  # 
  # fill  NA values of lower and upper CI with mean value
  sr_group_gmt_plotdata$lower[is.na(sr_group_gmt_plotdata$lower)] <- sr_group_gmt_plotdata$logtiter[is.na(sr_group_gmt_plotdata$lower)]
  sr_group_gmt_plotdata$upper[is.na(sr_group_gmt_plotdata$upper)] <- sr_group_gmt_plotdata$logtiter[is.na(sr_group_gmt_plotdata$upper)]
  
  facet_labeller <- function(x){
    x <- gsub("NA-", "", x)
    x <- gsub("NA|<65|>65|non_inf|V3|V1|inf", "", x)
    x <- gsub("-", " ", x)
    x <- gsub("B\\+O, B\\+O", "2x(B+O)", x)
    x <- arm_plot_names[x]
    x
  }
  
  # facet_labeller <- function(x){
  #   strsplit(x, "-")[[1]][2]
  # }
  sub_visit <- FALSE
  if(length(color_by) > 1){
    if("visit_code" %in% color_by | length(unique(data$visit_code)) > 1){
      sub_visit <- TRUE
    }
    color_by <- "sr_group"
  }
  
  data_long$visit_code <- factor(data_long$visit_code, levels = visit_code_order)
  sr_group_gmt_plotdata$visit_code <- factor(sr_group_gmt_plotdata$visit_code, levels = visit_code_order)
  print(unique(sr_group_gmt_plotdata$visit_code))
  data_long %>%
    ggplot(
      aes_string(
        x = "ag_name",
        y = "logtiter",
        color = color_by,
        fill = color_by,
        shape = shape_by
      )
    ) + 
    geom_line(
      aes(group = sr_name),
      linetype = "solid",
      alpha = 0.4
    ) + 
    geom_point(
      alpha = 0.4
    ) +
    geom_line(
      data = sr_group_gmt_plotdata,
      aes_string(group = "sr_group",
                 linetype = line_by),
      size = 1.3
    ) +
    geom_pointrange(
      data = sr_group_gmt_plotdata,
      aes(ymin = lower,
          ymax = upper)
    ) +
    geom_point(
      data = sr_group_gmt_plotdata,
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
    geom_text(data = sr_group_gmt_plotdata,
              mapping = aes(x = ag_name, y = y, label = label),
              color = "black",
              size = 2.5
    ) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    scale_linetype_manual(values= line_vals) +
    guides(linetype="none") + 
    scale_y_titer(
      ymin = log2(titer_thresh/10)
    ) + 
    scale_x_discrete(
      limits = antigens,
      labels = ag_plot_names[antigens],
      name = "Antigen variant"
    ) + 
    coord_cartesian(
      ylim = shrinkrange(c(-0.5, max_y+1), 0.05)
    ) +
    titerplot_theme() + 
    theme(
      legend.position = "none",
      strip.text.x = element_text(size = 10)
    ) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = -1,
      ymax = log2(titer_thresh/10),
      fill = "grey50",
      color = NA,
      alpha = 0.3
    ) -> gp
  
  
  visits <- unique(sr_group_gmt_plotdata$visit_code)
  
  if(line_by != "v_manuf_code" & "V3" %in% visits & "V4" %in% visits){
    line_by <- "visit_code"
  }
  # plot gmts of different serum groups together
  max_y <- 13
  sr_group_gmt_plotdata$y <- max_y #- 0.5
  
  if(length(unique(sr_group_gmt_plotdata$visit_code)) <=4 & dodge_group == "visit_code"){
    sr_group_gmt_plotdata$visit_code <- factor(sr_group_gmt_plotdata$visit_code, levels = c("D29", "D1", "D91", "D15"))
  }

  if(dodge_group == "arm_code"){
    sr_group_gmt_plotdata$arm_code <- factor(as.character(sr_group_gmt_plotdata$arm_code), levels = c("P+O", "D+O", "P", "B+O","O", "B+O, B+O", "P+B", "B"))
  }

  p_arm <- sr_group_gmt_plotdata %>% filter(arm_code == "P") %>% filter(visit_code == "D1")
  
  sr_group_gmt_plotdata %>%
    ggplot(
      aes_string(
        x = "x",
        y = "logtiter",
        color = color_by,
        fill = color_by,
        shape = shape_by
      )
    ) + 
    geom_line(
      aes_string(group = "sr_group",
                 linetype = line_by),
      size = 1
    ) +
    # pointrange if I want them dodged, x position by. If I just want it dodged then I can add position = position_dodge(width = 1)
    # geom_pointrange(
    #   aes(ymin = lower,
    #       ymax = upper)
    # ) +
    # errorbar if I want shifted bars
    geom_errorbar(
      aes_string(ymin = "lower",
          ymax = "upper",
          group = dodge_group), # usually visit code
      width = 0,
      position = position_dodge(width =0.4)
    ) +
    geom_point(
      size = 2
    ) +
    geom_point(
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
    geom_line(
      data = p_arm,
      aes(group = sr_group,
                 linetype = inf_code),
      size = 1
    ) +
     geom_point(
       data = p_arm,
      size = 2
    ) +
    geom_point(
      data = p_arm,
      size = 1.5,
      fill = "white",
      shape = 21
    ) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    scale_linetype_manual(values= line_vals) +
    guides(linetype="none") + 
    scale_y_titer(
      ymin = log2(titer_thresh/10)
    ) + 
    scale_x_continuous(
      breaks = ag_order,
      labels = ag_plot_names[names(ag_order)],
      limits = c(min(data_long$x, na.rm = T)-0.5, max(max(data_long$x, na.rm = T), max(ag_order, na.rm = T))+0.5),
      name = "Antigen variant"
    ) + 
    coord_cartesian(
      ylim = shrinkrange(c(-0.5, max_y+0.5), 0.05)
    ) +
    titerplot_theme() + 
    theme(
      legend.position = "top",
      plot.title = element_text(size = 6),
      strip.text.x = element_text(size = 11),
      strip.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 10)
    ) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = -1,
      ymax = log2(titer_thresh/10),
      fill = "grey50",
      color = NA,
      alpha = 0.3
    ) -> gp_gmt
  
  
  if(length(unique(sr_group_gmt_plotdata$v_manuf_code)) > 1 | gmt_grid_row != "v_manuf_code"){
    gp <- gp + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col)),
        labeller = as_labeller(facet_labeller)
      ) 
    gp_gmt <- gp_gmt + 
      facet_grid(
        as.formula(paste(gmt_grid_row, "~", gmt_grid_col)),
        labeller = as_labeller(facet_labeller)
      ) 
    
    count_label <- sr_group_gmt_plotdata %>% 
      ungroup() %>%
      select(all_of(c(gmt_facetter,"sr_group", "arm_code", "v_manuf_code", "count", "y", "age_code", color_by, shape_by))) %>%
      mutate(x = max(sr_group_gmt_plotdata$x),
             label = paste0("n = ", count)) %>%
      unique()
    
    if(sub_visit){
      visits <- unique(sr_group_gmt_plotdata$visit_code)
      
      target_visit <- visits[visits != "D15"]
      
      count_label <- count_label[grepl(target_visit[1], count_label$sr_group),]
    }
    
    count_label <- count_label %>%
      group_by(sr_group,arm_code, v_manuf_code, age_code) %>%
      summarize(total_count = sum(count),
                label = paste0("n (", age_code, ") = ", total_count),
                x = x,
                y = ifelse(age_code == "combined", y, ifelse(age_code == "<65", y - 1, y-2)),
                inf_code = inf_code)
    
    count_label <- count_label %>%
      mutate(y = ifelse(inf_code == "inf", y-0.5, y-1))
    
    count_label$label <- gsub("\\(combined\\)", "", count_label$label)
    
    
    #  count_label <- count_label[!duplicated(count_label$arm_code),]
    
    
  } else {
    gp <- gp + 
      facet_wrap(
        vars(sr_group),
        labeller = as_labeller(facet_labeller),
        nrow = facet_n_row
      ) 
    
    gp_gmt <- gp_gmt + facet_wrap(
      as.formula(paste("~", gmt_facetter)),
      labeller = as_labeller(facet_labeller),
      nrow = nrow_gmt
    ) 
    
    count_label <- sr_group_gmt_plotdata %>% 
      ungroup() %>%
      select(all_of(c(gmt_facetter,"sr_group", "arm_code", "count", "y", "age_code", color_by, shape_by))) %>%
      mutate(x = max(sr_group_gmt_plotdata$x),
             label = paste0("n = ", count)) %>%
      unique()
    
    if(sub_visit){
      visits <- unique(sr_group_gmt_plotdata$visit_code)
      
      target_visit <- visits[!("D15" == visits)]
      
      count_label <- count_label[grepl(target_visit[1], count_label$sr_group),]
    }
    
    count_label <- count_label %>%
      group_by(sr_group, arm_code, age_code) %>%
      summarize(total_count = sum(count),
                label = paste0("n (", age_code, ") = ", total_count),
                x = x,
                y = ifelse(age_code == "combined", y, ifelse(age_code == "<65", y - 1, y-2)),
                inf_code = inf_code)
    
    count_label <- count_label %>%
      mutate(y = ifelse(inf_code == "inf", y, y-1))
    
    count_label$label <- gsub("\\(combined\\)", "", count_label$label)
    
    # count_label <- count_label[!duplicated(count_label$arm_code),]
  }
  
  if(shape_by == "inf_code"){
    gp <- gp + 
      scale_shape_manual(values = inf_code_shapes)
    
    gp_gmt <- gp_gmt + 
      scale_shape_manual(values = inf_code_shapes)
    
  }
  
  if(show_gmt_label){
    gp_gmt <- gp_gmt +
      geom_text(data = sr_group_gmt_plotdata,
                #mapping = aes(x = ag_order[ag_name], y = ifelse(age_code == "<65", y, y-1.5), label = label),
                mapping = aes(x = ag_order[ag_name], y = y, label = label),
                color = "black",
                size = 2.5
      )
  }
  
  if(show_group_count){
    
    gp_gmt <- gp_gmt +
      geom_text(data = count_label,
                #mapping = aes(x = ag_order[ag_name], y = ifelse(age_code == "<65", y, y-1.5), label = label),
                mapping = aes(x = x, y = y+0.5, label = label, hjust = 1),
                color = "black",
                size = 4
      )
  }
  
  if(show_mean_line){
   
    gp <- gp +
      geom_line(aes_string(y = "gmt_arm", group = "sr_group", linetype = line_by), color = mean_line_color, 
                size = 1.2, alpha = 0.3)
    
    gp_gmt <- gp_gmt + 
      geom_line(aes(x = ag_order[ag_name], y = gmt_arm, group = sr_group)+aes_string(linetype = line_by), 
                color = mean_line_color, size = 1.1, alpha = 0.3)
    
  }
  
  return(list("all" = gp, "gmt" = gp_gmt))
  
  
}

