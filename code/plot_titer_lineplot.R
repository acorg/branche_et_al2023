#  Do line plots per serum group
# do B+O, B+O arm v8 as day 29 (= equal to V4 of other arms)
# and day 57 as day 1, as this corresponds to the pre of the 2nd B+O dose
# drop normal day 29 of B+O, B+O arm V4
# and make copy of V8 and put it as V5 (day 91) as it is last recorded time point
rm(list = ls())
library(meantiter)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(dplyr)


source(file.path("functions", "format_data.R"))
source(file.path("functions","gmt_calculation.R"))
source(file.path("functions","gmt_fold_change.R"))
source(file.path("functions","scales.R"))
source(file.path("functions","titerlineplot.R"))
source(file.path("functions","sr_group_color_functions.R"))


correct_magnitude <- F

date <- "14DEC2022"
day <- "91"

inf_names <- c("inf" = "inf", non_inf = "non_inf")

BO_D_ajudstment <- c("wo_2dose")
day_visno <- c("D1","D29", "D91")

data <- read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Moderna_D", day,".csv")), sep= ",") %>%
    plyr::rbind.fill(read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Pfizer_D", day,".csv")), sep= ",")) %>%
    plyr::rbind.fill(read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Sanofi_D", day,".csv")), sep= ",")) 
  
  plot_antigens <- c("D614G", "B.1.617.2", "B.1.351", "BA.1", "BA.4/5")

# format the data
data <- format_data(data, sr_group_code, combine_b_o_arm = FALSE, remove_all_bfl = FALSE, remove_all_subj_entries = FALSE)

data_b <- data

for(d_adj in BO_D_ajudstment){
  
  # now drop B+O, B+O arm day 29 make copy of B+O arm day 
  data <- BO_format_visit_code(data_b, d_adj)
  
  
  
  # filter out non responders
  non_responders <- non_responders <- data %>% filter(D614G == 20) %>% filter(B.1.617.2 == 20) %>%
    filter(BA.1 == 20) %>% filter(B.1.351 == 20) %>% filter(`BA.4/5` == 20) %>% pull(subjid)
  
  data <- data %>% filter(!(subjid %in% non_responders))
  warning(paste0(length(unique(non_responders)), " subjects were filtered out because they were non responders"))
  
  # add info to group by
  # default grouping is by study arm 
  sr_group_data <- group_data_sr_group(data, by_arm = T, by_visit = T, by_strata = F, by_age = F, by_infection = T)
    
  sr_group_data$age_code <- "combined"
  sr_group_data_b <- sr_group_data
  
  # cycle through infected & non infected
  for(inf_stat in inf_names){
    
  figure_dir <- file.path("figures",date, inf_names[inf_stat], "titer_lineplots", paste0(d_adj,"_b+o_adj"))
      suppressWarnings(dir.create(figure_dir, recursive = TRUE))
    # combine sr groups by age and and filter for inf_stat
    sr_group_data <- sr_group_data_b %>% 
      filter(inf_code == inf_stat)
    
    sr_group_data$arm_code <- factor(sr_group_data$arm_code, levels = arm_code_order)
    sr_group_data$age_code <- factor(sr_group_data$age_code, levels = age_code_order)
    sr_group_data$inf_code <- factor(sr_group_data$inf_code, levels = inf_code_order)
    
    # subset to idvls for which we have V3 data
    # v3_ids <- sr_group_data %>% filter(visit_code == day_visno[day]) %>% pull(subjid) %>% unique()
    
    # sr_group_data <- sr_group_data %>% filter(subjid %in% v3_ids)
    # create sr_group order for idvl panels
    sr_groups <- unique(sr_group_data$sr_group)
    sr_groups <- c(sr_groups[grepl("-inf", sr_groups)], sr_groups[grepl("non_inf", sr_groups)])
    
    sr_order_df <- sr_group_data[order(sr_group_data$arm_code, sr_group_data$age_code),]
    sr_order_arm_age <- unique(sr_order_df$sr_group)
    
    sr_order_arm_age <- c(sr_order_arm_age[grepl("<65", sr_order_arm_age)], sr_order_arm_age[grepl(">65", sr_order_arm_age)],
                          sr_order_arm_age[grepl(paste0("NA-", inf_names[inf_stat]), sr_order_arm_age)])
    
    
    # do big plot with facet of D1, D29, D91
    sr_group_order <- sr_order_arm_age
    
    unique_arms <- length(unique(sr_group_data$arm_code))
    nrow_facet <- 4
    plot_width <- 0.5+2*unique_arms
    
    # ---------------------------- V1, V3 comparison -----------------------------------
    
    sr_groups <- unique(sr_group_data$sr_group)
    
   # add p arm to superimpose 
    sr_group_data %>% 
      data_to_long(antigens =  plot_antigens) %>%
      arm_code_wo_vacc() %>%
      group_by(ag_name, visit_code, v_manuf_code) %>%
      mutate(gmt_arm = meantiter::mean_titers(titer[arm_code == "P"], method ="bayesian", dilution_stepsize = 0)$mean)->sr_group_data_long
    
    # x_position_by = "visit_code",
    plot_height <- ifelse(length(unique(sr_group_data_long$v_manuf_code)) > 1, 2.5*length(unique(sr_group_data_long$v_manuf_code)), 3.5)
    
    plot_width <- 1+2*5
    
    sr_group_colors["D91", "Color"] <- NA
    sr_group_colors["D1", "Color"] <- NA
    sr_group_colors["D29", "Color"] <- NA
    # do above but show p -arm
    for(v1 in day_visno){
      
      sr_group_order <- unique(c(sr_groups[grepl(v1, sr_groups)]))
      
      v_manuf_length <- length(unique(sr_group_data_long %>%  filter(visit_code == v1) %>% pull(v_manuf_code)))
      unique_arms <- length(unique(sr_group_data_long %>%  filter(visit_code == v1) %>% pull(arm_code)))
      nrow_facet <- ceiling(unique_arms*v_manuf_length/7)
      
      plot_width <- 0.5+2*unique_arms
      
      
      if(v1 != last(day_visno)){
        for(v2 in day_visno[c((grep(v1, day_visno)+1):length(day_visno))]) {
          
          
          sr_group_data_temp <- sr_group_data_long %>%
            filter(visit_code %in% c(v1, v2))
          
          v_manuf_length <- length(unique(sr_group_data_temp %>% pull(v_manuf_code)))
          unique_arms <- length(unique(sr_group_data_temp %>% pull(arm_code)))
          nrow_facet <- ceiling(unique_arms*v_manuf_length/7)
          
          plot_width <- 0.5+2*unique_arms
          
          sr_group_order <- unique(c(sr_groups[grepl(v1, sr_groups)], sr_groups[grepl(v2, sr_groups)]))
          
          sr_group_data_temp <- sr_group_data_temp[order(sr_group_data_temp$visit_code),]
          
          plots <- titerlineplot_dodge(sr_group_data_temp, sr_group_colors, titer_thresh = 40, antigens = plot_antigens,
                                 facet_n_row = nrow_facet, sr_group_order = sr_group_order, gmt_facetter = "arm_code", color_by = c("arm_code", "visit_code"),
                                 x_position_by = "age_code", cols_to_keep = c("arm_code", "gmt_arm", "visit_code", "inf_code",
                                                                                "age_code", "v_manuf_code"), show_group_count = TRUE,
                                 show_mean_line = T, mean_line_color = sr_group_colors["P", "Color"], to_long = F,
                                 nrow_gmt = nrow_facet,
                                 dodge_group = "visit_code")$gmt + theme(legend.position = "none")

          ggsave(file.path(figure_dir, paste0(inf_stat, "_", v1, "_", v2, "_gmts_age_all_by_arm_show_p_arm.png")), plot = plots, dpi = 300, width =plot_width, height = 2+2*nrow_facet)
          
          
        }

        # do all three visits combined
        sr_group_data_temp <- sr_group_data_long %>%
            filter(visit_code %in% c("D1", "D29", "D91"))
          
          sr_group_order <- unique(c(sr_groups[grepl("D1", sr_groups)], sr_groups[grepl("D29", sr_groups)], sr_groups[grepl("D91", sr_groups)]))
        
          sr_group_data_temp <- sr_group_data_temp[order(sr_group_data_temp$visit_code),]
          
          sr_group_colors_day <-read.csv("./data/metadata/sr_group_colors_day.csv", sep = ";", row.names = "Serum.group")
    
          sr_group_colors_day_small <- sr_group_colors_day$Color
          names(sr_group_colors_day_small) <- rownames(sr_group_colors_day)
          plots <- titerlineplot_dodge(sr_group_data_temp, sr_group_colors_day, titer_thresh = 40, antigens = plot_antigens,
                                 facet_n_row = nrow_facet, sr_group_order = sr_group_order, gmt_facetter = "arm_code", color_by = c("visit_code"),
                                 x_position_by = "age_code", cols_to_keep = c("arm_code", "gmt_arm", "visit_code", "inf_code",
                                                                                "age_code", "v_manuf_code"), show_group_count = TRUE,
                                 show_mean_line = FALSE, to_long = F,
                                 nrow_gmt = nrow_facet,
                                 dodge_group = "visit_code",
                                 show_lines = FALSE)$gmt +
                                 scale_fill_manual(values = sr_group_colors_day_small) + 
                                 guides(shape="none",
                                  fill = "none") +
                                 scale_color_manual(values = sr_group_colors_day_small,
                                    name = "", breaks = c("D1", "D29", "D91")) + 
                                theme(legend.position = c(.95, .9)) 

          ggsave(file.path(figure_dir, paste0(inf_stat, "_D1_D29_D91_gmts_age_all_by_arm_colour_visit.png")), plot = plots, dpi = 300, width =plot_width, height = 2+2*nrow_facet)
          
        
        
      }
    } 
    
  }
  
}

