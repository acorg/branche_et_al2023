# Setup workspace
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(r3js)
library(ablandscapes)
library(htmlwidgets)
library(webshot2)
library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(patchwork)

source(file.path("functions", "format_data.R"))
source(file.path("functions", "landscape_functions.R"))

set.seed(100)

exclude_2dose <- TRUE

date <- "14DEC2022"
day <- "91"

inf_names <- c(non_inf = "non_inf", inf = "inf")

#set day, visit numbers
day_visno <- c("D1", "D91")

BO_D_ajudstment <- c("wo_2dose")

only_gmt_landscapes <- TRUE

# set optimization nr to fit BA.4/5 position; optimization nr 56 gives upper BA.4/5 for presubmission map; 57 for map with new delta sera
opti_nr <- 1

ag_plot_names <- c("D614G" = "D614G", "B.1.617.2" = "Delta", "B.1.351" = "Beta", "BA.1" = "BA.1", "BA.2.12.1" = "BA.2.12.1", "BA.4/BA.5" = "BA.4/BA.5")


path_to_data <- file.path("data", "metadata")
sr_group_code <- read.csv(file.path(path_to_data, "sr_group_code.csv"), sep = ";")
sr_group_colors <- read.csv("./data/metadata/colour_scheme_emmes.csv", sep = ";", row.names = "Serum.group")

# this map is the one without the new delta sera, October 20202
map <- read.acmap(file.path("data", "maps", "map_ndsubset_no_outliers_slope_adjusted.ace"))
agShown(map)[agNames(map) %in% c("BA.1.1", "BA.1+A484K","BA.3", "BA.2.12.1")] <- FALSE
#align pre submission map
map <- realignMap(map, map)


map <- keepSingleOptimization(map, optimization_number = opti_nr)
lims <- Racmacs:::mapPlotLims(map, sera = FALSE)

#---------- Check for revision comments
srOutline(map)[srGroups(map) %in% c("2x mRNA-1273 new", "3x mRNA-1273", "3x mRNA-1273 BD01", "3x mRNA-1273 BD29", "3x mRNA-1273 (6 month)")] <- "grey30"
srSize(map)[srGroups(map) %in% c("2x mRNA-1273 new", "3x mRNA-1273", "3x mRNA-1273 BD01", "3x mRNA-1273 BD29", "3x mRNA-1273 (6 month)")] <- 5
#-----------------------------
png("figures/base_map_dark.png", 6, 5, units = 'in', res=300, pointsize = 12)
par(mar = rep(0.5, 4))
plot(map, xlim = lims$xlim, ylim = lims$ylim)
dev.off()

# adjust limits for pre resubmission map
lims$xlim[2] <- lims$xlim[2] + 1
lims$ylim[1] <- lims$ylim[1] + 1

if(opti_nr != 1 & length(map$optimizations) > 1){
  map <- keepSingleOptimization(map, optimization_number = opti_nr)
}


# read the data
if(date == "29JUN2022"){
  data <- read.csv(file.path("data","titer_data", date, "COVAIL_Landscape_data_Moderna_D15.csv"), sep= ",")
  plot_antigens <- c("D614G", "B.1.617.2", "B.1.351", "BA.1")
  
  arm_code_order_no_vacc <- c("P", "P+O", "D+O", "O", "B+O", "B+O, B+O")
  
} else if(date %in% c("19SEP2022", "14DEC2022")) {
  
  data <- read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Moderna_D", day,".csv")), sep= ",") %>%
    plyr::rbind.fill(read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Pfizer_D", day,".csv")), sep= ",")) %>%
    plyr::rbind.fill(read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Sanofi_D", day,".csv")), sep= ",")) 
  
  plot_antigens <- c("D614G", "B.1.617.2", "B.1.351", "BA.1", "BA.4/5")

  arm_code_order_no_vacc <- c("P","O","P+O", "D+O","B+O", "B+O, B+O", "P+B", "B")
} else {
  data <- read.csv(file.path("data","titer_data", date, paste0("COVAIL_Landscape_data_Monogram_Moderna_D", day,".csv")), sep= ",")
  plot_antigens <- c("D614G", "B.1.617.2", "B.1.351", "BA.1", "BA.4/5")
  
  arm_code_order_no_vacc <- c("P","O","P+O", "D+O","B+O", "B+O, B+O")
}


# format the data
data <- format_data(data, sr_group_code, combine_b_o_arm = FALSE, remove_all_bfl = FALSE, remove_all_subj_entries = FALSE)


if(exclude_2dose){
  arm_code_order_no_vacc <- arm_code_order_no_vacc[arm_code_order_no_vacc != "B+O, B+O"]
  data <- data %>%
    filter(arm_code != "B+O, B+O:M")
}

data_b <- data

for(d_adj in BO_D_ajudstment){
  
  # now drop B+O, B+O arm day 29 make copy of B+O arm day 
  data <- BO_format_visit_code(data_b, d_adj)
  
# filter out non responders
non_responders <- data %>% filter(visit_code %in% c("D1", "D91")) %>%
  filter(D614G == 20) %>% filter(B.1.617.2 == 20) %>%
  filter(BA.1 == 20) %>% filter(B.1.351 == 20) %>% filter(`BA.4/5` == 20) %>% pull(subjid)

non_responders_table <- data %>%
  filter(subjid %in% non_responders) %>%
  select(subjid, inf_code) %>%
  unique() %>%
  table() 

sum(non_responders_table[,"inf"])
sum(non_responders_table[,"non_inf"])

data <- data %>% filter(!(subjid %in% non_responders))
warning(paste0(length(unique(non_responders)), " subjects were filtered out because they were non responders"))

# now combine v1, v3 but group by arm and infection
sr_group_data <- group_data_sr_group(data, by_arm = T, by_visit = T, by_infection = T, by_age = F) %>%
  arm_code_wo_vacc()

if(only_gmt_landscapes & d_adj %in% c("D57D85")){
  sr_group_data <- sr_group_data %>%
    filter(!(arm_code == "B+O, B+O" & visit_code == "D1"))
}

# Highlight certain antigens
highlighted_ags <- c(plot_antigens, "BA.4/BA.5(2)")
highlighted_ags <- gsub("BA.4/5", "BA.4/BA.5", highlighted_ags)
arm_cols <- sr_group_colors[unique(sr_group_data$arm_code), "Color"]
names(arm_cols) <- unique(sr_group_data$arm_code)

sr_group_data_b <- sr_group_data

for(inf_stat in inf_names){
  
figure_dir <- file.path("figures",date,inf_stat, "landscapes", paste0(d_adj,"_b+o_adj"), paste0("optimization_", opti_nr))
suppressWarnings(dir.create(figure_dir, recursive = T))



sr_group_data <- sr_group_data_b %>% 
  filter(inf_code == inf_stat)

titers_adjusted <- sr_group_data
colnames(titers_adjusted) <- gsub("BA.4/5", "BA.4/BA.5", colnames(titers_adjusted))
titers_adjusted["BA.4/BA.5(2)"] <- titers_adjusted[,"BA.4/BA.5"]

titerdata <- titers_adjusted
titerdata %>%
  pivot_longer(
    cols = highlighted_ags,
    names_to = "variant",
    values_to = "titer"
  ) -> titerdata

titerdata %>%
  group_by(
    arm_code,
    visit_code,
    v_manuf_code
  ) -> titerdata

#titerdata <- titerdata %>% filter(!(variant %in% c("BA.4/BA.5", "BA.4/BA.5(2)")))
titerdata %>%
  group_map(
    get_titertable
  ) -> titertables

lndscp_fits <- lapply(
  titertables,
  function(titertable) {
    
    ablandscape.fit(
      titers = titertable[,c("BA.4/BA.5", "BA.1", "B.1.617.2", "B.1.351", "D614G"), drop = FALSE],
     # titers = titertable[,c("BA.1", "D614G"), drop = FALSE],
      bandwidth = 1,
      degree = 1,
      method = "cone",
      error.sd = 1,
      acmap = map,
      control = list(
        optimise.cone.slope = TRUE
      )
    )
    
  }
)


titertables_groups <- group_data(titerdata)

# Add impulses
titerdata %>%
  group_by(
    visit_code,
    arm_code,
    v_manuf_code,
    variant
  ) %>%
  summarize(gmt =meantiter::mean_titers(titer, method ="bayesian", dilution_stepsize = 0)$mean)-> gmt_data

# angle for html page
angle <- list(
  rotation = c(-1.5, -0.0093, -0.171), #rotation = c(-1.5335, -0.0093, -0.171),
  translation = c(0, 0,0), #translation = c(0.0344, 0.0459, 0.1175),
  zoom = 1.45
)

data3js <- base_plot_data3js(map, lndscp_fits, highlighted_ags, lims, ag_plot_names)

for(v_manuf in unique(titertables_groups$v_manuf_code)){
  
  for(v1 in day_visno){
    
    visit_rows <- c(grep("TRUE", v1 == titertables_groups$visit_code))
        manuf_rows <- grep(v_manuf, titertables_groups$v_manuf_code)
        
        #remove B+O, B+O D1 landscape
        
        lndscp_fits_t <- lndscp_fits[intersect(manuf_rows, visit_rows)]
        titertables_groups_t <- titertables_groups[intersect(manuf_rows, visit_rows),]
        
    if(v1 != last(day_visno)){
      
      for(v2 in day_visno[c((grep("TRUE", v1 == day_visno)+1):length(day_visno))]) {
        # get matching lndscp_fits
        visit_rows <- c(grep("TRUE", v1 == titertables_groups$visit_code), grep("TRUE", v2 == titertables_groups$visit_code))
        manuf_rows <- grep(v_manuf, titertables_groups$v_manuf_code)
        
        #remove B+O, B+O D1 landscape
        
        lndscp_fits_t <- lndscp_fits[intersect(manuf_rows, visit_rows)]
        titertables_groups_t <- titertables_groups[intersect(manuf_rows, visit_rows),]
        
        lndscp_3js <- plot_landscapes_from_list(data3js, titertables_groups_t, lndscp_fits_t, map, gmt_data, highlighted_ags, ag_plot_names,
                                                arm_cols = arm_cols)
        lndscp <-r3js(
          lndscp_3js,
          rotation = angle$rotation,
          zoom = angle$zoom
        )
        
        save_name <- file.path(figure_dir,paste0(v1, "_", v2, "_", v_manuf, "_gmt_landscapes"))
        screenshot_html_landscape_to_png(lndscp, save_name)
       # p <- plot_single_landscape_panel(lndscp, label = "", save_name = save_name, delete_html = FALSE)
       # ggsave(paste0(save_name, ".png"), plot = p, dpi = 300, width = 4, height = 4)
      }
      
    }
    
  }
  
}


if(!only_gmt_landscapes){
  
  # now do landscapes of arms, both visits
  for(v_manuf in unique(titertables_groups$v_manuf_code)){
    
    for(v1 in day_visno){
      
      if(v1 != last(day_visno)){
        
        for(v2 in day_visno[c((grep("TRUE", v1 == day_visno)+1):length(day_visno))]) {
          # get matching lndscp_fits
          visit_rows <- c(grep("TRUE", v1 == titertables_groups$visit_code), grep("TRUE", v2 == titertables_groups$visit_code))
          manuf_rows <- grep(v_manuf, titertables_groups$v_manuf_code)
          
          titertables_groups_temp <- titertables_groups[intersect(manuf_rows, visit_rows),]
          lndscp_fits_temp <- lndscp_fits[intersect(manuf_rows, visit_rows)]
          
          individual_landscape_plots <- list()
          
          for(arm in arm_code_order_no_vacc){
            
            arm_rows <- arm==titertables_groups_temp$arm_code
            
            if(!(TRUE %in% arm_rows)){
              x <- ggplot() + theme_void()
              
              individual_landscape_plots <- c(individual_landscape_plots, list(x))
            } else {
              
              lndscp_fits_t <- lndscp_fits_temp[arm_rows]
              titertables_groups_t <- titertables_groups_temp[arm_rows,]
              
              data3js <- base_plot_data3js(map, lndscp_fits, highlighted_ags, lims, ag_plot_names)
              lndscp_3js <- plot_landscapes_from_list(data3js, titertables_groups_t, lndscp_fits_t, map, gmt_data, highlighted_ags, ag_plot_names,
                                                      arm_cols = arm_cols)
              
              x <-r3js(
                lndscp_3js,
                rotation = angle$rotation,
                zoom = angle$zoom
              )
              individual_landscape_plots <- c(individual_landscape_plots, list(plot_single_landscape_panel(x, label = "", label_size = 4,
                                                                                                           show_border = T)))
            }
            
          }
          
          plotlist <- wrap_plots(individual_landscape_plots, nrow = 1)  
          ggsave(file.path(figure_dir, paste0(v1, "_", v2, "_", v_manuf, "_arm_landscapes.png")), plot = plotlist, dpi = 300, width = 10*length(arm_code_order_no_vacc)/8, height = 2)
        }
        
      }
      
    }
    
  }
}



}
}
