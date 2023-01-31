
sr_group_gmt_calc <- function(plotdata, thresh, half_thresh = FALSE, cols_to_keep = c(),
  pre_adj = FALSE) {
  plotdata$below_thresh <- plotdata$titer == paste0("<",thresh)
 
  cols_to_keep <- c("sr_group","count", "ag_name", "logtiter", "lower", "upper", "titer", "all_below_thresh", cols_to_keep)
  
  if(pre_adj){
    sr_group_gmt <- sr_group_gmt_pre_adj(plotdata, thresh, half_thresh, cols_to_keep)
    return(sr_group_gmt)
  }
  
  if(TRUE %in% (plotdata$logtiter < 0)){
   
    sr_group_gmt_plotdata <- sr_group_gmt_lin_boost(plotdata, thresh, half_thresh)
  } else {
    
    
    if(!(half_thresh)) {
      # Get gmts by serum group
      plotdata %>%
        group_by(
          sr_group,
          ag_name
        ) %>%
        mutate(
          count = length(titer),
          logtiter = meantiter::mean_titers(
            titers = titer,
            method = "bayesian",
            dilution_stepsize = 0
          )$mean,
          lower = meantiter::mean_titers(
            titers = titer,
            method = "bayesian",
            dilution_stepsize = 0
          )$mean_lower,
          upper = meantiter::mean_titers(
            titers = titer,
            method = "bayesian",
            dilution_stepsize = 0
          )$mean_upper,
          titer = 2^logtiter*10,
          all_below_thresh = !(FALSE %in% below_thresh)
        ) -> sr_group_gmt_plotdata
      
      sr_group_gmt_plotdata$titer <- unlist(lapply(1:nrow(sr_group_gmt_plotdata), function(x) {
        if(sr_group_gmt_plotdata$all_below_thresh[x]) {
          paste0("<",thresh)
        } else {
          sr_group_gmt_plotdata$titer[x]
        }
      }))
      
      sr_group_gmt_plotdata$logtiter <- unlist(lapply(1:nrow(sr_group_gmt_plotdata), function(x) {
        if(sr_group_gmt_plotdata$all_below_thresh[x]) {
          log2(thresh/10/2)
        } else {
          sr_group_gmt_plotdata$logtiter[x]
        }
      }))
    } else {
      
      plotdata$titer[plotdata$below_thresh] <- thresh/2
      
      plotdata %>%
        group_by(
          sr_group,
          ag_name
        ) %>%
        mutate(
          count = length(titer),
          logtiter = meantiter::mean_titers(
            titers = titer,
            method = "bayesian",
            dilution_stepsize = 0
          )$mean,
          lower = meantiter::mean_titers(
            titers = titer,
            method = "bayesian",
            dilution_stepsize = 0
          )$mean_lower,
          upper = meantiter::mean_titers(
            titers = titer,
            method = "bayesian",
            dilution_stepsize = 0
          )$mean_upper,
          titer = 2^logtiter*10,
          all_below_thresh = !(FALSE %in% below_thresh)
        ) -> sr_group_gmt_plotdata
      
    }
    
  }
  
  
  sr_group_gmt_plotdata <- sr_group_gmt_plotdata %>%
    select(all_of(cols_to_keep)) %>%
    unique()
  
  return(sr_group_gmt_plotdata)
}

sr_group_gmt_lin_boost <- function(plotdata, thresh, half_thresh = FALSE) {
  
  if(!(half_thresh)) {
    # Get gmts by serum group
    plotdata %>%
      group_by(
        sr_group,
        ag_name
      ) %>%
      mutate(
        count = length(titer),
        lower = Rmisc::CI(
          logtiter
        )[["lower"]],
        upper = Rmisc::CI(
          logtiter
        )[["upper"]],
        logtiter = Rmisc::CI(
          logtiter
        )[["mean"]],
        titer = 2^logtiter*10,
        all_below_thresh = !(FALSE %in% below_thresh)
      ) -> sr_group_gmt_plotdata
    
    sr_group_gmt_plotdata$titer <- unlist(lapply(1:nrow(sr_group_gmt_plotdata), function(x) {
      if(sr_group_gmt_plotdata$all_below_thresh[x]) {
        paste0("<",thresh)
      } else {
        sr_group_gmt_plotdata$titer[x]
      }
    }))
    
    sr_group_gmt_plotdata$logtiter <- unlist(lapply(1:nrow(sr_group_gmt_plotdata), function(x) {
      if(sr_group_gmt_plotdata$all_below_thresh[x]) {
        log2(thresh/10/2)
      } else {
        sr_group_gmt_plotdata$logtiter[x]
      }
    }))
  } else {
    
    plotdata$titer[plotdata$below_thresh] <- thresh/2
    
    plotdata %>%
      group_by(
        sr_group,
        ag_name
      ) %>%
      mutate(
        count = length(titer),
        logtiter = Rmisc::CI(
          logtiter
        )[["mean"]],
        lower = Rmisc::CI(
          logtiter
        )[["lower"]],
        upper = Rmisc::CI(
          logtiter
        )[["upper"]],
        titer = 2^logtiter*10,
        all_below_thresh = !(FALSE %in% below_thresh)
      ) -> sr_group_gmt_plotdata
    
  }
  
  
  return(sr_group_gmt_plotdata)
}


sr_group_gmt_pre_adj <- function(plotdata, thresh, half_thresh = FALSE, cols_to_keep = c()) {


  # get here arm difference from mean
  plotdata %>%
    group_by(inf_code, ag_name) %>%
    mutate(mean_pre = 2^mean(logtiter[visit_code == "D1"])*10) %>%
    ungroup() %>%
    group_by(arm_code, inf_code, ag_name) %>%
      mutate(mean_arm_pre = 2^mean(logtiter[visit_code == "D1"])*10,
            diff_to_pre = mean_pre - mean_arm_pre) %>%
    ungroup() %>%
    # group_by(ag_name, visit_code, v_manuf_code) %>%
    #   mutate(gmt_arm_log = mean(logtiter[arm_code == "P"]),
    #         gmt_titer = 2^gmt_arm_log*10 + diff_to_pre,
    #         gmt_arm = log2(gmt_titer/10)) %>%
    # ungroup() %>%
    group_by(
        sr_group,
        ag_name
      ) %>%
      mutate(
        count = length(titer),
        lower_lin = 2^Rmisc::CI(
          logtiter
        )[["lower"]]*10 + diff_to_pre,
        upper_lin = 2^Rmisc::CI(
          logtiter
        )[["upper"]]*10 + diff_to_pre,
        logtiter = Rmisc::CI(
          logtiter
        )[["mean"]],
        titer = 2^logtiter*10 + diff_to_pre,
        logtiter = log2(titer/10),
        lower = log2(lower_lin/10),
        upper = log2(upper_lin/10),
        all_below_thresh = !(FALSE %in% below_thresh)
      )  -> sr_group_gmt_plotdata

if("gmt_arm" %in% cols_to_keep){
  sub <- sr_group_gmt_plotdata %>% filter(arm_code == "P:M" | arm_code == "P")
  sr_group_gmt_plotdata$gmt_arm <- sub$logtiter[match(interaction(sr_group_gmt_plotdata$visit_code, sr_group_gmt_plotdata$ag_name), interaction(sub$visit_code, sub$ag_name))]
}
    sr_group_gmt_plotdata <- sr_group_gmt_plotdata %>%
    select(all_of(cols_to_keep)) %>%
    unique()
    
    return(sr_group_gmt_plotdata)
}
