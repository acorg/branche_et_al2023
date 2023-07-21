# function to format csv data, to add sr group and sr name
day_visno_code <- read.csv(file.path("data", "metadata", "visno_study_day.csv"), sep =";")
sr_group_code <- read.csv("./data/metadata/sr_group_code.csv", sep = ";")
sr_group_colors <- read.csv("./data/metadata/sr_group_colors_emmes.csv", sep = ";", row.names = "Serum.group")
sr_group_colors <- read.csv("./data/metadata/sr_group_colors.csv", sep = ";", row.names = "Serum.group")

format_data <- function(data, sr_group_code, lab = "Monogram", combine_b_o_arm = T, 
  remove_all_bfl = TRUE, remove_all_subj_entries = TRUE){
  rownames(sr_group_code) <- sr_group_code$full
 
  # create age column
  data$age_code <- "<65"
  data$age_code[grepl("OLDER", toupper(data$stratalab))] <- ">65" #grepl("001|002|009|010|005|006", data$STRATA)] <- ">65"
  
  # add infected or not infected
  if(lab == "Monogram"){
    data$inf_code <- sr_group_code[data$infstat, "code"]
  } else {
    data$inf_code <- "inf"
    data$inf_code[data$nprotein %in% c("NEGATIVE", "N")] <- "non_inf"
    data$inf_code[data$infstat %in% c("N")] <- "non_inf"
    data$inf_code[grepl("without", data$stratalab)] <- "non_inf"
  }
  
  
  # fill code names
  data$arm_code <- sr_group_code[data$trttrue, "code"]
  data$arm_code_no_vacc<- unlist(lapply(sr_group_code$arm_code, function(x) {
    strsplit(x, ":")[[1]][[1]]
    })
  )
  
  if("nprotein" %in% colnames(data)){
    data$nprot_code <- sr_group_code[data$nprotein, "code"]
    data$nprot_code[data$nprotein == ""] <- "ng"
  } else {
    data$nprotein_code <- data$infstat
  }
  
  data$vacc_code <- data$vacgrp
 
  
  if(combine_b_o_arm){
  #  data$arm_code[data$arm_code == "B+O, B+O"] <- "B+O"
    data$arm_code[data$arm_code == "B+O, B+O:M"] <- "B+O:M"
  }
  # have to add it once variant data comes in
  data$variant_code <- "ng"
  
  # ag name
  data$antigen <- sr_group_code[data$TESTMTHD, "code"]
  
  # by visit
  data$visit_code <- unlist(lapply(1:nrow(data), function(x){
    temp_visno <- data[x, "VISNO"]
    
    if(data[x, "trttrue"] == "B+O, B+O:M"){
     
      paste0("D",day_visno_code$DAY_DOUBLE[day_visno_code$VISNO == temp_visno])
    } else {
      paste0("D",day_visno_code$DAY_SINGLE[day_visno_code$VISNO == temp_visno])
    }
    }))
    
  
  # add lab that did triations
  data$lab_code <- lab
  
  # add titer method
  if(lab != "Monogram"){
    data$titer_type <- data$D2TYPE
  } else {
    data$titer_type <- "PNEUT50"
  }
  
  # now add sr name
  data$sr_name <- unlist(lapply(1:nrow(data), function(x){
    paste(data$STRATA[x], data$arm_code[x], data$age_code[x], data$inf_code[x], 
          data$vacc_code[x], data$variant_code[x], data$nprot_code[x], data$visit_code[x],
          data$subjid[x], data$titer_type[x], data$lab_code[x], sep = "-")
  }))
  
  # potentially make numeric in here 
  if("resultA" %in% colnames(data)){
    data <- data %>%
      select(!resultA)
  }
  
  # add antigen columns
  if(lab == "Monogram"){
    data %>%
      select(!TESTMTHD:RESULTA) %>%
      pivot_wider(names_from = antigen, values_from = ISORRES) -> data
  } else {
    data %>%
      select(!ISSTDTC) %>%
      select(!TESTMTHD:ISLLOQ) %>%
      pivot_wider(names_from = antigen, values_from = ISORRES) -> data
  }
 
  # filter out subjects with breakthrough infections
  if(remove_all_bfl){
    breakthrough_infs <- apply(data,1, function(x){
  
    if(TRUE %in% (x[grep("btfl", colnames(data))] == "Y")){
      return(x["subjid"])
    }
    })
    breakthrough_infs <- unique(unlist(breakthrough_infs))

    data <- data %>%
      filter(!(subjid %in% breakthrough_infs))
  } else {
    # create function here to only remove time points before bfl
    data <- remove_data_after_breakthrough(data)
  }
  
  # subjects to remove
  remove_subj <- c("COV.01897", "COV.01552", "COV.01935", "COV.00667")
  remove_d91 <- c("COV.01178")

  if(remove_all_subj_entries){
    data <- data %>%
      filter(!(subjid %in% c(remove_subj, remove_d91)))
  } else {
    data <- data %>%
      filter(!(subjid %in% remove_subj)) %>%
      filter(!(subjid == remove_d91 & visit_code == "D91"))
  }

  
  return(data)
}

# function to filter out data after a breakthrough infection
remove_data_after_breakthrough <- function(data){
  
  # test filtering breakthrough infections
  sub_df <- data[,c(c("subjid","VISNO"), colnames(data)[grep("btfl", colnames(data))])]

  data$btfl_visno <- apply(sub_df,1, function(x){
  
    if(TRUE %in% (x[grep("btfl", colnames(sub_df))] == "Y")){
      min_bf <- colnames(sub_df)[x == "Y"]
      min_bf <- gsub("btfl", "", min_bf)
      min_bf <- min(as.numeric(min_bf), na.rm =T)
      return(min_bf)
    } else {
      return(1000)
    }
    })

  data <- data %>% filter(VISNO < btfl_visno)

  return(data)
}

group_data_sr_group <- function(data, by_strata = F, by_arm = T, by_age = F, by_infection = F, by_variant =F, by_vaccine = F, by_nprotein = F,
                                by_visit = F) {
  
  # add sr group with the pattern:
  # strata-arm-age-infected-variant-vaccine
  data$sr_group <- ""
  
  if(by_strata) {
    data$sr_group <- paste(data$sr_group, data$STRATA, sep ="")
  } else {
    data$sr_group <- paste(data$sr_group, "NA", sep = "")
  }
  
  if(by_arm) {
    data$sr_group <- paste(data$sr_group, data$arm_code, sep = "-")
  } else {
    data$sr_group <- paste(data$sr_group, "NA", sep = "-")
  }
  
  if(by_age) {
    data$sr_group <- paste(data$sr_group, data$age_code, sep = "-")
  } else {
    data$sr_group <- paste(data$sr_group, "NA", sep = "-")
  }
  
  if(by_infection) {
    data$sr_group <- paste(data$sr_group, data$inf_code, sep = "-")
  } else {
    data$sr_group <- paste(data$sr_group, "NA", sep = "-")
  }
  
  if(by_variant){
    data$sr_group <- paste(data$sr_group, data$variant_code, sep ="-")
  } else {
    data$sr_group <- paste(data$sr_group, "NA", sep = "-")
  }
  
  if(by_vaccine){
    data$sr_group <- paste(data$sr_group, data$vacc_code, sep ="-")
  } else {
    data$sr_group <- paste(data$sr_group, "NA", sep = "-")
  }
  
  if(by_nprotein){
    data$sr_group <- paste(data$sr_group, data$nprot_code, sep ="-")
  } else {
    data$sr_group <- paste(data$sr_group, "NA", sep = "-")
  }
  
  if(by_visit){
    data$sr_group <- paste(data$sr_group, data$visit_code, sep ="-")
  } else {
    data$sr_group <- paste(data$sr_group, "NA", sep = "-")
  }
  
  return(data)
}


filter_data_sr_group <- function(data,target_strata = NA, target_arm = NA, target_age = NA, target_infection = NA, 
                                 target_variant =NA, target_vaccine = NA, target_nprot = NA, target_visit = NA,
                                 antigens = c("D614G","B.1.351", "B.1.617.2", "B.1.1.529")) {
  
  # filter target frame
  filtered <- data
  
  if(!is.na(target_strata)) {
    filtered %>% filter(STRATA == target_strata) -> filtered
  } 
  if(!is.na(target_arm)) {
    filtered %>% filter(arm_code == target_arm) -> filtered
  } 
  if(!(is.na(target_age))) {
    filtered %>% filter(age_code == target_age) -> filtered
  } 
  if(!(is.na(target_infection))) {
    filtered %>% filter(inf_code == target_infection) -> filtered
  }
  if(!(is.na(target_variant))){
    filtered %>% filter(variant_code == target_variant) -> filtered
  }
  if(!(is.na(target_vaccine))){
    filtered %>% filter(vacc_code == target_vaccine) -> filtered
  }
  if(!(is.na(target_nprot))){
    filtered %>% filter(nprot_code == target_nprot) -> filtered
  }
  
  if(!(is.na(target_visit))){
    filtered %>% filter(visit_code == target_visit) -> filtered
  }
 
  sr_group <- paste(target_strata, target_arm, target_age, target_infection, target_variant, target_vaccine, target_nprot, target_visit, sep = "-")
  
  filtered$sr_group <- sr_group
  
  return(filtered)
}

sr_group_log_titers_from_data <- function(data, target_strata = NA, target_arm = NA, target_age = NA, 
                                          target_infection = NA, target_variant =NA, target_vaccine = NA, target_nprot = NA,
                                          target_visit = NA,
                                          antigens = c("D614G","B.1.351", "B.1.617.2", "B.1.1.529")){
  
  # filter target frame
  filtered <- filter_data_sr_group(data, target_strata = target_strata, target_arm = target_arm, target_age = target_age,
                                   target_infection = target_infection, target_variant=target_variant, target_vaccine =  target_vaccine,
                                   target_nprot = target_nprot, target_visit = target_visit, antigens = antigens)
  
  filtered_long <- filtered
  
 if(dim(filtered)[1] == 0) {
   return(NULL)
 }
  
  filtered %>% 
    select("sr_name", antigens) %>%
    column_to_rownames("sr_name") -> filtered
  
  filtered <- t(filtered)
  filtered <- apply(filtered, 1:2, function(x) as.numeric(x))
  filtered <- log2(filtered/10)
 
  return(list(log_titers = filtered, full_sr_group = filtered_long))
}


sr_group_to_column <- function(data){
  
    # add sr group with the pattern:
    # strata-arm-age-infected-variant-vaccine
  
  data <- data %>%
    mutate(STRATA = NA,
           arm_code = NA,
           age_code = NA,
           inf_code = NA,
           variant_code = NA,
           vacc_code = NA,
           nprot_code = NA,
           visit_code = NA)
  
  for(x in 1:nrow(data)) {
    
    sr_group_split <- str_split(data$sr_group[x],pattern = "-")[[1]]
    names(sr_group_split) <- c("strata", "arm", "age", "infection", "variant", "vaccine", "nprot", "visit_no")
    
    data[x, "STRATA"] <- sr_group_split["strata"]
    data[x, "arm_code"] <- sr_group_split["arm"]
    data[x, "age_code"] <- sr_group_split["age"]
    data[x, "inf_code"] <- sr_group_split["infection"]
    data[x, "variant_code"] <- sr_group_split["variant"]
    data[x, "vaccine_code"] <- sr_group_split["vaccine"]
    data[x, "nprot_code"] <- sr_group_split["nprot"]
    data[x, "visit_code"] <- sr_group_split["visit_no"]
    
  }
  
  return(data)

}

data_to_long <- function(data, titer_thresh = 40, antigens = c("D614G", "B.1.617.2", "B.1.351", "BA.1")){
  data %>%
    pivot_longer(cols = all_of(antigens), names_to = "ag_name", values_to = "titer") %>%
    mutate(titer = as.numeric(titer),
        logtiter = log2(titer/10),
          logtiter = ifelse(titer < titer_thresh, log2((titer_thresh/20)), logtiter))
}

arm_code_wo_vacc <- function(sr_group_data){
  sr_group_data$v_manuf_code <- sapply(sr_group_data$arm_code, function(x){
    str_split(x, ":")[[1]][2]
  })
  sr_group_data$arm_code <- sapply(sr_group_data$arm_code, function(x){
    str_split(x, ":")[[1]][1]
  })  
  
  return(sr_group_data)
}

BO_format_visit_code <- function(data, to_adjust){
  
  if(to_adjust == "D91"){
    
    # make D85 data to D91
    data <- data %>%
      filter(!(visit_code %in% c("D57"))) %>%
      mutate(visit_code = ifelse(visit_code == "D85", "D91", visit_code))
    
  } else if(to_adjust == "D29") {
    
    # make D85 data to D29
    data <- data %>%
      filter(!(visit_code %in% c("D57"))) %>%
      filter(!(arm_code == "B+O, B+O:M" & visit_code == "D29")) %>%
      mutate(visit_code = ifelse(visit_code == "D85", "D29", visit_code))
    
  } else if(to_adjust == "D57D85"){
    
    # copy B+O, B+O day 85
    d85_data <- data %>%
      filter(arm_code == "B+O, B+O:M" & visit_code == "D85") %>%
      mutate(visit_code = "D29")
    
    d91_data <- d85_data %>%
      mutate(visit_code = "D91")
    
    d57_data <- data %>%
      filter(arm_code == "B+O, B+O:M" & visit_code == "D57") %>%
      mutate(visit_code = "D1")
    
    # now drop B+O, B+O arm day 29 make copy of B+O arm day 
    data <- data %>%
      filter(!(arm_code == "B+O, B+O:M" & visit_code %in% c("D1","D29"))) %>%
      filter(!(visit_code %in% c("D57", "D85"))) 
    
    data <- rbind(data, d85_data, d57_data, d91_data)
    
  } else if(to_adjust == "wo_2dose"){
    
    data <- data %>%
      filter(arm_code != "B+O, B+O:M")
    
  } else{
    warning("B+O, B+O:M visits were not adjusted.")
  }
  
  return(data)
}

# adjust titers to the mean pre titers and adjust post titers by difference in pre landscapes
adjust_titers_to_mean_pre <- function(data, antigens, groupings = c("inf_code", "ag_name"), linear_scale = TRUE){
  
  dots = sapply(groupings, . %>% {as.formula(paste0('~', .))})
  
  data_long <- data_to_long(data, titer_thresh = 40, antigens = antigens)
 
  # pre titers adjusted on log to fit GMT, but linear diff to mean is added to post titers
  # this keeps the linear boost the same
  if(linear_scale){
      data_long %>%
    group_by(.dots = dots) %>%
    mutate(mean_pre = mean(logtiter[visit_code == "D1"])) %>%
    ungroup() %>%
    group_by(arm_code, inf_code, ag_name) %>%
    mutate(mean_arm_pre =  mean(logtiter[visit_code == "D1"]),
            diff_to_pre = mean_pre - mean_arm_pre) %>%
    ungroup() %>%
    mutate(logtiter_pre = logtiter,
           titer_pre = titer,
           logtiter = logtiter_pre + diff_to_pre,
           titer_post = 2^logtiter*10) %>%
    group_by(subjid, ag_name) %>%
    mutate(titer_d1 = titer_pre[visit_code == "D1"],
          titer_d1_post = titer_post[visit_code == "D1"],
          abs_diff = titer_d1_post - titer_d1) %>%
    ungroup() %>%
    mutate(titer = titer_pre + abs_diff,
          logtiter = log2(titer/10)) %>%
    select(!mean_pre:abs_diff) %>%
    select(!logtiter) %>%
    pivot_wider(names_from = ag_name, values_from = titer)-> pre_adjusted
  } else {
    # here adjustment on log scale
    data_long %>%
    group_by(.dots = dots) %>%
    mutate(mean_pre = mean(logtiter[visit_code == "D1"])) %>%
    ungroup() %>%
    group_by(arm_code, inf_code, ag_name) %>%
    mutate(mean_arm_pre = mean(logtiter[visit_code == "D1"]),
            diff_to_pre = mean_pre - mean_arm_pre) %>%
    ungroup() %>%
    mutate(logtiter_pre = logtiter,
           titer_pre = titer,
           logtiter = logtiter_pre + diff_to_pre,
           titer = 2^logtiter*10) %>%
    select(!mean_pre:titer_pre) %>%
    select(!logtiter) %>%
    pivot_wider(names_from = ag_name, values_from = titer)-> pre_adjusted
  }

  return(pre_adjusted)
}


adjust_titers_to_mean_pre_idvl <- function(data, antigens, groupings = c("inf_code", "ag_name"), linear_scale = TRUE){
  
  dots = sapply(groupings, . %>% {as.formula(paste0('~', .))})
  
  data_long <- data_to_long(data, titer_thresh = 40, antigens = antigens)
 
  if(linear_scale){
      data_long %>%
    group_by(.dots = dots) %>%
    mutate(mean_pre = 2^mean(logtiter[visit_code == "D1"])*10) %>%
    ungroup() %>%
    group_by(subjid, ag_name) %>%
    mutate(diff_to_pre = mean_pre - titer[visit_code == "D1"]) %>%
    ungroup() %>%
    mutate(logtiter_pre = logtiter,
           titer_pre = titer,
           titer = titer_pre + diff_to_pre) %>%
    select(!mean_pre:titer_pre) %>%
    select(!logtiter) %>%
    pivot_wider(names_from = ag_name, values_from = titer)-> pre_adjusted
  } else {
    # here adjustment on log scale
    data_long %>%
    group_by(.dots = dots) %>%
    mutate(mean_pre = mean(logtiter[visit_code == "D1"])) %>%
    ungroup() %>%
    group_by(subjid, ag_name) %>%
    mutate(diff_to_pre = mean_pre - logtiter[visit_code == "D1"]) %>%
    ungroup() %>%
    mutate(logtiter_pre = logtiter,
           titer_pre = titer,
           logtiter = logtiter_pre + diff_to_pre,
           titer = 2^logtiter*10) %>%
    select(!mean_pre:titer_pre) %>%
    select(!logtiter) %>%
    pivot_wider(names_from = ag_name, values_from = titer)-> pre_adjusted
  }

  return(pre_adjusted)
}

adjust_titers_by_amount <- function(data, antigens, titer_adjust){
  
  data_long <- data_to_long(data, titer_thresh = 40, antigens = antigens)
 
    # here adjustment on log scale
    data_long %>%
      mutate(logtiter = logtiter + titer_adjust[ag_name],
           titer = 2^logtiter*10) -> data_long
    
    # find subjects with negative titers
    neg_titers <- data_long %>%
      filter(titer < 0) %>%
      pull(subjid) %>%
      unique()

    data_long %>%
    filter(!(subjid %in% neg_titers)) %>%
    select(!logtiter) %>%
    pivot_wider(names_from = ag_name, values_from = titer)-> pre_adjusted
  
  return(list("data" = pre_adjusted, "neg_titer_subjects" = neg_titers))
}


# adjust titers to the mean pre titers and adjust post titers by difference in pre landscapes
adjust_titers_to_mean_pre_neg_titer_inv <- function(data, antigens, groupings = c("inf_code", "ag_name")){
  
  dots = sapply(groupings, . %>% {as.formula(paste0('~', .))})
  
  data_long <- data_to_long(data, titer_thresh = 40, antigens = antigens)
 
  
    data_long %>%
    group_by(.dots = dots) %>%
    mutate(mean_pre = 2^mean(logtiter[visit_code == "D1"])*10) %>%
    ungroup() %>%
    group_by(arm_code, inf_code, ag_name) %>%
      mutate(mean_arm_pre = 2^mean(logtiter[visit_code == "D1"])*10,
            diff_to_pre = mean_pre - mean_arm_pre) %>%
    ungroup() %>%
    mutate(logtiter_pre = logtiter,
           titer_pre = titer,
           titer = titer_pre + diff_to_pre,
           logtiter = ifelse(titer > 0, log2(titer/10), log2(1/abs(titer/10)))) -> pre_adj

    neg_titer_ids <- pre_adj %>%
      filter(titer < 0) %>%
      pull(subjid)
 
    pre_adj %>% filter(subjid %in% neg_titer_ids)-> pre_adjusted
  

  return(pre_adjusted)
}


adjust_titers_to_mean_pre_pos_titer_inv <- function(data, antigens, groupings = c("inf_code", "ag_name")){
  
  dots = sapply(groupings, . %>% {as.formula(paste0('~', .))})
  
  data_long <- data_to_long(data, titer_thresh = 40, antigens = antigens)
 
  
    data_long %>%
    group_by(.dots = dots) %>%
    mutate(mean_pre = 2^mean(logtiter[visit_code == "D1"])*10) %>%
    ungroup() %>%
    group_by(arm_code, inf_code, ag_name) %>%
      mutate(mean_arm_pre = 2^mean(logtiter[visit_code == "D1"])*10,
            diff_to_pre = mean_pre - mean_arm_pre) %>%
    ungroup() %>%
    mutate(logtiter_pre = logtiter,
           titer_pre = titer,
           titer = titer_pre + diff_to_pre,
           logtiter = ifelse(titer > 0, log2(titer/10), log2(1/abs(titer/10)))) -> pre_adj

    neg_titer_ids <- pre_adj %>%
      filter(titer < 0) %>%
      pull(subjid)
 
    pre_adj %>% filter(!(subjid %in% neg_titer_ids))-> pre_adjusted
  

  return(pre_adjusted)
}
