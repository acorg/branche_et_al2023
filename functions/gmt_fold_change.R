# calculate fold change
calculate_fold_change_from_ag <- function(data, target_ag = "D614G") {
  
  data$fold_change <- unlist(lapply(1:nrow(data), function(x) {
    
    sr_group <- as.character(data[x, "sr_group"]$sr_group)
    
    hom_ag <- target_ag
    
    hom_gmt <- unique(data[data$ag_name == hom_ag & data$sr_group == sr_group, "logtiter"]$logtiter)
   
    hom_gmt <- (2^hom_gmt*10)
    
    this_titer <- (2^data[x, "logtiter"]$logtiter*10)
  
    if(length(hom_gmt) > 1){
      hom_gmt <- hom_gmt[!is.na(hom_gmt)]
    }

    if(length(hom_gmt) < 1 | length(this_titer) <1 ) {
      ""
    } else if(is.na(hom_gmt) | is.na(this_titer)) {
      ""
    } else if(this_titer > hom_gmt) {
      paste0("+", round(this_titer/hom_gmt,1))
    } else if(this_titer == hom_gmt) {
      "1"
    } else {
      paste0("-", round(hom_gmt/this_titer,1))
    }
    
  }))
  
  return(data)
}