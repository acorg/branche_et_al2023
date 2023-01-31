# create mix color for serum groups

blend_colors <- function(colors) {
  blend <- apply(col2rgb(colors), 1, mean)/255
  return(blend)
  
}

blend_hex <- function(colors){
    blended <- blend_colors(colors)
    blended_hex <- rgb(blended[1], blended[2], blended[3])
    return(blended_hex)
}

blend_sr_group_colors <- function(sr_group, sr_group_colors, split_by = NULL) {
  
  # split up sr group
  #target_strata, target_arm, target_age, target_infection, target_variant, target_vaccine, target_nprot
#  sr_group <- gsub(":M|:Pf|:S", "", sr_group)
  sr_group <- gsub(":", "-", sr_group)

  sr_group_split <- str_split(sr_group,pattern = "-")[[1]]

  names(sr_group_split) <- c("STRATA", "arm_code","v_manuf_code", "age_code", "inf_code", "variant_code", "vaccine_code", "nprot_code", "visit_code")
  
  if(length(sr_group_split) > 9){
    sr_group_split["visit_diff"] <- paste0(sr_group_split[10], "-", sr_group_split[11])
  }
  
  # select target identifiers
  if(is.null(split_by)) {
    target_identifiers <- sr_group_split[sr_group_split!= "NA"]
  } else {
    if("sr_group" %in% split_by) {
      target_identifiers <- sr_group_split[sr_group_split!= "NA"]
    } else {
      target_identifiers <- sr_group_split[split_by]
    }
  }
 
  # select colors
  colors <- sr_group_colors[target_identifiers, "Color"]
  colors <- colors[!is.na(colors)]
  
  # blend them
  blended_hex <- blend_hex(colors)
 
  return(blended_hex)
}