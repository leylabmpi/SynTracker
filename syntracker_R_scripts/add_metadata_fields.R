# add metadata: a function to add metadata to the list of comparisons
# input:
# 1. list of tables, each holding pairwise comparisons at different sampling depth
# 2. metadata file.
# output: modified synteny score tables
add_metadata<-function(grouped_list,metadata) {
  #metadata<-as.data.frame(metadata)
  #cat("\nStarting metadata function")
  #print(colnames(metadata))
  for(i in colnames(metadata)) { 
    print(i)
    if(i=="Sample") {
      next
    }
    varname1<-paste0(i,".1")
    varname2<-paste0(i,".2")
    cat(varname1)
    metadata_reduced<-metadata %>% select(Sample, !!i)
    grouped_list<-grouped_list %>% 
      left_join(metadata_reduced, by=c("sample1" = "Sample")) %>% 
      dplyr::rename(!!varname1 := !!i) %>%
      left_join(metadata_reduced, by=c("sample2" = "Sample")) %>% 
      dplyr::rename(!!varname2 := !!i) 
  }
  #cat("\nEnding metadata function")
  return(grouped_list)
}

