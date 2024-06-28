
# add metadata: a function to add metadata to the list of comparisons
# input:
# 1. list of tables, each holding pairwise comparisons at different sampling depth
# 2. metadata file.
# output: modified synteny score tables
add_metadata<-function(grouped_list, metadata) {
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
  return(grouped_list)
}


###########################################################
#  functions to create synteny-based distance tree  of samples
##########################################################

# dist_mats: a function to turn the final dataframe (with APSS values) to a symmetric distance matrix. 

# inputs:
# 1. dfss = list of dataframe with (at least) sample info and APSS values

# output: symetric distance matrix

dist_mats<-function(dfss) { #input is a list of tibbles

  dfs<-dfss %>% ungroup() %>% #select relevant columns
    select(sample1,sample2,average_score)  
  dfs<-as.data.frame(dfs) #create a second df, with the same values but the samples are in reverse order (to create a symetric matrix)
  dfs2<-data.frame(dfs$sample2,dfs$sample1, dfs$average_score)
  colnames(dfs2)<-c("sample1","sample2", "average_score")
  
  newdfs<-rbind(dfs, dfs2) %>% #bind the two dataframes
    spread(key = sample2, value = average_score) #convert to wide format
  row.names(newdfs)<-newdfs$sample1
  newdfs$sample1<-NULL #remove column without sample names
  for (i in 1:nrow(newdfs)) { #add 1 for each comparison of a sample against itself
    for (j in 1:ncol(newdfs)) {
      if (i==j) {
        newdfs[i,j]<-1
      }
    } 
  }
  return(newdfs)
}

#############
# tree making: afunction to create a phylogenetic trees based on symmetric distance matrix
# input: synteny_matrix: a symmetirc distance matrix. 
# requires a csv file with identifiers, to annotate the tree. 
# output: a phylo tree (either object or a plot)

tree_making<-function(synteny_matrix) {
  #modify the synteny matrix so that synteny score of 1 equals 0, and less similar distances will have bigger values 
  synteny_matrix<-abs(synteny_matrix-1)

  tree_upgma<-upgma(as.matrix(synteny_matrix))
  # Read csv file I created with the merged identifiers
  merged_identifiers<-read.table(file = "merged_ids.csv", header=TRUE, sep=",")
  #visualize trees
  p <- ggtree(tree_upgma) 
  p<-p %<+% merged_identifiers + geom_tiplab(aes(label = GroupA, col=GroupA)) + geom_nodepoint(color="black", alpha=1/2, size=1.5) +theme_tree2() + 
    xlab("synteny distance") + xlim(0,0.7)
  
  #return(tree_upgma) # put a # here if want to return a plot instead!!!!!!
}      
