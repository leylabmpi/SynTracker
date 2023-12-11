#requirements
library(DECIPHER)

############
# fucntions
###########

# synteny_analysis: function to run Decipher analysis for each gene
# inputs: 1. inpath = path to fasta file(s)
#         2. gene name = central region identifier 
#         3. tmp_folder = temporary folder, to store the DECIPHER data. 

# output: synteny object
synteny_analysis<-function(inpath, gene_name, tmp_folder) { 
  #cat("   starting first function - synteny analysis")
  fas<-inpath
  seqs <- readDNAStringSet(fas)
  #print(fas)
  #print(gene_name)
  db<-paste0(tmp_folder,gene_name)
  flag<-0
  for (i in seq_along(seqs)) {
      # The following line adds sequences to the DECIPHER DB: 
      ###########################
      #   IMPORTANT:
      # SPLITTING THE names(seqs[i]) yields the "sample.xxx" from the fasta header of the sequence - therefore only one sequence per sample will be used. 
      # If the str_split is removed, potentially more than one contig per sample can be compared - i.e,, substrains in the same sample. 
      Seqs2DB(seqs[i], "XStringSet", db, str_split(names(seqs[i]), "_")[[1]][1]) 
      ###########################      
      flag<-flag+1
  }
  
  if (flag>1) { # only run DECIPHER analysis if we have more than 1 sequence matching this specific region. 
    synteny <- FindSynteny(db, maxSep=15, maxGap = 15)
  } else {
    synteny<-matrix((1:4),nrow=1)
  }
}


# synteny_scores: a function to do the downstream analysis (i.e., assign per-region pairwise synteny scores) 
# for each region that was analysed with DECIPHER syntenty.
# inputs: 
# 1. synteny_object = (filtered output of "synteny_analysis")
# 2. metadata = metadata file (sea documentation for formatting instructions)

# output: dataframe (columns specified within the function). 

synteny_scores<- function(synteny_object) {
  
  #cat("   starting second function - synteny scores")
  
  # Define a data frame that will hold all the data: change number of Groups according to number of Fields in metadata (i.e., GroupA, GroupB, etc.): here, for example, we have 4
  per_region_table<-data.frame("sample1" = character(),
                             "sample2"= character(),
                             "length1" = integer(),
                             "length2" = integer() , 
                             "overlap" = integer(), #accomulative length of overlapping regions
                             "Blocks" = integer(), # number of synteny blocks per pairwise comparison
                             "synteny_score" = integer(), stringsAsFactors = FALSE)

  # for each two samples, create the values to be kept in the dataframe (assing to one row).
  for (i in 2:ncol(synteny_object)) {
    for (j in 1:(i-1)) {
      sample1<-as.character(str_split(colnames(synteny_object)[j], "_")[[1]][1]) # pay attention to the str_split: could (and should) be changed to suit other naming formats !
      sample2<-as.character(str_split(colnames(synteny_object)[i], "_")[[1]][1]) # pay attention to the str_split: could (and should) be changed to suit other naming formats !
      length1<-synteny_object[j,j][[1]]
      length2<-synteny_object[i,i][[1]]
      overlap<-sum(synteny_object[j,i][[1]][,4])
      blocks<-nrow(synteny_object[i,j][[1]])
      
      # calculate the synteny score
      if (length1>length2) {
        syn_score<-1+log10((overlap/length2)/blocks)
      } else {
        syn_score<-1+log10((overlap/length1)/blocks)
      }
       temprow<-data.frame(sample1,sample2,length1,length2,overlap,blocks,syn_score)
      
       per_region_table<-rbind(per_region_table, temprow)
    }
  }
  return(per_region_table)
}



# add_names: a function to add the ID of the central regions that it's homologs are compared (and other information) to each dataframe outputed by "synteny_scores"
# inputs:
# 1. dfs = list of dataframes returned by "synteny_scores"
# 2.names = namas(dfs)

# output: the same list of dataframes, but with addtional columns

add_names<-function(dfs, names) {
  newnames<-str_extract_all(names, "[[:digit:]]+") # extract the numeric parts of the region name
  nth<-length(newnames[[1]]) #find the number of the last element
  newernames<-sapply(newnames,'[', nth) #extract the last element
  contig<-str_split(names, "_")[[1]][9]
  dfs %>% mutate(contig=contig, position_counter=newernames, region=paste(contig, position_counter, sep="_" )) # add to the dataframe
}


# subsample_regions: a function to sample x regions per pairwise comparison and calculate APSS (Average Pairwise Synteny Scores)
# inputs:
# 1. big_organized_dfs = dataframe that holds per-region pairwise comparisons (multiple regions, multiple pairs)
# 2. subsampling_value = how many regions/pairwise to sample. 

# output: dataframe with APSS values (col name is "average_score"), based on "subsampling_value" regions per pair. Pairs with <"subsampling_value" regions are filtererd out. 
subsample_regions<-function(big_organized_dfs, subsampling_value, set_seed_arg) {
  if(set_seed_arg != 0) {
    set.seed(set_seed_arg)
    cat("Using set.seed for reproducable subsampling. Seed: ", set_seed_arg, sep = " " )
    cat("\n\n")
  }
  newdf<-big_organized_dfs %>% 
    # pay attention to the grouping variable below - should match the groups specified in the "synteny_score" function
    group_by(sample1, sample2) %>% 
    filter(n() > subsampling_value-1) %>%
    sample_n(subsampling_value) %>% #subsample "subsampling_value" regions from each group
    mutate(regions = n()) %>%
    summarise(average_score=mean(syn_score), compared_regions=mean(regions)) %>%
    #mutate(ref_genome=genome_name, .before = sample1)
  return(newdf)
}


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









