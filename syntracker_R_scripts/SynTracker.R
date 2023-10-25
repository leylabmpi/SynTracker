#requirements
library(DECIPHER)
library(tidyverse)
library(parallel)
library(tools)

###### functions file ######
source("syntracker_R_scripts/SynTracker_functions.R")
source("syntracker_R_scripts/add_metadata_fields.R")

# Getting the arguments from the python script
args <- commandArgs(trailingOnly = TRUE)
genome_name <- args[1]
old_new_names <- args[2] # Sample names dictionary file
blastdbcmd_output_path <- args[3] # The directory in which the fasta sequences (the output of blastdbcmd) are located
output_folder <- args[4] # The destination folder for the output tables
tmp_folder <- args[5] # A temporary folder for the db files - should be deleted in the end
intermediate_file_folder <- args[6] # If not empty - the folder for saving R intermediate objects (if empty - do not save them)
set_seed_arg <- as.integer(args[7]) # an integer to set the seed (if 0 - do not use seed) -> Currently not active - fix it!!!
core_number <- as.integer(args[8])
metadata_file <- args[9] # If not empty - teh path of the metadata file (if empty - there is no metadata)

cat("\ngenome name: ", genome_name, sep = "\n" )
cat("\nold_new_names name: ", old_new_names, sep = "\n" )
cat("\nblastdbcmd_output_path name: ", blastdbcmd_output_path, sep = "\n" )
cat("\noutput_folder name: ", output_folder, sep = "\n" )
cat("\ntmp_folder name: ", tmp_folder, sep = "\n" )
cat("\nintermediate_file_folder: ", intermediate_file_folder, sep = "\n" )
cat("\nset_seed_arg: ", set_seed_arg, sep = "\n" )
cat("\nNumber of cores: ", core_number, sep = "\n" )
cat("\nmetadata_file: ", metadata_file, sep = "\n" )

######################################################################################################
# If the user gave a metadata file, read it
if(metadata_file=='NA') {
    metadata="NA"
    cat("\nRunning analysis without Metadata\n")
} else {
    metadata<-read.csv(file=metadata_file, sep=";", header = TRUE)
    cat("\nRunning analysis with Metadata file\n")
}

# List the fasta files from the blastdbcmd output
filepaths <-list.files(path=blastdbcmd_output_path, full.names=TRUE)
gene_names<-""
for (i in 1:length(filepaths)) {gene_names[i]<-basename(filepaths[i])} # extract the file names from the full path
gene_names<-file_path_sans_ext(gene_names) # remove file extensions

# Run the Decipher synteny analysis (in multi-core)
objects<-mcmapply(synteny_analysis, filepaths, gene_names, tmp_folder, SIMPLIFY = F, mc.preschedule=F, mc.cores=core_number)
names(objects)<-gene_names

# identify iterations of synteny_anlysis that failed for some reason and filter these elements out...
bad_objects_elements <- sapply(objects, inherits, what = "try-error")
objects<-objects[!bad_objects_elements]
narrow<-Filter(function(x) nrow(x) > 1, objects) #filter elements with comparisons of less than 2 valid seqs, if this happened.
#print(length(objects))
rm(objects)

# If the user asked to save intermediate objects - save the narrow ds to the right folder
if(intermediate_file_folder != 'NA') {
#if(!is.na(intermediate_file_folder)){
    saveRDS(narrow, file = paste0(intermediate_file_folder, "narrow.rds"))
}

# second part: Process synteny objects (multi-core processing)
dfs<-mcmapply(synteny_scores, narrow, SIMPLIFY = F, mc.preschedule=F, mc.cores=core_number)
bad_dfs_elements <- sapply(dfs, inherits, what = "try-error") #identify iterations of synteny scores that failed for some reason. Mostly (although very rare), those are two hits for the same region
dfs<-dfs[!bad_dfs_elements] # and filter these elements out...

if(intermediate_file_folder != 'NA') {
#if(!is.na(intermediate_file_folder)){
    saveRDS(dfs, file = paste0(intermediate_file_folder,"dfs.rds"))
}

#third part: add names to each table in a new column, merge to one big dataframe, arrange it.
improved_dfs<-map2(dfs, names(dfs), add_names)
big_dfs<-bind_rows(improved_dfs)  # bind to one dataframe

#####
# add the original sample names to the table (changed in the second step of "find_overlapping_regions.sh")
#####
print("adding new columns now")
old_new_names_minimal<-old_new_names %>% select(new.sample.name, old.sample.name) %>% distinct(new.sample.name, old.sample.name)
big_dfs<-left_join(big_dfs, old_new_names_minimal, by=c("sample1"= "new.sample.name")) %>%
dplyr::rename("temp sample1"="sample1", "sample1" := "old.sample.name" )
big_dfs<-left_join(big_dfs, old_new_names_minimal, by=c("sample2"= "new.sample.name")) %>%
dplyr::rename("temp sampl2"="sample2", "sample2" := "old.sample.name" )

# change order of sample1 and sample2, according to some rules, so that the order will be uniform throughout the table
# i.e. sampleX-sampleY will always be like that and not sampleY-sampleX ==> if the order is not uniform it will be treated as two different comparisons
big_dfs<-big_dfs %>% mutate(replaced = ifelse(sample2>sample1, "yes", "no")) # add a column specifing if the order of sample 1 and 2 should be replaced (for the sake of Grouping correctly in the next lines)

#change the order of sample specific fields
big_organized_dfs<-big_dfs %>%
mutate(temp=ifelse(replaced == "yes", as.character(sample1), "no need"),  #if replaced == yes: temp will hold sample1
sample1 = ifelse(replaced == "yes", as.character(sample2), as.character(sample1)), #sample1 will hold sample2
sample2=ifelse(replaced == "yes", temp, as.character(sample2)))  #sample2 will hold temp (the original sample1...)

# Create a big summary table containing all the comparisons and write it to a file under the main output folder
big_organized_dfs<- big_organized_dfs %>% arrange(sample1,sample2,position_counter) #arrange the order of rows in the table: This is done as in different Operating systems the order is different => leads to different subsampling, even if set.seed() is the same
genome_big_table_path<-paste0(output_folder, genome_name , "_summary.tab")
write.table(big_organized_dfs,file = genome_big_table_path, sep=",", row.names = FALSE)

##################
# filter and subsample regions
##################
# find the maximal number of regions/pairwise-comparison:
biggest_group<-max(big_organized_dfs %>%
                   group_by(sample1, sample2) %>%
                   mutate(regions = n()) %>%
                   ungroup() %>%
                   pull(regions))

# create a list of data frames, with different regions subsampling values (i.e., subsample x regions per pair-wise comparison)
# the conditons is used to avoid subsampling more regions than there are in the biggest group (could result in an error)
grouped_list<-list()
regions_sampled<-c(20,30,40,60,80,100,200)

#run the subsampling function for values below the maximal nuber of regions/pair:
for (i in 1:length(regions_sampled)) {
    ifelse(biggest_group >= regions_sampled[i],
           grouped_list[i]<-mapply(subsample_regions, list(big_organized_dfs), regions_sampled[i], SIMPLIFY = F), #add the set.seed argument to the function
           next)
}

# give names to the elements in the newly filled list
for (i in 1:length(grouped_list)) {names(grouped_list)[i]<-paste0(path_names, "_",regions_sampled[i])}
#print("number of regions found:")
#print(names(grouped_list))
#print("length of grouped list:")
#print(length(grouped_list))
unlink(tmp_folder, recursive = T)

if (exists("metadata")) {
    cat("adding metadata fields\n")
    grouped_list<-mapply(add_metadata,grouped_list, MoreArgs = list(metadata), SIMPLIFY = FALSE)
}

cat("Finished synteny analysis for:", genome_name, sep = "\n" )

mapply(
     write.table,
     x=grouped_list, file=paste(output_folder, names(grouped_list), ".txt", sep=""),
     MoreArgs=list(row.names=FALSE, sep=",")
)
return(grouped_list)

