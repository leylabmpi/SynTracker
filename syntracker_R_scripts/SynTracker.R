#requirements
library(DECIPHER)
library(tidyverse)
library(parallel)
library(tools)

###### functions file ######
source("syntracker_R_scripts/SynTracker_functions.R")

# Getting the arguments from the python script
args <- commandArgs(trailingOnly = TRUE)
genome_name <- args[1]
old_new_names_file <- args[2] # Sample names dictionary file
blastdbcmd_output_path <- args[3] # The directory in which the fasta sequences (the output of blastdbcmd) are located
output_folder <- args[4] # The  destination folder for the per-genome output tables
output_summary_folder <- args[5] # The destination folder for the summary output tables (all genomes together)
tmp_folder <- args[6] # A temporary folder for the db files - should be deleted in the end
intermediate_file_folder <- args[7] # If not empty - the folder for saving R intermediate objects (if empty - do not save them)
set_seed_arg <- as.integer(args[8]) # an integer to set the seed (if 0 - do not use seed)
avg_all_regions <- args[9] # If it's True, add an output table without subsampling
core_number <- as.integer(args[10])
metadata_file <- args[11] # If not empty - the path of the metadata file (if empty - there is no metadata)

######################################################################################################

# If the user gave a metadata file, read it
if(metadata_file=='NA') {
    metadata=NA
    print("Running analysis without Metadata")
} else {
    metadata<-read.csv(file=metadata_file, sep=";", header = TRUE)
    print("Running analysis with Metadata file")
}

# List the fasta files from the blastdbcmd output
filepaths <-list.files(path=blastdbcmd_output_path, full.names=TRUE)
region_names<-""
for (i in 1:length(filepaths)) {region_names[i]<-basename(filepaths[i])} # extract the file names from the full path
region_names<-file_path_sans_ext(region_names) # remove file extensions

# Run the Decipher synteny analysis (in multi-core)
print("Running synteny analysis using Decipher...")

# Run synteny analysis, including the calculation of the synteny scores, for each region
# Returns an object containing all the scores for all the regions
synteny_scores_dfs<-mcmapply(synteny_analysis_per_region, filepaths, region_names, tmp_folder, intermediate_file_folder, SIMPLIFY = F, mc.preschedule=F, mc.cores=core_number)

names(synteny_scores_dfs)<-region_names

# Filter out empty data-frames of regions in which synteny has failed or returned an invalid object
synteny_scores_dfs_filtered<-Filter(function(x) nrow(x) > 0, synteny_scores_dfs)

if(intermediate_file_folder != 'NA') {
    saveRDS(synteny_scores_dfs_filtered, file = paste0(intermediate_file_folder,"synteny_scores_tables.rds"))
}

#third part: add names to each table in a new column, merge to one big dataframe, arrange it.
improved_dfs<-map2(synteny_scores_dfs_filtered, names(synteny_scores_dfs_filtered), add_names)
big_dfs<-bind_rows(improved_dfs)  # bind to one dataframe

###################################
# Read the sample names dictionary file and add the original sample names to the table
old_new_names_minimal <- read.table(file=old_new_names_file, header = T, sep="\t")
big_dfs<-left_join(big_dfs, old_new_names_minimal, by=c("sample1"= "new.sample.name")) %>%
dplyr::rename("temp sample1"="sample1", "sample1" := "old.sample.name" )
big_dfs<-left_join(big_dfs, old_new_names_minimal, by=c("sample2"= "new.sample.name")) %>%
dplyr::rename("temp sample2"="sample2", "sample2" := "old.sample.name" )

# change order of sample1 and sample2, according to some rules, so that the order will be uniform throughout the table
# i.e. sampleX-sampleY will always be like that and not sampleY-sampleX ==> if the order is not uniform it will be treated as two different comparisons
big_dfs<-big_dfs %>% mutate(replaced = ifelse(sample1>sample2, "yes", "no")) # add a column specifying if the order of sample 1 and 2 should be replaced (for the sake of Grouping correctly in the next lines)

#change the order of sample specific fields
big_organized_dfs<-big_dfs %>%
mutate(temp=ifelse(replaced == "yes", as.character(sample1), "no need"),  #if replaced == yes: temp will hold sample1
sample1 = ifelse(replaced == "yes", as.character(sample2), as.character(sample1)), #sample1 will hold sample2
sample2=ifelse(replaced == "yes", temp, as.character(sample2)))  #sample2 will hold temp (the original sample1...)

# Create a big summary table containing all the comparisons with the original sample names, sorted by sample1, sample2, region
# (The sorting is important for the reproducibility of the subsampling given the same seed)
big_organized_dfs<-big_organized_dfs %>% arrange(sample1,sample2, ref_genome_region)

# Edit this table for the user to have only the necessary information
big_organized_dfs_final<-big_organized_dfs %>% select(sample1, sample2, ref_genome_region, length1, length2, overlap, blocks, syn_score)

genome_big_table_path<-paste0(output_folder, genome_name , "_synteny_scores_per_region.tab")
write.table(big_organized_dfs_final, file=genome_big_table_path, sep=",", row.names = FALSE)

# Add a column for the ref genome in the final table
summary_organized_dfs<-big_organized_dfs_final %>% mutate(ref_genome=genome_name, .before = sample1)

# Append the per-genome big table to the main table, that contains all the ref-genomes together
summary_table_path<-paste0(output_summary_folder, "synteny_scores_per_region.tab")
write.table(summary_organized_dfs, file=summary_table_path, sep=",", row.names = FALSE, col.names=FALSE, append=TRUE)

##################
# filter and subsample regions
##################
# find the maximal number of regions/pairwise-comparison:
biggest_group<-max(big_organized_dfs_final %>%
                   group_by(sample1, sample2) %>%
                   mutate(regions = n()) %>%
                   ungroup() %>%
                   pull(regions))

# create a list of data frames, with different regions subsampling values (i.e., subsample x regions per pair-wise comparison)
# the conditons is used to avoid subsampling more regions than there are in the biggest group (could result in an error)
regions_sampled<-c(20,30,40,60,80,100,200)

for (i in 1:length(regions_sampled)) {
    if(biggest_group >= regions_sampled[i]){
        grouped_df<-as.data.frame(mapply(subsample_regions, list(big_organized_dfs_final), regions_sampled[i], set_seed_arg, SIMPLIFY = F))

        write.table(grouped_df, file=paste(output_folder, genome_name, "_avg_synteny_scores_", regions_sampled[i], "_regions.txt", sep=""), row.names=FALSE, sep=",")
        grouped_df_with_genome<-grouped_df %>% mutate(ref_genome=genome_name, .before = sample1)
        write.table(grouped_df_with_genome, file=paste(output_summary_folder, "avg_synteny_scores_", regions_sampled[i], "_regions.txt", sep=""),
        row.names=FALSE, col.names=FALSE, append=TRUE, sep=",")
    }
}

# Create an additional output table without subsampling if the user has requested for it
if (avg_all_regions == 'True'){
    grouped_df_avg_all<-big_organized_dfs_final %>%
    group_by(sample1, sample2) %>%
    mutate(regions = n()) %>%
    summarise(average_score=mean(syn_score), compared_regions=mean(regions))

    write.table(grouped_df_avg_all, file=paste(output_folder, genome_name, "_avg_synteny_scores_all_regions.txt", sep=""), row.names=FALSE, sep=",")
    grouped_df_all_with_genome<-grouped_df_avg_all %>% mutate(ref_genome=genome_name, .before = sample1)
    write.table(grouped_df_all_with_genome, file=paste(output_summary_folder, "avg_synteny_scores_all_regions.txt", sep=""),
    row.names=FALSE, col.names=FALSE, append=TRUE, sep=",")
}

unlink(tmp_folder, recursive = T)

cat("\nFinished synteny analysis for:", genome_name, sep = "\n" )

