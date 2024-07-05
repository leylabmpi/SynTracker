#requirements
library(DECIPHER)
library(tidyverse)
library(parallel)
library(tools)

###### functions file ######
#source("syntracker_R_scripts/SynTracker_functions.R")

############
# fucntions
###########

# A function that combines two stages and is applied on each region:
# First stage: Synteny analysis using DECIPHER - calculate all vs. all hits
# Second stage: Calculate synteny score for each comparison
synteny_analysis_per_region<-function(inpath, region_name, tmp_folder, intermediate_file_folder, logfile) {

    # Define a data frame that will hold the synteny scores for all the comparisons in the region:
    per_region_table<-data.frame("sample1" = character(),
                               "sample2"= character(),
                               "length1" = integer(),
                               "length2" = integer() ,
                               "overlap" = integer(), #accomulative length of overlapping regions
                               "Blocks" = integer(), # number of synteny blocks per pairwise comparison
                               "synteny_score" = integer(), stringsAsFactors = FALSE)


    cat("\nRunning synteny analysis using Decipher for region ", region_name, sep = "\n" )
    seqs <- readDNAStringSet(inpath)
    db<-paste0(tmp_folder,region_name)
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

    # only run DECIPHER analysis if we have more than 1 sequence matching this specific region.
    if (flag>1) {

        try_catch <- function(exprs) {!inherits(try(eval(exprs)), "try-error")}

        # Synteny run for this region finished successfully and returned a synteny object
        if (try_catch(synteny_object_region<-FindSynteny(db, maxSep=15, maxGap = 15))) {
            cat("\nSynteny analysis for region", region_name, "finished successfully\n", sep = " " )
            cat("\nSynteny analysis for region", region_name, "finished successfully\n", sep = " ", file=logfile, append=TRUE)

        # Synteny analysis for this region returned an error - return an empty matrix
        } else{
            cat("\nSynteny analysis for region", region_name, "could not be completed\n", sep = " ")
            cat("\nSynteny analysis for region", region_name, "could not be completed\n", sep = " ", file=logfile, append=TRUE)

            unlink(db) # Remove the DECIPHER db file
            # Add a '_done' suffix to the region-specific file in the blastdbcmd folder
            # (so this region will not be analysed again in case the user repeats the run in continue mode)
            file.rename(inpath, paste0(inpath,"_done"))
            return(per_region_table)
        }

    # Otherwise, return an emtpy matrix
    } else {
        cat("\nCannot perform synteny analysis for region", region_name, "- only one hit was found\n", sep = " ")
        cat("\nCannot perform synteny analysis for region", region_name, "- only one hit was found\n", sep = " ", file=logfile, append=TRUE)
        # Remove the DECIPHER db file
        unlink(db)
        # Add a '_done' suffix to the region-specific file in the blastdbcmd folder
        # (so this region will not be analysed again in case the user repeats the run in continue mode)
        file.rename(inpath, paste0(inpath,"_done"))
        return(per_region_table)
    }

    # The synteny object was created but something is wrong and it contains no real information
    if (nrow(synteny_object_region[2,1][[1]]) == 0){

        unlink(db) # Remove the DECIPHER db file
        # Add a '_done' suffix to the region-specific file in the blastdbcmd folder
        # (so this region will not be analysed again in case the user repeats the run in continue mode)
        file.rename(inpath, paste0(inpath,"_done"))
        return(per_region_table)
    }

    # Remove the DECIPHER db file
    unlink(db)
    # Add a '_done' suffix to the region-specific file in the blastdbcmd folder
    # (so this region will not be analysed again in case the user repeats the run in continue mode)
    file.rename(inpath, paste0(inpath,"_done"))

    # Continue to calculating the syntey scores only if the synteny analysis returned a valid object
    cat("\nCalculating synteny scores for region ", region_name, sep = "\n" )

    # for each two samples, create the values to be kept in the dataframe (assing to one row).
    for (i in 2:ncol(synteny_object_region)) {
        for (j in 1:(i-1)) {
            sample1<-as.character(str_split(colnames(synteny_object_region)[j], "_")[[1]][1]) # pay attention to the str_split: could (and should) be changed to suit other naming formats !
            sample2<-as.character(str_split(colnames(synteny_object_region)[i], "_")[[1]][1]) # pay attention to the str_split: could (and should) be changed to suit other naming formats !
            length1<-synteny_object_region[j,j][[1]]
            length2<-synteny_object_region[i,i][[1]]
            overlap<-sum(synteny_object_region[j,i][[1]][,4])
            blocks<-nrow(synteny_object_region[i,j][[1]])

            # calculate the synteny score
            if (length1>length2) {
                syn_score<-1+log10((overlap/length2)/blocks)
            }
            else {
                syn_score<-1+log10((overlap/length1)/blocks)
            }
            temprow<-data.frame(sample1,sample2,length1,length2,overlap,blocks,syn_score)

            per_region_table<-rbind(per_region_table, temprow)
        }
    }

    cat("\nSynteny scores calculation for region", region_name, "finished successfully\n", sep = " ")
    cat("\nSynteny scores calculation for region", region_name, "finished successfully\n", sep = " ", file=logfile, append=TRUE)

    # Save the synteny scores object of the current region to be used in a continuing run when needed
    synteny_scores_object_file_name = paste0(region_name, ".rds")
    saveRDS(per_region_table, file = paste0(intermediate_file_folder, synteny_scores_object_file_name))

    return(per_region_table)
}


# add the names (ref_genome_header + position) to the dataframe
add_names<-function(dfs, names) {
  dfs %>% mutate(ref_genome_region=names)
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
    group_by(Sample1, Sample2) %>%
    filter(n() > subsampling_value-1) %>%
    sample_n(subsampling_value) %>% #subsample "subsampling_value" regions from each group
    mutate(regions = n()) %>%
    summarise(Average_score=mean(Synteny_score), Compared_regions=mean(regions)) %>%
    #mutate(ref_genome=genome_name, .before = sample1)
  return(newdf)
}
###########################################################################################################################

# Getting the arguments from the python script
args <- commandArgs(trailingOnly = TRUE)
genome_name <- args[1]
old_new_names_file <- args[2] # Sample names dictionary file
blastdbcmd_output_path <- args[3] # The directory in which the fasta sequences (the output of blastdbcmd) are located
output_folder <- args[4] # The  destination folder for the per-genome output tables
output_summary_folder <- args[5] # The destination folder for the summary output tables (all genomes together)
tmp_folder <- args[6] # A temporary folder for the db files - should be deleted in the end
intermediate_file_folder <- args[7] # In continue mode, it contains the previously processed regions (if any)
set_seed_arg <- as.integer(args[8]) # an integer to set the seed (if 0 - do not use seed)
avg_all_regions <- args[9] # If it's True, add an output table without subsampling
core_number <- as.integer(args[10])
metadata_file <- args[11] # If not empty - the path of the metadata file (if empty - there is no metadata)
logfile <- args[12] # The path of the logfile

######################################################################################################

# If the user gave a metadata file, read it
if(metadata_file=='NA') {
    metadata=NA
    print("Running analysis without Metadata")
} else {
    metadata<-read.csv(file=metadata_file, sep=";", header = TRUE)
    print("Running analysis with Metadata file")
}

# Define an empty list for the synteny scores per region dfs
synteny_scores_dfs_total <- list()
total_region_names <- character()

# Check if there are files in the intermediate_file_folder. If so, it is a continue mode run:
# 1. Read the already processed regions to a list of dfs
# 2. Run the synteny analysis only for the regions that have not been processed yet
# 3. Unite the two lists of dfs for final results
checkfiles <-list.files(path=intermediate_file_folder, full.names=TRUE)
if(length(checkfiles)>0) {
    cat("\nFound", length(checkfiles), "regions that were already processed in a previous run\n", sep = " ", file=logfile, append=TRUE)

    # Create a list of region names
    processed_region_names<-""
    for (i in 1:length(checkfiles)) {processed_region_names[i]<-basename(checkfiles[i])} # extract the file names from the full path
    processed_region_names<-file_path_sans_ext(processed_region_names) # remove file extensions
    total_region_names <- c(processed_region_names)

    # Read the saved RDS score-per-region objects
    for (i in 1:length(checkfiles)) {
        scores_per_region_df <- readRDS(checkfiles[i])
        scores_per_region_df_list <- list(scores_per_region_df)
        synteny_scores_dfs_total <- c(synteny_scores_dfs_total, scores_per_region_df_list)
    }
}

# List the fasta files from the blastdbcmd output
# (files that have already been processed will have '_done' suffix and shouldn't be listed)
filepaths <-list.files(path=blastdbcmd_output_path, full.names=TRUE, pattern="fasta$")
# Make sure that there are files that haven't been processed yet in order to run the DECIPHER analysis
if(length(filepaths)>0) {
    print("\nFound regions that haven't been processed yet\n")
    cat("\nFound", length(filepaths), "regions that haven't been processed yet\n", sep = " ", file=logfile, append=TRUE)
    region_names<-""
    for (i in 1:length(filepaths)) {region_names[i]<-basename(filepaths[i])} # extract the file names from the full path
    region_names<-file_path_sans_ext(region_names) # remove file extensions

    # Run the Decipher synteny analysis (in multi-core)
    print("Running synteny analysis using Decipher...")
    cat("\nGoing to run synteny analysis using Decipher for", length(filepaths), "regions...\n", sep = " ", file=logfile, append=TRUE)

    # Run synteny analysis, including the calculation of the synteny scores, for each region
    # Returns an object containing all the scores for all the regions
    synteny_scores_dfs<-mcmapply(synteny_analysis_per_region, filepaths, region_names, tmp_folder, intermediate_file_folder, logfile, SIMPLIFY = F, mc.preschedule=F, mc.cores=core_number)
    cat("\nNumber of processed regions (not necessarily valid:)", length(synteny_scores_dfs), "\n", sep = " ", file=logfile, append=TRUE)

    synteny_scores_dfs_total <- c(synteny_scores_dfs_total, synteny_scores_dfs)
    total_region_names <- c(total_region_names, region_names)
}

cat("\nTotal number of processed regions in current run and previous runs (not necessarily valid:)", length(synteny_scores_dfs_total), "\n", sep = " ", file=logfile, append=TRUE)

names(synteny_scores_dfs_total)<-total_region_names

# Filter out empty data-frames of regions in which synteny has failed or returned an invalid object
synteny_scores_dfs_filtered<-Filter(function(x) nrow(x) > 0, synteny_scores_dfs_total)

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
big_organized_dfs_final<-big_organized_dfs %>% select(sample1, sample2, ref_genome_region, length1, length2, overlap, blocks, syn_score) %>%
dplyr::rename(Sample1="sample1", Sample2="sample2",  Region="ref_genome_region",  Length1="length1",  Length2="length2",  Overlap="overlap", Blocks="blocks", Synteny_score="syn_score") %>%
dplyr::distinct()

genome_big_table_path<-paste0(output_folder, genome_name , "_synteny_scores_per_region.csv")
write.table(big_organized_dfs_final, file=genome_big_table_path, sep=",", row.names = FALSE)

# Add a column for the ref genome in the final table
summary_organized_dfs<-big_organized_dfs_final %>% mutate(ref_genome=genome_name, .before = Sample1)

# Append the per-genome big table to the main table, that contains all the ref-genomes together
summary_table_path<-paste0(output_summary_folder, "synteny_scores_per_region.csv")
write.table(summary_organized_dfs, file=summary_table_path, sep=",", row.names = FALSE, col.names=FALSE, append=TRUE)

##################
# filter and subsample regions
##################
# find the maximal number of regions/pairwise-comparison:
biggest_group<-max(big_organized_dfs_final %>%
                   group_by(Sample1, Sample2) %>%
                   mutate(regions = n()) %>%
                   ungroup() %>%
                   pull(regions))

# create a list of data frames, with different regions subsampling values (i.e., subsample x regions per pair-wise comparison)
# the conditons is used to avoid subsampling more regions than there are in the biggest group (could result in an error)
regions_sampled<-c(40,60,80,100,200)

for (i in 1:length(regions_sampled)) {
    if(biggest_group >= regions_sampled[i]){
        grouped_df<-as.data.frame(mapply(subsample_regions, list(big_organized_dfs_final), regions_sampled[i], set_seed_arg, SIMPLIFY = F))

        write.table(grouped_df, file=paste(output_folder, genome_name, "_avg_synteny_scores_", regions_sampled[i], "_regions.csv", sep=""), row.names=FALSE, sep=",")
        grouped_df_with_genome<-grouped_df %>% mutate(ref_genome=genome_name, .before = Sample1)
        write.table(grouped_df_with_genome, file=paste(output_summary_folder, "avg_synteny_scores_", regions_sampled[i], "_regions.csv", sep=""),
        row.names=FALSE, col.names=FALSE, append=TRUE, sep=",")
    }
}

# Create an additional output table without subsampling if the user has requested for it
if (avg_all_regions == 'True'){
    grouped_df_avg_all<-big_organized_dfs_final %>%
    group_by(Sample1, Sample2) %>%
    mutate(regions = n()) %>%
    summarise(Average_score=mean(Synteny_score), Compared_regions=mean(regions))

    write.table(grouped_df_avg_all, file=paste(output_folder, genome_name, "_avg_synteny_scores_all_regions.csv", sep=""), row.names=FALSE, sep=",")
    grouped_df_all_with_genome<-grouped_df_avg_all %>% mutate(ref_genome=genome_name, .before = Sample1)
    write.table(grouped_df_all_with_genome, file=paste(output_summary_folder, "avg_synteny_scores_all_regions.csv", sep=""),
    row.names=FALSE, col.names=FALSE, append=TRUE, sep=",")
}

unlink(tmp_folder, recursive = T)

cat("\nFinished synteny analysis for:", genome_name, sep = "\n" )
