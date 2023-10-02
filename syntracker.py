import os
import subprocess
import shutil
import re
import csv
import time
import config
import syntracker_first_stage.parser as parser
import syntracker_first_stage.central_regions as cr
import syntracker_first_stage.target_genomes as tr
import syntracker_first_stage.blast as blast

before = time.time()

config.working_dir = os.getcwd()

# Parse the command-line arguments
error = parser.parse_arguments()

# Exit if there was a problem with the user input
if error != "":
    print(error)
    exit()

# Add a slash to the output dir name
if not re.search("^(\S+)\/$", config.output_dir):
    config.output_dir += "/"

# Set the absolute path of the main output dir
if os.path.isabs(config.output_dir):
    config.main_output_path = config.output_dir
else:
    config.main_output_path = config.working_dir + "/" + config.output_dir

# Set the full path of the running params file
config.conf_file_path = config.main_output_path + config.conf_file

###############################################################################
# The current run is new -> start everything from scratch
if config.running_mode == "new":

    print("\nStarting a new SynTracker run\n")

    # If the output dir already exists - delete its content
    if os.path.exists(config.main_output_path):
        print("\nDirectory " + config.main_output_path + " already exists - deleting its content")
        shutil.rmtree(config.main_output_path)
    # Create a new output dir
    try:
        os.makedirs(config.main_output_path)
    except OSError:
        print("\nmkdir " + config.main_output_path + "has failed")
        exit()

    # Create the running parameters file
    out_param = open(config.conf_file_path, "w")
    out_param.write("Running Parameters:\n")
    out_param.write("--------------------\n")
    out_param.write("\nReference genomes directory: " + config.input_ref_dir + "\n")
    out_param.write("Target genomes directory: " + config.input_target_dir + "\n")
    out_param.write("Output directory: " + config.main_output_path + "\n")
    out_param.write("\nRegion length: " + str(config.region_length) + "\n")
    out_param.write("Flanking regions length: " + str(config.flanking_length) + "\n")
    out_param.write("\nMinimal coverage: " + str(config.minimal_coverage) + "\n")
    out_param.write("Minimal identity: " + str(config.minimal_identity) + "\n")

    ################################################################
    # Take care of the naming issues of the target genomes:
    # a. change sample names (i.e., assemblies/genomes) to Sample.xxx
    # b. change fasta headers to contig.xxx
    # c. merge fasta files to one file
    # d. create a table with old and new names.
    print("\nStart merging metagnome-assemblies / genome files\n")

    # Create a folder for the combined targets
    config.combined_output_path = config.main_output_path + config.combined_output_dir
    try:
        os.makedirs(config.combined_output_path)
    except OSError:
        print("\nmkdir " + config.combined_output_path + "has failed")
        exit()
    config.combined_renamed_genomes_file_path = config.combined_output_path + config.combined_renamed_genomes
    config.dictionary_table_path = config.combined_output_path + config.dictionary_table

    tr.create_unique_names()

    print("\nMerging metagnome-assemblies / genome files is complete")
    out_param.write("\nMerging metagnome-assemblies stage is complete\n")
    config.complete_target_merge = True

    #########################################################
    # Make a blast database from all the contigs.
    config.blast_db_path = config.main_output_path + config.blast_db_dir
    try:
        os.makedirs(config.blast_db_path)
    except OSError:
        print("\nmkdir " + config.blast_db_path + "has failed")
        exit()

    config.blast_db_file_path = config.blast_db_path + config.blast_db_file
    blast.make_blast_db()

    ##############################################################################
    # Extract the reference genomes from the user-defined directory
    # and save the names and file paths in a dictionary and in the conf file
    print("\nStarting reference genomes loop...\n")
    out_param.write("\nReference genomes:\n")
    out_param.write("--------------------\n")

    for ref_genome_file in os.listdir(config.input_ref_dir):
        if re.search(r"^\.", ref_genome_file):
            continue

        # Extract the reference genome's file name without suffix and save it in a dictionary
        ref_genome = os.path.splitext(ref_genome_file)[0]
        ref_genome_file_path = config.input_ref_dir + ref_genome_file
        # Add the genome name and file path to the main dict as new genmome that hasn't been processed yet (=0)
        config.genomes_dict[ref_genome] = dict()
        config.genomes_dict[ref_genome]['input_file'] = ref_genome_file_path
        config.genomes_dict[ref_genome]['processed'] = 0
        # Write the entries to the config file
        out_param.write(ref_genome + "\t" + ref_genome_file_path + "\n")

        # Add the genome name to the list of genomes that should be processed in the current run
        config.run_genomes_list.append(ref_genome)

        out_param.write("\nProcessed reference genomes (BLAST + blastdbcmd):\n\n")

    out_param.close()

#######################################################################################################
# The current run continues a previous run that was stopped -> read the parameters from the conf file
else:
    print("\nContinuing a previous SynTracker run\n")

    # Read the parameters and genome list from the conf file
    error = parser.read_conf_file()

    # Exit if some important parameters were missing from the config file
    if error != "":
        print(error)
        print("Cannot continue the previous run - please rerun your dataset from the beginning (-mode 'new')")
        exit()

    # Assign all the global paths
    config.combined_output_path = config.main_output_path + config.combined_output_dir
    config.combined_renamed_genomes_file_path = config.combined_output_path + config.combined_renamed_genomes
    config.dictionary_table_path = config.combined_output_path + config.dictionary_table
    config.blast_db_path = config.main_output_path + config.blast_db_dir
    config.blast_db_file_path = config.blast_db_path + config.blast_db_file

    # Add the genomes that were not processed to the list of genomes that should be processed in the current run
    for genome in config.genomes_dict:
        if config.genomes_dict[genome]['processed'] == 0:
            config.run_genomes_list.append(genome)
            print("Found reference genome that has not been processed yet: " + genome)

#######################################################################################################################
# From now on: the same operations for both a new run and a continuing run

# Open the config file in append mode
out_param = open(config.conf_file_path, "a")

# A loop over the ref-genomes that should be processed in the current run (all or part - depending on the running mode)
for ref_genome in config.run_genomes_list:

    print("\nProcessing reference genome " + ref_genome)

    # Create a directory for the reference genome under the main output dir
    ref_genome_output_dir = config.main_output_path + ref_genome + "/"

    # If the genome dir already exists - delete its content
    if os.path.exists(ref_genome_output_dir):
        print("\nDirectory " + ref_genome_output_dir + " already exists - deleting its content")
        shutil.rmtree(ref_genome_output_dir)
    # Create a new ref genome dir
    else:
        try:
            os.makedirs(ref_genome_output_dir)
        except OSError:
            print("\nmkdir " + ref_genome_output_dir + "has failed")
            exit()

    #####################################
    # step 1: find the "central_regions"

    # Create a folder for the central regions
    genome_central_regions_dir = ref_genome_output_dir + config.central_regions_dir
    try:
        os.makedirs(genome_central_regions_dir)
    except OSError:
        print("\nmkdir " + genome_central_regions_dir + "has failed")
        exit()

    cr.find_central_regions(ref_genome, genome_central_regions_dir)
    print("\nFound central regions. They are located in: " + genome_central_regions_dir)

    ###########################################
    # Step 2: run blast per central region and extract the flanking sequences for each hit
    print("\nRunning BLAST search for each region of the central regions...")

    # Create a main blast output folder
    genome_blast_out_dir = ref_genome_output_dir + config.blast_out_dir
    try:
        os.makedirs(genome_blast_out_dir)
    except OSError:
        print("\nmkdir " + genome_blast_out_dir + "has failed")
        exit()

    # Create a blastdbcmd output folder (for the hits sequences includeing the flanking regions)
    genome_blastdbcmd_out_dir = ref_genome_output_dir + config.blastdbcmd_out_dir
    try:
        os.makedirs(genome_blastdbcmd_out_dir)
    except OSError:
        print("\nmkdir " + genome_blastdbcmd_out_dir + "has failed")
        exit()

    # List the central regions files for the current ref-genome
    for region_file in os.listdir(genome_central_regions_dir):
        if re.search(r"^.+\.fasta", region_file):  # List only fasta files

            region_name = os.path.splitext(region_file)[0]
            full_path_region_file = genome_central_regions_dir + region_file
            blast_region_outfile = genome_blast_out_dir + region_name + ".tab"
            blastdbcmd_region_outfile = genome_blastdbcmd_out_dir + region_name + ".fasta"

            # Run blast for each region
            blast.run_blastn(full_path_region_file, blast_region_outfile)

            # Count the number of hits in the blast output file
            file = open(blast_region_outfile)
            hits_num = len(file.readlines())  # Count the number of lines (=num of hits) in the output file

            # If the blast output file has less than 2 hits - delete it
            if hits_num < 2:
                try:
                    os.remove(blast_region_outfile)
                except OSError:
                    print("\nRemoving " + blast_region_outfile + "has failed")

            # At least two hits -> read the file
            else:
                with open(blast_region_outfile) as file:
                    reader = csv.reader(file, delimiter='\t')
                    for row in reader:
                        sample_name = row[0]
                        start = row[1]
                        end = row[2]
                        strand = row[3]

                        # Plus strand
                        if strand == "plus":
                            flank_start = int(start) - int(config.flanking_length)
                            flank_end = int(end) + int(config.flanking_length)
                        # Minus strand
                        else:
                            flank_start = int(end) - int(config.flanking_length)
                            flank_end = int(start) + int(config.flanking_length)

                        # If the flanking-start is a negative value
                        if flank_start < 0:
                            # If the upstream length is at least 1800, take this sequence
                            # and change the start position to 1
                            if abs(flank_start) <= config.flanking_length - config.minimal_flanking_length:
                                flank_start = 1
                            # The overall upstream length is shorter than 1800 -> ignore the sequence
                            else:
                                continue

                        # Run blastdbcmd to get the hit including flanking sequences from the database
                        blast.run_blastdbcmd(sample_name, str(flank_start), str(flank_end), strand,
                                             blastdbcmd_region_outfile)

    print("BLAST search and blastdbcmd was completed successfully\n")
    out_param.write(ref_genome + "\n")

out_param.close()

after = time.time()
duration = (after - before)
#print("SynTracker first stage took "+str(duration)+" seconds")





