import os
import shutil
import re
import argparse
import config
import syntracker_first_stage.target_genomes as tr
import syntracker_first_stage.blast as blast


def parse_arguments():

    error = ""
    input_target_dir = ""
    output_dir = "blastDB_output/"  # Default output dir

    parser = argparse.ArgumentParser()

    parser.add_argument("-target", metavar="target_directory_path",
                        help="Path of the target directory which contains metagenome assemblies or genomes", type=str)

    parser.add_argument("-out", metavar="output_directory_path",
                        help="The path to the output directory which will contain the uniquely renamed target genomes "
                             "and the created blastDB.",
                        type=str, default="blastDB")

    # Parse the given arguments
    args = parser.parse_args()

    if args.target is not None:
        input_target_dir = args.target

        # Not absolute path -> turn it into absolute
        if not os.path.isabs(input_target_dir):
            input_target_dir = os.path.abspath(input_target_dir) + "/"
        # Absolute path
        else:
            # Add ending slash
            if not re.search(r"^(\S+)\/$", input_target_dir):
                input_target_dir += "/"

    else:
        error = "Error: you must provide a path to the target folder which contains metagenome assemblies or " \
                "genomes (using -target)\n"
        return error, input_target_dir, output_dir

    # Set the output directory
    if args.out is not None:
        output_dir = args.out

    # Not absolute path -> turn it into absolute
    if not os.path.isabs(output_dir):
        output_dir = os.path.abspath(output_dir) + "/"
    # Absolute path
    else:
        # Add ending slash
        if not re.search(r"^(\S+)\/$", output_dir):
            output_dir += "/"

    return error, input_target_dir, output_dir


def main():

    # Parse the command-line arguments
    error, config.input_target_dir, config.main_output_path = parse_arguments()

    # Exit if there was a problem with the user input
    if error != "":
        print(error)
        exit()

    # If the output dir already exists - delete its content
    if os.path.exists(config.main_output_path):
        print("\nDirectory " + config.main_output_path + " already exists - deleting its content")
        shutil.rmtree(config.main_output_path)
    # Create a new output dir
    try:
        os.makedirs(config.main_output_path)
    except OSError:
        print("\nmkdir " + config.main_output_path + " has failed")
        exit()

    # Set the full path of the logfile
    config.logfile_path = config.main_output_path + config.logfile
    # Open the logfile
    logfile = open(config.logfile_path, "w")

    ################################################################
    # Take care of the naming issues of the target genomes:
    # a. change sample names (i.e., assemblies/genomes) to Sample.xxx
    # b. change fasta headers to contig.xxx
    # c. merge fasta files to one file
    # d. create a table with old and new names.
    print("\nStart merging metagnome-assemblies / genome files\n")
    logfile.write("\nStart merging metagnome-assemblies / genome files\n")

    # Create a folder for the combined targets
    config.combined_output_path = config.main_output_path + config.combined_output_dir
    try:
        os.makedirs(config.combined_output_path)
    except OSError:
        print("\nmkdir " + config.combined_output_path + "has failed")
        exit()
    config.combined_renamed_genomes_file_path = config.combined_output_path + config.combined_renamed_genomes
    config.dictionary_table_full_path = config.combined_output_path + config.dictionary_table_full
    config.sample_dictionary_table_path = config.combined_output_path + config.sample_dictionary_table

    tr.create_unique_names(logfile)

    print("\nMerging metagnome-assemblies / genome files is complete")
    logfile.write("Merging metagnome-assemblies / genome files is complete\n")
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
    blast.make_blast_db(logfile)
    logfile.close()


if __name__ == '__main__':
    main()
