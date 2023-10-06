"""
This script takes the target genomes/assemblies, assigns new generic names to the samples (i.e., file name = sample name)
and to the contigs (replaces the original fasta header with contig.xxx).

This is done to avoid:
	a. Naming issues, mostly when special characters are included in the file name/fasta headers.
	b. Too long fasta headers, which results in abortion of the blast search  

All the input target genomes/assemblies are written to one combined file. 
A second file holds a table with old and new Sample names, old and new contig names.

* user should provide the path to the folder with the target genomes
* user should provide the paths of the output and table files
* execution: python old_to_new_names.py (target_genomes_folder) (output_file) (table_file)
"""

import re
import os
import config


def create_unique_names():

    file_counter = 0

    target_dir = config.input_target_dir
    output_file = config.combined_renamed_genomes_file_path
    dictionary_table = config.dictionary_table_path

    output = open(output_file, "w")
    table = open(dictionary_table, "w")
    table.write("old sample name" + "\t" + "new sample name" + "\t" + "old contig header" + "\t" +
                "new contig header" + "\n")

    for input_file in os.listdir(target_dir):

        input_path = target_dir + input_file

        if re.search(r"^\.", input_file) or os.path.isfile(input_path) is False:
            continue

        file_counter += 1
        newfile_name = 'Sample.' + str(file_counter)

        old_sample_name = os.path.splitext(input_file)[0]

        contig_counter = 0  # zero the contig counter
        with open(input_path) as infile:
            for line in infile:
                if re.search(r"^>", line):
                    contig_counter += 1
                    line = line.rstrip('\n')
                    # Replace the target string
                    newline = ">" + newfile_name + "_contig." + str(contig_counter) + "\n"
                    output.write(newline)
                    table.write(old_sample_name + "\t" + newfile_name + "\t" + line + "\t" + "contig." +
                                str(contig_counter) + "\n")
                else:
                    output.write(line)

    table.close()
    output.close()
