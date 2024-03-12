"""
This script takes the target genomes/assemblies, assigns new generic names to the samples (i.e., file name = sample name)
and to the contigs (replaces the original fasta header with contig.xxx).

This is done to avoid:
a. Naming issues, mostly when special characters are included in the file name/fasta headers.
b. Too long fasta headers, which results in abortion of the blast search

All the input target genomes/assemblies are written to one combined file. 
A second file holds a table with old and new Sample names, old and new contig names.
"""

import re
import os
import config


def create_unique_names():

    sample_counter = 0

    target_dir = config.input_target_dir
    output_file = config.combined_renamed_genomes_file_path
    dictionary_table_full = config.dictionary_table_full_path
    sample_only_dictionary = config.sample_dictionary_table_path

    output = open(output_file, "w")
    table = open(dictionary_table_full, "w")
    table.write("old sample name" + "\t" + "new sample name" + "\t" + "old contig header" + "\t" +
                "new contig header" + "\n")
    sample_table = open(sample_only_dictionary, "w")
    sample_table.write("old sample name" + "\t" + "new sample name" + "\n")

    for input_file in os.listdir(target_dir):

        input_path = target_dir + input_file

        if re.search(r"^\.", input_file) or os.path.isfile(input_path) is False:
            continue

        contigs = []
        contig_counter = 0

        with open(input_path) as infile:
            for line in infile:

                # Header line
                if re.search(r"^>", line):
                    contig_counter += 1

                    # check the length of the previous contig
                    if contig_counter > 1:
                        # Add the contig only if it's longer than the 'jump length' (5000 by default)
                        if len(seq) >= config.full_length:
                            temp_dict = dict()
                            temp_dict['header'] = header
                            temp_dict['seq'] = seq + "\n"
                            contigs.append(temp_dict)

                    header = line
                    seq = ""

                # Sequence line
                else:
                    seq += line.rstrip('\n')

        # Check the last contig - add it only if it's longer than the 'jump length' (5000 by default)
        if contig_counter > 0:
            if len(seq) >= config.full_length:
                temp_dict = dict()
                temp_dict['header'] = header
                temp_dict['seq'] = seq + "\n"
                contigs.append(temp_dict)

        # The file contains no contigs
        else:
            print("\nThe file " + input_path + " contains no contig.\n")

        # Verify that there is at least one contig which meets the length criteria in order to include the sample
        if len(contigs) > 0:
            sample_counter += 1
            newfile_name = 'Sample.' + str(sample_counter)
            old_sample_name = os.path.splitext(input_file)[0]

            # Write the minimal dictionary, with the old/new sample names only
            sample_table.write(old_sample_name + "\t" + newfile_name + "\n")

            # A loop over the valid contigs
            for i in range(len(contigs)):
                new_header = ">" + newfile_name + "_contig." + str(i+1) + "\n"

                # Write the contig sequence to the pre-database file
                output.write(new_header)
                output.write(contigs[i]['seq'])

                # Write the contig new and old name to the dictionary table
                table.write(old_sample_name + "\t" + newfile_name + "\t" + contigs[i]['header'] + "\t" + "contig." +
                            str(i+1) + "\n")

            print("Found " + str(contig_counter) + " contigs, from which " + str(len(contigs)) + " are length-valid\n")

        else:
            print("\nThe file " + input_path + " contains no valid contig.\n")

    table.close()
    sample_table.close()
    output.close()
