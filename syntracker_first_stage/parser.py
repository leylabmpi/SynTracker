import argparse
import os
import re
import config


def parse_arguments():

    error = ""

    parser = argparse.ArgumentParser()

    parser.add_argument("-target", metavar="target_directory_path",
                        help="Path of the target directory which contains metagenome assemblies or genomes", type=str)
    parser.add_argument("-ref", metavar="ref_directory_path",
                        help="Path of the references folder containing the reference genomes", type=str)
    parser.add_argument("-out", metavar="output_directory_path",
                        help="path of the output directory (optional, by default a folder named "
                             + str(config.output_dir) + " will be created under the current directory). "
                             "IMPORTANT: if this directory already exists, it will be written over!!!",
                        type=str, default=config.output_dir)
    parser.add_argument("-mode", metavar="[new|continue]",
                        help="The running mode: 'new' or 'continue' (default='new') "
                             "(Start a new run or continue a previous run that has been stopped).",
                        type=str, default=config.running_mode)
    parser.add_argument("-identity", metavar="blast_identity",
                        help="Minimal blast identity (optional, default=" + str(config.minimal_identity) + ")",
                        type=int, default=config.minimal_identity)
    parser.add_argument("-coverage", metavar="blast_coverage",
                        help="Minimal blast coverage (optional, default=" + str(config.minimal_coverage) + ")",
                        type=int, default=config.minimal_coverage)
    parser.add_argument("-length", metavar="flanking_sequences_length",
                        help="The length of the flanking sequences (from both sides of the BLAST hit). "
                             "(Optional, default=" + str(config.flanking_length) + ")",
                        type=int, default=config.flanking_length)

    # Parse the given arguments
    args = parser.parse_args()

    # Verify that the user provided a target directory
    if args.target is not None:
        input_target_dir = args.target
        if not os.path.isabs(input_target_dir):
            config.input_target_dir = os.path.abspath(input_target_dir) + "/"
        if not re.search(r"^(\S+)\/$", config.input_target_dir):
            config.input_target_dir += "/"
    else:
        error = "Error: you must provide a path to the target folder which contains metagenome assemblies or " \
                "genomes (using -target)\n"
        return error

    # Verify that the user provided a reference directory
    if args.ref is not None:
        input_ref_dir = args.ref
        if not os.path.isabs(input_ref_dir):
            config.input_ref_dir = os.path.abspath(input_ref_dir) + "/"
        if not re.search(r"^(\S+)\/$", config.input_ref_dir):
            config.input_ref_dir += "/"
    else:
        error = "Error: you must provide a path to the references folder containing the reference genomes(using -ref)\n"
        return error

    if args.mode is not None:
        if args.mode == "new" or args.mode == "continue":
            config.running_mode = args.mode
        else:
            error = "-mode ['new' | 'continue'] (Start a new run or continue a previous run that has been stopped).\n"
            return error

    # Set the global variables according to the user's input
    config.output_dir = args.out
    config.minimal_identity = args.identity
    config.minimal_coverage = args.coverage
    config.flanking_length = args.length
    config.jump_length = config.region_length + config.flanking_length * 2

    # Verify that identity is an integer between 0-100
    if config.minimal_identity < 0 or config.minimal_identity > 100:
        error = "Error: the minimal identity for BLAST should be an integer between 0 and 100\n"
        return error

    # Verify that coverage is an integer between 0-100
    if config.minimal_coverage < 0 or config.minimal_coverage > 100:
        error = "Error: the minimal coverage for BLAST should be an integer between 0 and 100\n"
        return error

    return error


def read_conf_file():

    ref_dir = ""
    target_dir = ""
    output_dir = ""
    region_length = ""
    flanking_length = ""
    minimal_coverage = ""
    minimal_identity = ""
    error = ""

    with open(config.conf_file_path) as read_conf:

        in_ref_genomes_list = 0
        in_processed_genomes_list = 0
        for line in read_conf:
            if re.search("^Reference genomes directory:", line):
                m = re.search("^Reference.+\:\s(\S+)\n", line)
                if m:
                    ref_dir = m.group(1)

            elif re.search("^Target", line):
                m = re.search("^Target.+\:\s(\S+)\n", line)
                if m:
                    target_dir = m.group(1)

            elif re.search("^Output", line):
                m = re.search("^Output.+\:\s(\S+)\n", line)
                if m:
                    output_dir = m.group(1)

            elif re.search("^Region", line):
                m = re.search("^Region.+\:\s(\d+)\n", line)
                if m:
                    region_length = m.group(1)

            elif re.search("^Flanking", line):
                m = re.search("^Flanking.+\:\s(\d+)\n", line)
                if m:
                    flanking_length = m.group(1)

            elif re.search("^Minimal coverage", line):
                m = re.search("^Minimal coverage\:\s(\d+)\n", line)
                if m:
                    minimal_coverage = m.group(1)

            elif re.search("^Minimal identity", line):
                m = re.search("^Minimal identity\:\s(\d+)\n", line)
                if m:
                    minimal_identity = m.group(1)

            elif re.search("^Reference genomes:", line):
                in_ref_genomes_list = 1

            elif re.search("^\S+\t\S+\n", line) and in_ref_genomes_list and in_processed_genomes_list == 0:
                m = re.search("^(\S+)\t(\S+)\n", line)
                genome_name = m.group(1)
                genome_file = m.group(2)
                #print("Genome name: " + genome_name)
                config.genomes_dict[genome_name] = dict()
                config.genomes_dict[genome_name]['input_file'] = genome_file
                config.genomes_dict[genome_name]['processed'] = 0

            elif re.search("^Processed reference genomes", line):
                in_processed_genomes_list = 1

            elif re.search("^\S+\n$", line) and in_ref_genomes_list and in_processed_genomes_list:
                m = re.search("^(\S+)\n$", line)
                genome_name = m.group(1)
                #print("Processed Genome name: " + genome_name)
                config.genomes_dict[genome_name]['processed'] = 1

    # Verify that all the parameters were written in the file. If not, print error. If yes, save them in the config
    if ref_dir != "":
        config.input_ref_dir = ref_dir
    else:
        error = "The reference genomes directory is not written in the config file."
        return error

    if target_dir != "":
        config.input_target_dir = target_dir
    else:
        error = "The target genomes directory is not written in the config file."
        return error

    if output_dir != "":
        config.main_output_path = output_dir
    else:
        error = "The output directory is not written in the config file."
        return error

    if region_length != "":
        config.region_length = int(region_length)
    else:
        error = "The region length is not written in the config file."
        return error

    if flanking_length != "":
        config.flanking_length = int(flanking_length)
    else:
        error = "The flanking length is not written in the config file."
        return error

    if minimal_coverage != "":
        config.minimal_coverage = int(minimal_coverage)
    else:
        error = "The minimal coverage is not written in the config file."
        return error

    if minimal_identity != "":
        config.minimal_identity = int(minimal_identity)
    else:
        error = "The minimal identity is not written in the config file."
        return error

    # Verify that there is at least one reference genome and that input files exist
    genomes_counter = 0
    for genome in config.genomes_dict:
        genomes_counter += 1
        if config.genomes_dict[genome]['input_file'] == "":
            error = "No input file for reference genome " + genome
            return error
    # No reference genomes at all
    if genomes_counter == 0:
        error = "No reference genome was read during the previous run"
        return error

    return ""
