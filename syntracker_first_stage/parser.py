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
                        help="The path to the output directory\n. "
                             "When running in 'new' mode (the default), "
                             "this argument is optional. By default a folder named \'" + str(config.output_dir) +
                             "\' will be created under the current directory (if the given path already exists, "
                             "it will be written over).\n"
                             "When running in mode 'continue' or 'continue_all_genomes', it is mandatory to provide "
                             "the path to the output directory of the run that is requested to be continued.",
                        type=str, default=config.output_dir)
    parser.add_argument("-mode", metavar="'new'/'continue'/'continue_all_genomes'",
                        help="The running mode. Start a new run or continue a previous run that has been terminated "
                             "(default='new').\n"
                             "'continue' mode: continue from the last reference genome that was previously processed.\n"
                             "'continue_all_genomes' mode: process all the reference genomes again, without repeating "
                             "the stage in which a blast database is built from the target genomes.",
                        type=str, default=config.running_mode)
    parser.add_argument("-blastDB", metavar="blastDB_directory_path",
                        help="The path to the directory which was previously created by syntracker_makeDB.py "
                             "and contains the uniquely renamed target genomes and the blastDB.\n"
                             "This is an advanced optional argument to be used when the blastDB has already been "
                             "created by syntracker_makeDB.py. When using it, there is no need to provide the "
                             "'-target' argument.",
                        type=str)
    parser.add_argument("-cores", metavar="number_of_cores",
                        help="The number of cores to use for the parallelization of the BLAST-related stages. "
                             "(Optional, default is the number of computer available cores).",
                        type=int)
    parser.add_argument("-length", metavar="region_length",
                        help="The length of the compared region. "
                             "(Optional, default=" + str(config.full_length) + ")",
                        type=int, default=config.full_length)
    parser.add_argument("--identity", metavar="blast_identity",
                        help="Minimal blast identity (optional, default=" + str(config.minimal_identity) + ")",
                        type=int, default=config.minimal_identity)
    parser.add_argument("--coverage", metavar="blast_coverage",
                        help="Minimal blast coverage (optional, default=" + str(config.minimal_coverage) + ")",
                        type=int, default=config.minimal_coverage)
    parser.add_argument("--no_seed", help="Set no seed for the subsampling of n regions per pairwise. This means that "
                        "the average synteny scores may change between SynTracker runs. This is an optional parameter. "
                        "By default, a seed=1 is set to enable reproducibility between different runs.",
                        action='store_true', default=False)

    # Parse the given arguments
    args = parser.parse_args()

    # Set the running mode (new run or continue run)
    if args.mode is not None:
        if args.mode == "new" or args.mode == "continue" or args.mode == "continue_all_genomes":
            config.running_mode = args.mode
        else:
            error = "-mode 'new'/'continue'/'continue_all_genomes'\n" \
                    "(Start a new run or continue a previous run that has been stopped).\n"
            return error

    # Verify the mandatory arguments in 'new' mode
    if config.running_mode == "new":
        if args.blastDB is not None:
            blastDB_dir = args.blastDB

            # Not absolute path -> turn it into absolute
            if not os.path.isabs(blastDB_dir):
                config.syntracker_makeDB_dir = os.path.abspath(blastDB_dir) + "/"
            # Absolute path
            else:
                config.syntracker_makeDB_dir = blastDB_dir
                # Add ending slash
                if not re.search(r"^(\S+)\/$", config.syntracker_makeDB_dir):
                    config.syntracker_makeDB_dir += "/"

            config.is_syntracker_makeDB_dir = True

        else:
            # Verify that the user provided a target directory
            if args.target is not None:
                input_target_dir = args.target

                # Not absolute path -> turn it into absolute
                if not os.path.isabs(input_target_dir):
                    config.input_target_dir = os.path.abspath(input_target_dir) + "/"
                # Absolute path
                else:
                    config.input_target_dir = input_target_dir
                    # Add ending slash
                    if not re.search(r"^(\S+)\/$", input_target_dir):
                        config.input_target_dir += "/"
            else:
                error = "Error: you must provide a path to the target folder which contains metagenome assemblies or " \
                        "genomes (using -target)\n"
                return error

        # Verify that the user provided a reference directory
        if args.ref is not None:
            input_ref_dir = args.ref

            # Not absolute path -> turn it into absolute
            if not os.path.isabs(input_ref_dir):
                config.input_ref_dir = os.path.abspath(input_ref_dir) + "/"
            # Absolute path
            else:
                config.input_ref_dir = input_ref_dir
                # Add ending slash
                if not re.search(r"^(\S+)\/$", input_ref_dir):
                    config.input_ref_dir += "/"
        else:
            error = "Error: you must provide a path to the references folder containing the reference " \
                    "genomes(using -ref)\n"
            return error

        # Set the output directory
        if args.out is not None:
            config.output_dir = args.out

    # Verify that the user provided the previous output directory
    elif config.running_mode == "continue" or config.running_mode == "continue_all_genomes":
        if args.out is not None:
            config.output_dir = args.out
        else:
            error = "Error: in modes 'continue' and 'continue_all_genomes' you must provide a path to the output " \
                    "folder of the run that you wish to continue.\n"
            return error

    # Set the cpu number parameter
    if args.cores is not None and args.cores > 0:
        config.cpu_num = args.cores
    # Get the computer's number of cores
    else:
        config.cpu_num = os.cpu_count()
    print("\nNumber of computer cores to use: " + str(config.cpu_num))

    # Set the global variables according to the user's input
    config.minimal_identity = args.identity
    config.minimal_coverage = args.coverage
    config.full_length = args.length
    config.flanking_length = (config.full_length - config.region_length) / 2
    config.minimal_flanking_length = config.flanking_length * 0.9
    config.minimal_full_length = config.full_length * 0.9

    # Verify that identity is an integer between 0-100
    if config.minimal_identity < 0 or config.minimal_identity > 100:
        error = "Error: the minimal identity for BLAST should be an integer between 0 and 100\n"
        return error

    # Verify that coverage is an integer between 0-100
    if config.minimal_coverage < 0 or config.minimal_coverage > 100:
        error = "Error: the minimal coverage for BLAST should be an integer between 0 and 100\n"
        return error

    if args.no_seed:
        config.is_set_seed = False
        config.seed_num = 0

    return error


def read_conf_file(old_conf_file, new_conf_file, mode):

    current_genome_name = ""

    # In 'continue_all_genomes', the file should be written again without the list of processed genomes
    if mode == "continue_all_genomes":
        out_param = open(new_conf_file, "w")
        out_param.write("Running Parameters:\n")
        out_param.write("--------------------\n")

    with open(old_conf_file) as read_conf:

        in_ref_genomes_list = 0
        in_processed_genomes_list = 0
        for line in read_conf:
            if re.search("^Reference genomes directory:", line):
                m = re.search(r'^Reference.+: (\S+)\n', line)
                if m:
                    config.input_ref_dir = m.group(1)

                if mode == "continue_all_genomes":
                    out_param.write("\n" + line)

            elif re.search("^Target", line):
                m = re.search(r'^Target.+: (\S+)\n', line)
                if m:
                    config.input_target_dir = m.group(1)

                if mode == "continue_all_genomes":
                    out_param.write("\n" + line)

            elif re.search("^syntracker_makeDB", line):
                m = re.search(r'^syntracker_makeDB.+: (\S+)\n', line)
                if m:
                    config.syntracker_makeDB_dir = m.group(1)

                if mode == "continue_all_genomes":
                    out_param.write("\n" + line)

                config.is_syntracker_makeDB_dir = True

            elif re.search("^Output", line):
                m = re.search(r'^Output.+: (\S+)\n', line)
                if m:
                    config.main_output_path = m.group(1)

                if mode == "continue_all_genomes":
                    out_param.write("\n" + line)

            elif re.search("^Full", line):
                m = re.search(r'^Full.+: (\d+)\n', line)
                if m:
                    full_length = m.group(1)

                config.full_length = int(full_length)
                config.flanking_length = (config.full_length - config.region_length) / 2
                config.minimal_flanking_length = config.flanking_length * 0.9
                config.minimal_full_length = config.full_length * 0.9

                if mode == "continue_all_genomes":
                    out_param.write("\n" + line)

            elif re.search("^Minimal coverage", line):
                m = re.search(r'^Minimal coverage: (\d+)\n', line)
                if m:
                    minimal_coverage = m.group(1)
                config.minimal_coverage = int(minimal_coverage)

                if mode == "continue_all_genomes":
                    out_param.write("\n" + line)

            elif re.search("^Minimal identity", line):
                m = re.search(r'^Minimal identity: (\d+)\n', line)
                if m:
                    minimal_identity = m.group(1)
                config.minimal_identity = int(minimal_identity)

                if mode == "continue_all_genomes":
                    out_param.write(line)

            elif re.search("^No seed", line):
                config.is_set_seed = False
                config.seed_num = 0

                if mode == "continue_all_genomes":
                    out_param.write("\n" + line)

                if mode == "continue_all_genomes":
                    out_param.write("\n" + line)

            elif re.search("^Reference genomes:", line):
                in_ref_genomes_list = 1

                if mode == "continue_all_genomes":
                    out_param.write("\n" + line)
                    out_param.write("--------------------\n\n")

            elif re.search(r'^\S+\t\S+\n', line) and in_ref_genomes_list and in_processed_genomes_list == 0:
                m = re.search(r'^(\S+)\t(\S+)\n', line)
                genome_name = m.group(1)
                genome_file = m.group(2)
                config.genomes_dict[genome_name] = dict()
                config.genomes_dict[genome_name]['input_file'] = genome_file
                config.genomes_dict[genome_name]['processed'] = 0
                config.genomes_dict[genome_name]['finished_blast'] = 0
                config.genomes_dict[genome_name]['finished_R'] = 0

                if mode == "continue_all_genomes":
                    out_param.write(line)

            elif re.search("^Processed reference genomes", line):
                # In normal continue mode - read the list of processed genomes
                if mode == "continue":
                    in_processed_genomes_list = 1

                # In 'continue_all_genomes' mode - no need to read the rest of the file
                else:
                    out_param.write("\n" + line)
                    out_param.write("------------------------------\n\n")
                    out_param.close()
                    break

            elif re.search("^ref_genome:", line) and in_processed_genomes_list:
                m = re.search(r'^ref_genome:\s+(\S+)\n$', line)
                current_genome_name = m.group(1)

            elif re.search("BLAST finished", line) and in_processed_genomes_list:
                config.genomes_dict[current_genome_name]['finished_blast'] = 1

            elif re.search("Synteny finished", line) and in_processed_genomes_list:
                config.genomes_dict[genome_name]['finished_R'] = 1
                config.genomes_dict[current_genome_name]['processed'] = 1

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

