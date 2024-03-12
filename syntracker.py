import os
import subprocess
import multiprocessing
import shutil
import re
import time
import config
import syntracker_first_stage.parser as parser
import syntracker_first_stage.central_regions as cr
import syntracker_first_stage.target_genomes as tr
import syntracker_first_stage.blast as blast


def main():

    # Get the working directory
    config.working_dir = os.getcwd()

    # Parse the command-line arguments
    error = parser.parse_arguments()

    # Exit if there was a problem with the user input
    if error != "":
        print(error)
        exit()

    # Set the absolute path of the main output directory
    if not os.path.isabs(config.output_dir):
        config.main_output_path = os.path.abspath(config.output_dir) + "/"
    # The user gave an absolute path
    else:
        config.main_output_path = config.output_dir
        # Add ending slash
        if not re.search(r"^(\S+)\/$", config.main_output_path):
            config.main_output_path += "/"

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
            print("\nmkdir " + config.main_output_path + " has failed")
            exit()

        # Create the running parameters file
        out_param = open(config.conf_file_path, "w")
        out_param.write("Running Parameters:\n")
        out_param.write("--------------------\n")
        out_param.write("\nReference genomes directory: " + config.input_ref_dir + "\n")
        out_param.write("\nTarget genomes directory: " + config.input_target_dir + "\n")
        out_param.write("\nOutput directory: " + config.main_output_path + "\n")
        if config.is_metadata:
            out_param.write("\nMetadata file path: " + config.metadata_file_path + "\n")
        out_param.write("\nFull regions length: " + str(config.full_length) + "\n")
        out_param.write("\nMinimal coverage: " + str(config.minimal_coverage) + "\n")
        out_param.write("Minimal identity: " + str(config.minimal_identity) + "\n")
        if config.save_intermediate:
            out_param.write("\nSave intermediate: " + str(config.save_intermediate) + "\n")
        if config.is_set_seed is False:
            out_param.write("No seed\n")
        if config.avg_all:
            out_param.write("Average all regions\n")

        ############################################
        # Create a directory for the summary output (all genomes together)
        config.summary_output_path = config.main_output_path + config.summary_output_dir
        try:
            os.makedirs(config.summary_output_path)
        except OSError:
            print("\nmkdir " + config.summary_output_path + " has failed")
            exit()

        # Prepare R main summary output file
        config.output_summary_file_path = config.summary_output_path + config.output_summary_file
        out_summary = open(config.output_summary_file_path, "w")
        out_summary.write("\"Ref_genome\",\"Sample1\",\"Sample2\",\"Region\",\"Length1\",\"Length2\",\"Overlap\",\"Blocks\",\"Synteny_score\"\n")
        out_summary.close()

        # Prepare R subsampled-regions output files for all the genomes together
        for file_name in config.subsampled_regions_file_names:
            file_path = config.summary_output_path + file_name
            out_subsampled = open(file_path, "w")
            out_subsampled.write("\"Ref_genome\",\"Sample1\",\"Sample2\",\"Average_score\",\"Compared_regions\"\n")
            out_subsampled.close()

        # If the user added the --avg_all option, create also a file for the average-all-regions output
        if config.avg_all:
            file_path = config.summary_output_path + config.avg_all_file_name
            out_avg_all = open(file_path, "w")
            out_avg_all.write("\"Ref_genome\",\"Sample1\",\"Sample2\",\"Average_score\",\"Compared_regions\"\n")
            out_avg_all.close()

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
        config.dictionary_table_full_path = config.combined_output_path + config.dictionary_table_full
        config.sample_dictionary_table_path = config.combined_output_path + config.sample_dictionary_table

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
        out_param.write("\nReference genomes:\n")
        out_param.write("--------------------\n\n")

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
            config.genomes_dict[ref_genome]['finished_blast'] = 0
            config.genomes_dict[ref_genome]['finished_R'] = 0
            # Write the entries to the config file
            out_param.write(ref_genome + "\t" + ref_genome_file_path + "\n")

            # Add the genome name to the list of genomes that should be processed in the current run
            config.run_genomes_list.append(ref_genome)

        out_param.write("\nProcessed reference genomes:\n")
        out_param.write("------------------------------\n\n")

        out_param.close()

    #######################################################################################################
    # The current run continues a previous run that was stopped -> read the parameters from the conf file
    else:
        print("\nContinue a previous SynTracker run\n")

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
        config.dictionary_table_full_path = config.combined_output_path + config.dictionary_table_full
        config.sample_dictionary_table_path = config.combined_output_path + config.sample_dictionary_table
        config.blast_db_path = config.main_output_path + config.blast_db_dir
        config.blast_db_file_path = config.blast_db_path + config.blast_db_file
        config.summary_output_path = config.main_output_path + config.summary_output_dir

        # Add the genomes that were not processed to the list of genomes that should be processed in the current run
        counter = 0
        for genome in config.genomes_dict:
            if config.genomes_dict[genome]['processed'] == 0:
                config.run_genomes_list.append(genome)
                print("Found reference genome that has not been processed yet: " + genome)
                counter += 1

        # Found no genome to process (=all were already processed)
        if counter == 0:
            print("All reference genomes have already been processed!")
            exit(0)

    #######################################################################################################################
    # From now on: the same operations for both a new run and a continuing run

    # A loop over the ref-genomes that should be processed in the current run
    # (all or part - depending on the running mode)
    for ref_genome in config.run_genomes_list:

        before = time.time()

        # Set the directory for the reference genome under the main output dir
        ref_genome_output_dir = config.main_output_path + ref_genome + "/"

        # Set the temp folder
        genome_tmp_out_dir = ref_genome_output_dir + "tmp/"

        # Set the blastdbcmd output folder (for the hits sequences includeing the flanking regions)
        genome_blastdbcmd_out_dir = ref_genome_output_dir + config.blastdbcmd_out_dir

        # Run the directory creation and the BLAST parts only if it is a new run
        # or the BLAST part hasn't finished in a previous run (in case of 'continue' mode)
        if config.running_mode == "new" or config.genomes_dict[ref_genome]['finished_blast'] == 0:

            print("\nProcessing reference genome " + ref_genome)
            out_param = open(config.conf_file_path, "a")
            out_param.write("ref_genome: " + ref_genome + "\n")
            out_param.close()

            # If the genome dir already exists - delete it and its content
            if os.path.exists(ref_genome_output_dir):
                print("\nDirectory " + ref_genome_output_dir + " already exists - deleting its content")
                shutil.rmtree(ref_genome_output_dir)
            # Create a new ref genome dir
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
            try:
                os.makedirs(genome_blastdbcmd_out_dir)
            except OSError:
                print("\nmkdir " + genome_blastdbcmd_out_dir + "has failed")
                exit()

            # Create a tmp output folder
            try:
                os.makedirs(genome_tmp_out_dir)
            except OSError:
                print("\nmkdir " + genome_tmp_out_dir + "has failed")
                exit()

            # Create a list of the central regions files for the current ref-genome
            region_files_list = []
            for region_file in os.listdir(genome_central_regions_dir):
                if re.search(r"^.+\.fasta", region_file):  # List only fasta files
                    region_files_list.append(region_file)

            # A loop over batches of processes in the size of the available number of threads
            batch_counter = 0
            failed = 0
            for batch_index in range(0, len(region_files_list), config.cpu_num):

                batch_processes = []
                batch_counter += 1

                # A loop over the regions of a certain threads batch
                for region_index in range(batch_index, batch_index+config.cpu_num):

                    if region_index < len(region_files_list):

                        region_file = region_files_list[region_index]
                        region_name = os.path.splitext(region_file)[0]
                        full_path_region_file = genome_central_regions_dir + region_file
                        blast_region_outfile = genome_blast_out_dir + region_name + ".tab"
                        blastdbcmd_region_outfile = genome_blastdbcmd_out_dir + region_name + ".fasta"
                        blastdbcmd_region_outfile_tmp = genome_tmp_out_dir + region_name + ".fasta"

                        process = multiprocessing.Process(target=blast.blast_per_region_process,
                                                          args=(full_path_region_file, blast_region_outfile,
                                                                blastdbcmd_region_outfile,
                                                                blastdbcmd_region_outfile_tmp,
                                                                config.blast_db_file_path,
                                                                config.flanking_length, config.minimal_flanking_length,
                                                                config.minimal_full_length,
                                                                config.minimal_identity, config.minimal_coverage,
                                                                config.blast_num_threads))
                        # Start the process
                        process.start()

                        # Add the process to the list for later control
                        batch_processes.append(process)

                # wait until all the processes in the batch are finished
                success_counter = 0
                for proc in batch_processes:
                    proc.join()
                    exit_code = proc.exitcode
                    if exit_code == 0:
                        success_counter += 1

                if success_counter == len(batch_processes):
                    print("\nAll processes in batch number " + str(batch_counter) + " finished successfully")
                elif success_counter >= len(batch_processes) * 0.9:
                    print("\n" +str(success_counter) + " processes in batch number " + str(batch_counter)
                          + " finished successfully")
                else:
                    failed += 1
                    print("\nBatch number " + str(batch_counter) + " has failed or finished without any valid hits")
                    continue

            # All the batches failed or finished without results
            if failed == batch_counter:
                print("\nThe BLAST search has not finish successfully or found not enough valid hits - "
                      "skipping the current reference genome (" + ref_genome + ")...\n")
                continue

            print("\nBLAST search for all the regions finished successfully\n")
            config.genomes_dict[ref_genome]['finished_blast'] = 1
            out_param = open(config.conf_file_path, "a")
            out_param.write("- BLAST finished\n")
            out_param.close()

            after = time.time()
            duration = after - before
            print("\nThe BLAST stage took " + str(duration) + " seconds.\n")

            # Delete the blastdbcmd temporary output folder
            if os.path.exists(genome_tmp_out_dir):
                print("\nRemoving the temporary folder " + genome_tmp_out_dir + "\n")
                shutil.rmtree(genome_tmp_out_dir)

        # In 'continue' mode, where only the BLAST stage of the current ref-genome was finished successfully
        else:
            # Set the names of the folders that should hold the R per-genome output and intermediate files
            final_output_path = ref_genome_output_dir + config.final_output_dir
            r_temp_path = ref_genome_output_dir + config.r_temp_dir

            # It can happen that the R process continue running successfully until the end, but detached from
            # the parent python process, so it doesn't know that the R has finished (and wrote the results to the outfiles).
            # If this is the case - there is no need to run the synteny scores calculation again for the current ref-genome.

            # If the synteny calculation was finished successfully, the final_output directory with all the per-genome
            # output files should exist and the R_temp folder should have been removed.
            if os.path.exists(final_output_path) is True and os.path.exists(r_temp_path) is False:

                # Check the number of existing per-genome final output files
                outfiles_num = len(os.listdir(final_output_path))

                # According to the number of output files - R has finished successfully
                if outfiles_num == len(config.subsampling_lengths) + 1:
                    print("\nFound synteny calculation output files - no need to run this stage again")
                    config.genomes_dict[ref_genome]['finished_R'] = 1
                else:
                    print("\nNumber of output files doesn't match the expected - run the synteny calculation again")
                    config.genomes_dict[ref_genome]['finished_R'] = 0

            # Probably R hasn't finished successfully -> run it again
            else:
                print("\nFound no output files of the synteny calculation - run it again")
                config.genomes_dict[ref_genome]['finished_R'] = 0

        ####################################################################################
        # Step 3: Run synteny calculation using R
        if config.genomes_dict[ref_genome]['finished_R'] == 0:

            before = time.time()

            # Create a folder for R final output tables
            final_output_path = ref_genome_output_dir + config.final_output_dir

            # If the R final output dir already exists - delete it and its content
            if os.path.exists(final_output_path):
                print("\nDirectory " + final_output_path + " already exists - deleting its content")
                shutil.rmtree(final_output_path)
            # Create a new R final output dir
            try:
                os.makedirs(final_output_path)
            except OSError:
                print("\nmkdir " + final_output_path + "has failed")
                exit()

            # Create a folder for R temporary files (will be deleted in the end of the run)
            r_temp_path = ref_genome_output_dir + config.r_temp_dir

            # If the R temporary files dir already exists - delete it and its content
            if os.path.exists(r_temp_path):
                print("\nDirectory " + r_temp_path + " already exists - deleting its content")
                shutil.rmtree(r_temp_path)
            # Create a new R temporary files dir
            try:
                os.makedirs(r_temp_path)
            except OSError:
                print("\nmkdir " + r_temp_path + "has failed")
                exit()

            # Create a folder for R intermediate objects if the user has asked for it (for debugging purposes)
            if config.save_intermediate:
                intermediate_objects_path = ref_genome_output_dir + config.r_intermediate_objects_dir

                # If the R intermediate objects dir already exists - delete it and its content
                if os.path.exists(intermediate_objects_path):
                    print("\nDirectory " + intermediate_objects_path + " already exists - deleting its content")
                    shutil.rmtree(intermediate_objects_path)
                # Create a new R intermediate objects dir
                try:
                    os.makedirs(intermediate_objects_path)
                except OSError:
                    print("\nmkdir " + intermediate_objects_path + "has failed - cannot save R intermediate objects")
                    config.save_intermediate = False
            else:
                intermediate_objects_path = ""

            # Run the R script for the synteny analysis of the current reference genome
            print("\nStarting synteny analysis for genome " + ref_genome + "\n")
            if intermediate_objects_path == "":
                intermediate_objects_path = "NA"
            if config.metadata_file_path == "":
                metadata_file_path = "NA"
            else:
                metadata_file_path = config.metadata_file_path

            command = "Rscript syntracker_R_scripts/SynTracker.R" + " " + ref_genome + " " + \
                        config.sample_dictionary_table_path + " " + genome_blastdbcmd_out_dir + " " + \
                        final_output_path + " " + config.summary_output_path + " " + r_temp_path + " " + " " + \
                        intermediate_objects_path + " " + str(config.seed_num) + " " + str(config.avg_all) + " " + \
                        str(config.cpu_num) + " " + metadata_file_path
            print("\nRunning the following Rscript command:\n" + command + "\n")

            try:
                subprocess.run(["Rscript", "syntracker_R_scripts/SynTracker.R", ref_genome,
                                config.sample_dictionary_table_path,
                                genome_blastdbcmd_out_dir, final_output_path, config.summary_output_path,
                                r_temp_path, intermediate_objects_path, str(config.seed_num), str(config.avg_all),
                                str(config.cpu_num), metadata_file_path], check=True)
            except subprocess.CalledProcessError as err:
                print("\nThe following command has failed:")
                print(command)
                print(err)
                exit()
            except Exception as err:
                print("\nThe following command has failed:")
                print(command)
                print(err)
                exit()

            after = time.time()
            duration = after - before
            print("\nThe synteny scores calculation stage took " + str(duration) + " seconds.\n")

        # Mark this ref-genome as finished and write in the config.txt file that the synteny stage was finished
        print("\nThe processing of genome " + ref_genome + " completed successfully\n")
        out_param = open(config.conf_file_path, "a")
        out_param.write("- Synteny finished\n\n")
        out_param.close()
        config.genomes_dict[ref_genome]['finished_R'] = 1
        config.genomes_dict[ref_genome]['processed'] = 1


if __name__ == '__main__':
    main()





