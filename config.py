# Paths
working_dir = ""
input_target_dir = ""  # Should be given by the user
input_ref_dir = ""  # Should be given by the user
output_dir = "Syntracker_output/"  # Default output dir
main_output_path = ""
conf_file = "config.txt"
old_conf_file = "config_old.txt"
conf_file_path = ""
old_conf_file_path = ""
logfile = "SynTracker_log.txt"
logfile_path = ""
central_regions_dir = "central_regions/"
combined_output_dir = "combined_targets/"
summary_output_dir = "summary_output/"
summary_output_path = ""
output_summary_file = "synteny_scores_per_region.csv"
output_summary_file_path = ""
final_output_dir = "final_output/"
r_temp_dir = "R_temp/"
r_intermediate_objects_dir = "R_intermediate_objects/"
combined_output_path = ""
combined_renamed_genomes = "combined_renamed_genomes.fasta"
combined_renamed_genomes_file_path = ""
dictionary_table_full = "names_dictionary_full.tab"
dictionary_table_full_path = ""
sample_dictionary_table = "sample_names_dictionary.tab"
sample_dictionary_table_path = ""
blast_db_dir = "blastDB/"
blast_db_path = ""
blast_db_file = "GroupsDB"
blast_db_file_path = ""
blast_out_dir = "blast_output/"
blastdbcmd_out_dir = "blastdbcmd_output/"
is_metadata = False
metadata_file_path = ""

# Central regions related parameters
region_length = 1000
flanking_length = 2000
minimal_flanking_length = flanking_length * 0.9
full_length = region_length + flanking_length * 2
minimal_full_length = full_length * 0.9

# BLAST related parameters
minimal_coverage = 70
minimal_identity = 97
blast_num_threads = 2
minimal_hits_num = 2

# Job-related parameters
cpu_num = 8

# R - related parameters
save_intermediate = False
is_set_seed = True  # Whether to set a seed for the subsampling process (to have same results between different runs)
seed_num = 1
subsampling_lengths = [40, 60, 80, 100, 200]
subsampled_regions_file_names = []
for i in range(len(subsampling_lengths)):
    subsampled_regions_file_names.append("avg_synteny_scores_" + str(subsampling_lengths[i]) + "_regions.csv")
avg_all = False  # Whether to add non-subsampled output (average all the regions per pair of samples)
avg_all_file_name = "avg_synteny_scores_all_regions.csv"

# Run related parameters
running_mode = "new"  # Mode can be 'new' or 'continue'
complete_target_merge = False

# A dictionary containing all the ref-genomes names, input files paths and an indication whether they were finished
# being processed or not (0 - not processed yet, 1 - finished being processed)
# For example: genomes_dict['E_coli_K-12_MG1655']['input_file'] = 'full_path/E_coli_K-12_MG1655.fasta'
#              genomes_dict['E_coli_K-12_MG1655']['processed'] = 1
genomes_dict = {}
run_genomes_list = []  # a list of the genomes that should be processed in the current run
