# Paths
working_dir = ""
input_target_dir = ""  # Should be given by the user
input_ref_dir = ""  # Should be given by the user
output_dir = "Syntracker_output/"  # Default output dir
main_output_path = ""
conf_file = "config.txt"
conf_file_path = ""
central_regions_dir = "central_regions/"
combined_output_dir = "combined_targets/"
combined_output_path = ""
combined_renamed_genomes = "combined_renamed_genomes.fasta"
combined_renamed_genomes_file_path = ""
dictionary_table = "names_dictionary.tab"
dictionary_table_path = ""
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
minimal_flanking_length = 1800
jump_length = region_length + flanking_length * 2

# BLAST related parameters
minimal_coverage = 70
minimal_identity = 97
blast_num_threads = 2
minimal_hits_num = 2

# Job-related parameters
cpu_num = 8

# R - related parameters
save_intermediate = False
is_set_seed = False
seed_num = 1

# Run related parameters
running_mode = "new"  # Mode can be 'new' or 'continue'
complete_target_merge = False

# A dictionary containing all the ref-genomes names, input files paths and an indication whether they were finished
# being processed or not (0 - not processed yet, 1 - finished being processed)
# For example: genomes_dict['E_coli_K-12_MG1655']['input_file'] = 'full_path/E_coli_K-12_MG1655.fasta'
#              genomes_dict['E_coli_K-12_MG1655']['processed'] = 1
genomes_dict = {}
run_genomes_list = []  # a list of the genomes that should be processed in the current run
