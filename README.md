
# SynTracker: a pipeline for tracking closely related microbial strains using genome synteny
#### The SynTracker pipeline is designated to determine the biological relatedness of conspecific strains (microbial strains of the same species) using genome synteny. It consists of two main steps:
#### A. Fragmentation of the reference genomes and execution of BLASTn search against the target genomes/metagenomes.
#### B. Calculation of average pairwise synteny scores (APSS).  

## Requirements: 
NCBI BLAST+

All the other required packages are contained in the attached conda environment (.yml file).

## Installation:
### From source:
Download SynTrakcer latest release from: https://github.com/leylabmpi/SynTracker/releases.

Extract the tar.gz file into the desired working-directory.

Create a new conda environment using the ‘SynTracker_env.yml’ file (located in the root directory of SynTrakcer) using the following command:
      `conda env create -f SynTracker_env.yml`

Activate the newly created environment: 
      `conda activate SynTracker_1_2_0`


## Input:
SynTracker requires  three types of data as input:  
#### a.	Reference genomes: 
Reference genomes can be provided complete or as a collection of contigs. If using a number of contigs belonging to the same reference genome, all sequences should be placed in a single .fasta file. 

If using more than one reference genome (i.e., analyzing more than one species per run) all reference genome files should be located in the same directory.  

The directory `Sample_input/Reference_genomes/` contains an example reference genome.
#### b.	Metagenomic assemblies/genomes: 
These are the metagenomic assemblies (or assembled genomes, if those are studied) that would be compared.
These data should be organized in per-sample assembly files - i.e., all the contigs assembled from sample X would be kept in a single file. 

If genomes are to be compared, each genome will be stored in a single fasta file. 
All files should be stored in the same directory (refered below as the "target directory"). 

The directory `Sample_input/Target_genomes/` contains a collection of target genomes, for the purpuse of self testing the instalation.
   
#### c.	Metadata file (optional): 
The metadata file contains information regarding the genomes/assemblies to be compared. 
The metadata file should be a tab delimited file. One of the columns should contain the sample ID, which is identical to the naming of the fasta files in the "target folder".

## Usage: 

```
python syntracker.py [-h] [-target target_directory_path] [-ref ref_directory_path] [-out output_directory_path]
                     [-metadata metadata_file] [-mode 'new'/'continue'] [-cores number_of_cores] [--identity blast_identity]
                     [--coverage blast_coverage] [--length flanking_sequences_length] [--save_intermediate]
                     [--set_seed integer_for_seed]


options:
  -h, --help        show this help message and exit
  -target [target_directory_path]
                    Path of the target directory which contains metagenome assemblies or genomes
  -ref [ref_directory_path]
                    Path of the references folder containing the reference genomes
  -out [output_directory_path]
                    The path to the output directory . When running in 'new' mode (the default), this argument is optional. By
                    default a folder named 'Syntracker_output/' will be created under the current directory (if the given path
                    already exists, it will be written over). When running in 'continue' mode, it is mandatory to provide the
                    path to the output directory of the run that is requested to be continued.
  -metadata [metadata_file]
                    Path to a metadata file (optional). The file should be in CSV format and must include the sample ID.
  -mode ['new'/'continue']  
                    The running mode: 'new' or 'continue' (default='new') (Start a new run or continue a previous run that has been terminated).
  -cores [number_of_cores]
                    The number of cores to use for the parallelization of the BLAST-related stages. (Optional, default is the number of computer
                    available cores).
  --identity [blast_identity]
                    Minimal blast identity (optional, default=97)
  --coverage [blast_coverage]
                    Minimal blast coverage (optional, default=70)
  --length [flanking_sequences_length]
                    The length of the flanking sequences (from both sides of the BLAST hit). (Optional, default=2000)
  --save_intermediate   
                    Saves R intermediate data structures for debugging purposes (by default, they are not saved).
  --set_seed [integer_for_seed]
                    An integer number to set the seed for subsampling of n regions per pairwise (by default, the seed is 1).
  --no_seed         Set no seed for the subsampling of n regions per pairwise (by default, seed=1 is set).
```

### Usage examples using the provided sample data:

**A new run:**
```
python syntracker.py -target Sample_input/Target_genomes/ -ref Sample_input/Reference_genomes/ -out SynTracker_output/
```

**Continue a previous run that has been terminated:**
```
python syntracker.py -out SynTracker_output/ -mode continue
```

## Output:
SynTracker outputs two types of tables, both include pairwise specific information.
The output and intermediate files of a SynTracker run are organized by the given reference genomes.
Under the main output directory (provided by the user) there is a directory for each reference genome.

