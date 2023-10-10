
# SynTracker: a pipeline to track closely related microbial strains using genome synteny
### SynTracker is a pipeline to determine the biological relatedness of conspecific strains (microbial strains of the same species) using genome synteny.
 

## Requirements: 
Most required packages are contained in the attached conda environment (.yml file).
Additional requirements:
NCBI BLAST+

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
#### a.	Reference genome: 
Reference genomes could be provided complete or as a collection of contigs. If using a number of contigs belonging to the same reference genome, all sequences should be placed in a single .fasta file. 
If using more than one reference genome (i.e., analyzing more than one species per run) all reference genome files should be located in the same directory.  
**The directory `Sample_input/Reference_genomes/` contains an example reference genome.**
#### b.	Metagenomic assemblies/genomes: 
These data are the metagenomic assemblies (or assembled genomes, if those are studied) that would be compared.
These data should be organized in per-sample assembly files - i.e., all the contigs assembled from sample X would be kept in a single file. If genomes are to be compared, each genome will be stored in a single fasta file. 
All files should be stored in the same directory (refered below as the "target directory"). 
**The directory `Sample_input/Target_genomes/` contains a collection of target genomes, for the purpuse of self testing the instalation.**
   
#### c.	Metadata file (optional): 
The metadata file contains information regarding the genomes/assemblies to be compared. 
The metadata file should be a tab delimited file. One of the columns should contain the sample ID, which is identical to the naming of the fasta files in the "target folder".

## Usage: 

The SynTracker pipeline is divided to two main steps:
a. Fragmentation of the reference genomes and execution of BLASTn search against the target genomes/metagenomes.
b. Calculation of average pairwise synteny scores (APSS).  
Each step is executed by a single command. 

### a. First command: 
```
python syntracker.py [-h] [-target target_directory_path] [-ref ref_directory_path] [-out output_directory_path] [-mode [new|continue]]
                     [-identity blast_identity] [-coverage blast_coverage] [-length flanking_sequences_length] [-cores number_of_cores] 


options:
  -h, --help        show this help message and exit
  -target target_directory_path
                    Path of the target directory which contains metagenome assemblies or genomes
  -ref ref_directory_path
                    Path of the references folder containing the reference genomes
  -out output_directory_path
                    path of the output directory (optional, by default a folder named Syntracker_output/ will be created under the current
                    directory). IMPORTANT: if this directory already exists, it will be written over!!!
  -mode [new|continue]  
                    The running mode: 'new' or 'continue' (default='new') (Start a new run or continue a previous run that has been stopped).
  -identity blast_identity
                    Minimal blast identity (optional, default=97)
  -coverage blast_coverage
                    Minimal blast coverage (optional, default=70)
  -length flanking_sequences_length
                    The length of the flanking sequences (from both sides of the BLAST hit). (Optional, default=2000)
  -cores number_of_cores
                    The number of cores to use for the parallelization of the BLAST-related stages. (Optional, default is the number of computer
                    available cores).
```

**Example using the provided input data:**
`python syntracker.py -target Sample_input/Target_genomes -ref Sample_input/Reference_genomes/ -out first_step_output/`

#### b.	Second command:
**From the folder in which the SynTracker scripts are located, execute:**
```
Rscript SynTracker.R  [input directory] [output directory] [number of cores] [save RDS] [set.seed] [metadata file]
```
```

input directory: The path of the output directory specified in step a. 
output directory: The path of the final output directory.
number of cores : 
save RDS: should RDS image of be saved for the run (yes/no : "--intermediate"/"--no_indermediate"). If not provided the script fails
set.seed : should the subsampling of n regions per pairwise be random or not ("--use.setseed", "--setseed.off").  If not provided the script fails. #IMPORTANT! This flag temporarily does not ifluence the run!
metadata : metadata file, should include the sample ID, and any other relevant fields. Optional.  
```
**Example using the provided input data:**

`Rscript SynTracker.R first_step_output/ final_syntracker_output/ 8 --no_indermediate --use.setseed`
