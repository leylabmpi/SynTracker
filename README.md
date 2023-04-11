
# SynTracker: a pipeline to track closely related microbial strains using genome synteny
### SynTracker is a pipeline to determine the biological relatedness of conspecific strains (microbial strains of the same species) using genome synteny.
 

## Requirements: 
Most required software/packages are contained in the attached conda environment (.yml file).
Additional software:
NCBI BLAST+


## Input:
SynTracker requires as input three types of data:  
#### a.	Reference genome: 
Reference genomes could be provided complete or as a collection of contigs. If using a number of contigs belonging to the same reference genome, all sequences should be placed in a single .fasta file. 
If using more than one reference genome (i.e., analyzing more than one species per run) all reference genome files should be located in the same directory.  
#### b.	Metagenomic assemblies/genomes: 
These data are the metagenomic assemblies (or assembled genomes, if those are studied) that would be compared.
These data should be organized in per-sample assembly files - i.e., all the contigs assembled from sample X would be kept in a single file. If genomes are studied, each genome will be stored in a unique .fasta file. 
Ideally, file naming should be the sample or isolate isolate ID, with no special characters. 
All files should be stored in the same directory (refered below as the "target directory").    
#### c.	Metadata file (optional): 
The metadata file contains information regarding the genomes/assemblies to be compared. 
The metadata file should be a tab delimited file. One of the columns should contains the sample ID, which is identical to the naming of the fasta files in the "target folder".

## Running: 

The Syntracker pipeline is composed of two main parts, each executed by running a dedicated script:  
The first part fragments the reference genomes to yield the "central regions" and run a BLASTn search, while that in the second part the average pairwise synteny score is calculated. 

### a. First part: 
Extracting 
```
./find_overlapping_regions.sh  -t -r -o -i -c -l 
```
        -t : path of the target directory, contains metagenome assemlies/genomes
        -r : path of the references folder containing the reference genomes 
        -o : path of the output directory. 
        -i : optional, minimal blast identity, default is 97%
        -c : optinal, minimal BLAST coverage, default is 70%
        -l : optinal, the length of flanking sequences (on the BLAST hit), default is 2000 bp. 

#### b.	Second part:
Synteny analysis and calculation of Average Pairwise Synteny Scores (APSS)
```
	This part is executed by running the R script:  
Rscript SynTracker_revised_testing_version.R  [working directory] [input library] [number of cores] [save RDS] [set.seed] [metadata file]
```
```
working dir : location of the directory where the R scripts are located
input library : The path of the /"blastcmddb_output" directory within the output directory specified in step a. 
number of cores : 
save RDS: should RDS image of be saved for the run (yes/no : "--intermediate"/"--no_indermediate"). If not provided the script fails
set.seed : should the subsampling of n regions per pairwise be random or not ("--use.setseed", "--setseed.off").  If not provided the script fails. For the sake of repreducibility it is recommended to specify --use.setseed
metadata : metadata file, should include the sample ID, and any other relevant fields. Optional.  
```