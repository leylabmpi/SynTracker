
# SynTracker: a pipeline for tracking closely related microbial strains using genome synteny
#### The SynTracker pipeline is designated to determine the biological relatedness of conspecific strains (microbial strains of the same species) using genome synteny. It consists of two main steps:
#### A. Fragmentation of the reference genomes and execution of BLASTn search against the target genomes/metagenomes.
#### B. Calculation of average pairwise synteny scores (APSS).  

More detailed information about SynTracker's pipeline and its usage is available here:
https://github.com/leylabmpi/SynTracker/blob/master/SynTracker_Manual.pdf

## Installation

### Requirements: 
NCBI BLAST+

All the other required packages are contained in the attached conda environment (.yml file).

### Installing from source:
Download SynTrakcer's latest release from: https://github.com/leylabmpi/SynTracker/releases.

Extract the tar.gz file into the desired working-directory.

Create a new conda environment using the ‘SynTracker_env.yml’ file (located in the root directory of SynTrakcer) using the following command:
      `conda env create -f SynTracker_env.yml`

Activate the newly created environment: 
      `conda activate SynTracker_1_3`


## Input
SynTracker requires two types of data as input:

#### a.	Reference genomes: 
Reference genomes can be provided complete or as a collection of contigs. If using a number of contigs belonging to the same reference genome, all sequences should be placed in a single .fasta file. 

If using more than one reference genome (i.e., analyzing more than one species per run) all reference genome files should be located in the same directory.  

The directory `Sample_Data/Input_example/Reference_genomes/` contains an example reference genome.
#### b.	Metagenomic assemblies/genomes: 
These are the metagenomic assemblies (or assembled genomes, if those are studied) that would be compared.
These data should be organized in per-sample assembly files - i.e., all the contigs assembled from sample X would be kept in a single file. 

If genomes are to be compared, each genome will be stored in a single fasta file. 
All files should be stored in the same directory (refered below as the "target directory"). 

The directory `Sample_Data/output_example/Target_genomes/` contains a collection of target genomes, for the purpuse of self testing the instalation.
   
#### Metadata file (optional): 
In addition to the mandatory input, the user can also provide a metadata file, which contains information regarding the genomes/assemblies that should be compared. 
The metadata file should be a tab delimited file. One of the columns should contain the sample ID, which is identical to the naming of the fasta files in the "target folder".

#### Sample input:

The directory ‘Sample_Data/Sample_Input/’, which is included in the SynTracker package, contains one reference genome and a collection of target genomes. It can be used to clarify the structure of the required input files and to test the installation.

## Usage

SynTracker has two main modes of execution: 
1. **'New' mode**: a new run. In this case the user must provide the path to the reference genomes and to the target genomes. 
All the other parameters are optional (including the output directory, which is created by default under 
the working directory with the name 'Syntracker_output'). 


2. **'Continue' mode**: continue a previous run that has been terminated for some reason without having to start the process from the beginning. 
In this case, the user must provide only the path to the output folder of the run that he wish to continue. 

      Continue a previous run can be done in two ways:
      
      a. Using mode = '**continue**': the run will continue from the point it stopped within the last reference genome that has been processed without finishing successfully.
      
      b. Using mode = '**continue_all_genomes**': process all the reference genomes again, without repeating the stage in which a blast database is built from the target genomes (which can be very time-consuming in case of many targets). 
It only makes sense to use this mode with when running more than one reference genome.

### Usage examples using the provided sample data
(With the minimal required mandatory input parameters)

**A new run:**
```
python syntracker.py -target Sample_Data/Input_example/Target_genomes/ -ref Sample_Data/Input_example/Reference_genomes/ -out SynTracker_output/
```

**Continue a previous run that has been terminated:**

1. Continue from the last reference genome that has been processed without finishing successfully:
```
python syntracker.py -out SynTracker_output/ -mode continue
```

2. Process all the reference genome again without repeating the blastDB building stage 
(relevant only for datasets containing more than one reference genome):
```
python syntracker.py -out SynTracker_output/ -mode continue_all_genomes
```

### A description of all SynTracker's possible command line arguments:

```
python syntracker.py [-h] [-target target_directory_path] [-ref ref_directory_path] 
                     [-out output_directory_path] [-metadata metadata_file] 
                     [-mode 'new'/'continue'] [-cores number_of_cores] [-length region_length] 
                     [--identity blast_identity] [--coverage blast_coverage] [--no_seed] [--avg_all]

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
  
  -mode ['new'/'continue'/'continue_all_genomes']  
                    The running mode: 'new' or 'continue' (default='new'). 
                    Start a new run or continue a previous run that has been terminated.
                    'continue' mode: continue from the last reference genome that was previously processed.
                    'continue_all_genomes' mode: process all the reference genomes again, without repeating the stage in which a blast database is built from the target genomes.
  
  -cores [number_of_cores]
                    The number of cores to use for the multi-processed stages of the calculation. 
                    (Optional, by default SynTracker uses the maximal number of available cores).
  
  -length [region_length]
                    The length of the compared region. (Optional, default=5000)
  
  --identity [blast_identity]
                    Minimal blast identity (optional, default=97)
  
  --coverage [blast_coverage]
                    Minimal blast coverage (optional, default=70)
  
  --no_seed         Set no seed for the subsampling of n regions per pairwise (optional). 
                    This means that the average synteny scores may change between SynTracker runs due to the subsampling. 
                    By default, a seed=1 is set to enable reproducibility between different runs.
 
  --avg_all   
                    Create an additional output table with APSS (Average Pairwise Synteny Scores), 
                    which are based on all the available regions per each pair of samples 
                    (in addition to the output tables, based on the subsampling of n regions).                 
```

## Output

#### Output per genome:
For each given reference genome, SynTracker outputs two types of tables, both include pairwise specific information.
All the output files for a certain reference genome are located under the directory `[genome_name]/final_output/`.

The table `[genome name]_synteny_scores_per_region.csv` contains the raw results obtained by the comparison of each two homologous genomic 
regions in each two metagenomes in which they were detected (or genomes, if those are being compared).

The second type of output tables, `[genome name]_avg_synteny_scores_[subsampling length]_regions.csv`, gives the APSS 
(Average Pairwise Synteny Score) that was calculated by subsampling N regions per pair of samples
from the overall regions that appear in the raw table (detailed above). 
By default, N equals to 40, 60, 80, 100, 200 regions per pair of samples.

In case the user has applied the --avg_all option, an additional table, 
named `[genome name]_avg_synteny_scores_all_regions.csv` is created too. In this table, the APSS are calculated 
using all the available regions per each pair of samples. 

#### Summary output (all genomes together):
Syntracker also creates the same output tables mentioned above for all the references genomes combined together. 
These summary output files are located under the directory `summary_output/`.

The raw (per-region synteny scores) table is called `synteny_scores_per_region.csv`. 

The tables containing the APSS in different subsampling lengths are called `avg_synteny_scores_[subsampling length]_regions.csv`.

The table containing the APSS using all regions (in case of applying the --avg_all option) 
is called `avg_synteny_scores_all_regions.csv`.

#### Sample output:

The directory ‘Sample_Data/Output_example/’, which is included in the SynTracker package, contains the output of a SynTracker’s run, using the sample input data with default parameters. It can facilitate the user in better understanding the structure of the output directories and files.

## Citation

If you use SynTracker please cite:  

**Strain tracking in complex microbiomes using synteny analysis reveals per-species modes of evolution.**  
Enav H, Paz I and Ley RE.  
Nature Biotechnology (2024). DOI: https://doi.org/10.1038/s41587-024-02276-2
