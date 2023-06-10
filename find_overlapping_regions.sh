#!/bin/bash
# Title 


#create variables for the flags:

length="${l:=2000}"
minimal_coverage="${c:=70}"
minimal_identity="${i:=97}"

while getopts t:r:o:c:i:l: flag
do
    case "${flag}" in
        t) target=${OPTARG%/};; #metagenome assemlies/genomes - used to create the blast DB
        r) references=${OPTARG%/};; # folder with reference genomes 
        o) output=${OPTARG%/};; #output and temp folders will be written to this path
        i) minimal_identity=${OPTARG};; #minimal blast identity
        c) minimal_coverage=${OPTARG};; #blast qcov_hsp_perc
        l) length=${OPTARG};; #length of flanking sequences (on the retrieved traget sequence)
    esac
done

output=$(readlink -f "$output") #get the full path of the output folder
references="$references/"
target="$target/"
# step 0: run the python script first, to find the "central_regions"
rm -rf "$output"
mkdir "$output"
mkdir "$output/tmp"
central_regions="$output/central_regions/"

python central_regions.py $references $central_regions


# step 1: run script to:
#	a. change sample names (i.e., assemblies/genomes) to Sample.xxx
#	b. change fasta headers to contig.xxx
#	c. merge fasta files to one file
#	d. create a table with old and new names. 

echo $'starting merging metagnome-assemblies/genome files\n'

combined_output="$output/combined_targets/"
python old_to_new_names.py $target $combined_output

echo $'\nMerging metagnome-assemblies/genome files is complete'


# step 2:  Make a blast database, from all the contigs. 
#          Rerunning this step overrides the previous blast DB - pay attention!
newdir="$output/blastDB"
rm -rf "$newdir"
mkdir "$newdir"
makeblastdb -in  "$combined_output/combined_renamed_genomes.fasta" -out "$newdir/GroupsDB" -title "GroupsDB" -parse_seqids -dbtype nucl

echo $'Making blast DB is complete\n'

# Step 3: #run blast per genome => per region
#a. make output folder
#b. for each genome, make an output folder
#c. for each region run blast

newdir1="$output/blast_output"
rm -rf "$newdir1"
mkdir "$newdir1"

#run blast 

for folder in "$central_regions"/*; do
	rm -rf "$newdir1/"$(basename "${folder%.*}")""
	mkdir "$newdir1/"$(basename "${folder%.*}")""
	for query in "$folder"/*.fasta; do
	blastn -query "$query" -db "$newdir/GroupsDB" -out "$newdir1/"$(basename "${folder%.*}")"/$(basename "${query%.*}").tab"  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand"  -max_target_seqs 10000 -perc_identity "$minimal_identity" -qcov_hsp_perc "$minimal_coverage" -num_threads 8;
	done
done

echo $'BLAST search is complete!\n'




# Step 4: delete blast output files with less than 2 hits
for folder in "$newdir1/"*; do
	for seqfolder in "$folder"/*; do
	find "$seqfolder" -type f -exec awk -v x=2 'NR==x{exit 1}' {} \; -exec rm -f {} \;
done
done

echo $'Deleting low-hit regions is complete\n'



# Step 5: Use the blast output file to:
# a. check the orientation of the hit (minus or plus strand). If minus, change the order of start and end of hit. 
# b. take the 2/9/10/13 columns of the output (query id/start/end/strand). add/remove x base pairs on each side of the hit (columns 9-10).
# c. pipe to xargs, execute "blastdbcmd", to fetch the sequence in each hit, with x basepairs on each side.
# d. write to file to be analyzed in R synteny pipeline

set --
newdir2="$output/blastcmddb_output"
rm -rf "$newdir2"
mkdir "$newdir2"
export blastDB="$newdir/GroupsDB"
#export length=$length

for folder in "$newdir1"/*; do
	LogFile="$newdir2"/$(basename "${folder%.*}").log 
	echo $"Retrieval of genomic regions: log file is $LogFile"
	mkdir "$newdir2/$(basename "${folder%.*}")"
	for tabular in "$folder"/*.tab; do 
	awk -v "len=$length" -v "outlog=$LogFile" '{if ($13 == "minus") print $2"\t"$10-len"\t"$9+len"\t"$13"\t"; else print $2"\t"$9-len"\t"$10+len"\t"$13}' "$tabular" | xargs -n 4 sh -c 'blastdbcmd -db "$blastDB" -entry "$0" -range "$1"-"$2" -strand "$3"  -outfmt %f -logfile outlog;'  >  "$newdir2"/$(basename "${folder%.*}")/$(basename "${tabular%.*}").fasta;
	done
done 

echo $'Retrieval of flanking regions is complete!\n'

echo $'First step is complete!\n'
