"""
This is a part of the SynTracker pipeline.
Title: old to new names
Author: Hagay Enav
Date: april 2023
This script takes the target genomes/assemlies, assignes new generic names to the samples (i.e., file name = sample name)
and to the contigs (replaces the original fasta header with contig.xxx).

This is done to avoid 
	a. Naming issues, mostly when special characters are included in the file name/fasta headers.
	b. Too long fasta headers, which results in abortion of the blast search  

All the input target genomes/assemblies are written to one combined file. 
A second file holds a table with old and new Sample names, old and new contig names. 


* user should provide the path to the folder with the target genomes
* user should provide the paths of the output and table files
* execution: python old_to_new_names.py (target_genomes_folder) (output_file) (table_file)
"""


import sys
import subprocess
import re
import os
from os import listdir
import shutil

#print(sys.argv)
if (len(sys.argv) != 3):
	print("\n\tUser input should be target_genomes_folder   Output_folder\n\tExit due to improper user input.\n")
	sys.exit()
else:
	print ("old to new names: user input is OK")
	input_folder= sys.argv[1]
	output_folder= sys.argv[2]

#input_folder="/Users/henav/GitHub/Sample_data_input/Target_genomes/"
#output_file="/Users/henav/GitHub/combined_renamed_genomes.fasta"
#dictionary_table="/Users/henav/GitHub/old_to_new_names.tab"

file_counter=0
#output folder contains the dictionary table (old file/contig names => new names)
if os.path.isdir(output_folder):
    shutil.rmtree(output_folder)
os.mkdir(output_folder)
 
output_file=str(output_folder)+"combined_renamed_genomes.fasta"
dictionary_table=str(output_folder)+"names_dictionary.tab"
output=open(output_file, "a")
table=open(dictionary_table, "a")
table.write("old sample name"+"\t"+ "new sample name"+"\t"+"old contig header"+"\t"+"new contig header"+"\n")


for input_file in listdir(input_folder):
    if re.search(r"^.DS", input_file):
        continue
    file_counter+=1
    newfile_name='Sample.' + str(file_counter) 
    
    input_path=input_folder+input_file
    old_sample_name=os.path.splitext(input_file)[0]

    #tmp_file.write(outfile + "\t" + input_file + "\n")
    contig_counter=0 #zero the contig counter
    infile=open(input_path, "r") #read the file
    
    for line in infile:
        if re.search(r"^>", line):
            contig_counter+=1   
            line=line.rstrip('\n')
            # Replace the target string
            newline = ">" + newfile_name + "_contig." + str(contig_counter) +"\n"    
            output.write(newline)
            table.write(old_sample_name+"\t"+newfile_name+"\t"+line+"\t"+"contig."+str(contig_counter)+"\n")
        else:
            output.write(line)
    infile.close()
table.close()
output.close()
