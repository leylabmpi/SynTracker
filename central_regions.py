"""
This is the first step in the SynTracker pipeline.
Title: Ref genome to central regions
Author: Hagay Enav
Date: October 2021
This script takes a number of reference genomes, and divides them to 1kbp "central regions".
Central regions are located 4 kbp apart. 

* All regions coming from the same reference genome are located in the same folder. 
* Ref genomes can contain multiple contigs per genome, as long as they are in the same fasta file.
* one file per genome, all genomes should be located in the same folder. 
* user should provide the path to the folder with the reference genome and to the output folder
* execution: python make_central_regions.py (reference_genomes_folder) (output_folder)
"""


import sys
import subprocess
import re
import os
from os import listdir

#print(sys.argv)
if (len(sys.argv) != 3):
	print("\n\tUser input should be Reference_genomes_folder_path Output_folder_path.\n\tExit due to improper user input.\n")
	sys.exit()
else:
	print ("user input is OK")
	Ref_folder= sys.argv[1]
	out_folder= sys.argv[2]


#make the output folder
os.mkdir(out_folder)
#Create a dictionaries of contigs => linearized DNA sequence, for each input file
#print 1kb regions for each ref to a matching folder. 
file_dict={}
    
for input_file in listdir(Ref_folder):
    if re.search(r"^.DS", input_file):
        continue
    input_file_path=Ref_folder + input_file
    read_file=open(input_file_path)
    acc=os.path.splitext(input_file)[0] 
    out_path=out_folder + acc
    os.mkdir(out_path)
    print(out_path)
    print(input_file)

    #generate and fill the dictionary
    contig_seqs={}
    for line in read_file:
        if re.search(r"^>", line):
            tmp_title=re.sub(' ', '_', line)
            tmp_title=re.sub('>','',tmp_title)
            tmp_title=tmp_title.rstrip("\n")
            contig_seqs[tmp_title]="" 
            #print(tmp_title)
        else:
            newline=line.rstrip("\n")        
            contig_seqs[tmp_title] = contig_seqs[tmp_title]+newline
                
    file_dict[acc]=contig_seqs


# IN each reference, in each contig, go over the sequence and create 1kb regions, 5kb apart (between start positions)
# Create file name, print sequence to file. 
for up_key in file_dict:
    #print(up_key)
    counter=0
    for key in file_dict[up_key]:
        #print(key)
        pos=0
        while pos<(len(file_dict[up_key][key])-1000):
            title= out_folder + up_key + "/" + key + "_" + str(pos) + "_" + str(pos+1000) +".fasta"
           # print(title)
            #print(contig_seqs[key][pos:pos+1000][1:50])
            outgene_file=open(title, "w")
            outgene_file.write(">"+ key + "_" + str(pos) + ":" + str(pos+1000) + "\n")
            outgene_file.write(file_dict[up_key][key][pos:pos+1000])
            outgene_file.write("\n")
            outgene_file.close()
            pos=pos+5000
            counter=counter+1
    #print(str(counter))

    