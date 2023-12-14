import re
import config


################################################################################################
# This function takes one reference genome, divides it into "central regions" of a pre-defined length (default=1000 bp).
# By default, central regions are located 4 kbp apart (flanking regions of 2000 bp from each side).
# Writes a separate fasta file for each region in each contig - all located under the 'central_regions' directory
# under the ref-genome output directory
def find_central_regions(genome_name, central_regions_dir):

    # Create a dictionaries of contigs => linearized DNA sequence, for each input file
    # print 1kb regions for each ref to a matching folder.
    contig_seqs = {}

    input_file_path = config.genomes_dict[genome_name]['input_file']
    with open(input_file_path) as read_file:
        for line in read_file:
            # Find a fasta header
            if re.search(r"^>", line):
                tmp_title = re.sub(' ', '_', line)
                tmp_title = re.sub('>', '', tmp_title)
                tmp_title = tmp_title.rstrip("\n")
                contig_seqs[tmp_title] = ""
            else:
                newline = line.rstrip("\n")
                contig_seqs[tmp_title] += newline

    # For each contig, go over the sequence and create 1kb regions, 5kb apart
    # (between start positions)
    # Create file name, print sequence to file.
    for contig in contig_seqs:
        pos = 0
        while pos < (len(contig_seqs[contig]) - config.region_length):
            contig_region_file = central_regions_dir + contig + "_" + str(pos) + "_" + \
                                 str(pos + config.region_length) + ".fasta"
            outgene_file = open(contig_region_file, "w")
            outgene_file.write(">" + contig + "_" + str(pos) + ":" + str(pos + config.region_length) + "\n")
            outgene_file.write(contig_seqs[contig][pos:pos + config.region_length])
            outgene_file.write("\n")
            outgene_file.close()
            pos += config.full_length

