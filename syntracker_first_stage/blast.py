from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
import os
import csv
import config


def make_blast_db():

    input_file = config.combined_renamed_genomes_file_path

    command = NcbimakeblastdbCommandline(input_file=input_file, out=config.blast_db_file_path, dbtype="nucl",
                                         parse_seqids=True, title=config.blast_db_file)
    print("\nExecuting the following BLAST command:")
    print(command)

    try:
        stdout, stderr = command()
        print(stdout)
    except Exception as err:
        print("\nmakeblastdb command has failed for input file: " + input_file)
        print(err)
        exit()

    print("...Done!")


def blast_per_region_process(full_path_region_file, blast_region_outfile, blastdbcmd_region_outfile, blast_db_file_path,
                             flanking_length, minimal_flanking_length, minimal_identity, minimal_coverage, num_threads):

    # Run blast for each region
    run_blastn(full_path_region_file, blast_region_outfile, blast_db_file_path, minimal_identity, minimal_coverage,
               num_threads)

    # Count the number of hits in the blast output file
    file = open(blast_region_outfile)
    hits_num = len(file.readlines())  # Count the number of lines (=num of hits) in the output file

    # If the blast output file has less than 2 hits - delete it
    if hits_num < 2:
        try:
            os.remove(blast_region_outfile)
        except OSError:
            print("\nRemoving " + blast_region_outfile + "has failed")

    # At least two hits -> read the file
    else:
        with open(blast_region_outfile) as file:
            reader = csv.reader(file, delimiter='\t')
            for row in reader:
                sample_name = row[0]
                start = row[1]
                end = row[2]
                strand = row[3]

                # Plus strand
                if strand == "plus":
                    flank_start = int(start) - int(flanking_length)
                    flank_end = int(end) + int(flanking_length)
                # Minus strand
                else:
                    flank_start = int(end) - int(flanking_length)
                    flank_end = int(start) + int(flanking_length)

                # If the flanking-start is a negative value
                if flank_start < 0:
                    # If the upstream length is at least 1800, take this sequence
                    # and change the start position to 1
                    if abs(flank_start) <= flanking_length - minimal_flanking_length:
                        flank_start = 1
                    # The overall upstream length is shorter than 1800 -> ignore the sequence
                    else:
                        continue

                # Run blastdbcmd to get the hit including flanking sequences from the database
                run_blastdbcmd(sample_name, str(flank_start), str(flank_end), strand,
                                     blastdbcmd_region_outfile, blast_db_file_path)


def run_blastn(query_file, outfile, blast_db_file_path, minimal_identity, minimal_coverage, num_threads):
    command = NcbiblastnCommandline(query=query_file, db=blast_db_file_path, out=outfile,
                                    outfmt="6 sseqid sstart send sstrand",
                                    max_target_seqs=10000, perc_identity=minimal_identity,
                                    qcov_hsp_perc=minimal_coverage, num_threads=num_threads)

    try:
        #print("Executing BLAST command:\n" + str(command))
        command()
    except Exception as err:
        print("\nThe following command has failed:")
        print(str(command))
        print(err)
        exit()


def run_blastdbcmd(entry, start, end, strand, outfile, blast_db_file_path):
    range = start + "-" + end
    args = "-db " + blast_db_file_path + " -entry " + entry + " -range " + range + " -strand " \
           + strand + " -outfmt %f"
    command = "blastdbcmd " + args + " >> " + str(outfile)

    try:
        os.system(command)
    except Exception as err:
        print("\nThe following command has failed:")
        print(command)
        print(err)
        exit()

