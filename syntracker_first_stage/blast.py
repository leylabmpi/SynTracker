from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
import os
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


def run_blastn(query_file, outfile):
    command = NcbiblastnCommandline(query=query_file, db=config.blast_db_file_path, out=outfile,
                                    outfmt="6 sseqid sstart send sstrand",
                                    max_target_seqs=10000, perc_identity=config.minimal_identity,
                                    qcov_hsp_perc=config.minimal_coverage, num_threads=config.num_threads)

    try:
        command()
    except Exception as err:
        print("\nThe following command has failed:")
        print(str(command))
        print(err)
        exit()


def run_blastdbcmd(entry, start, end, strand, outfile):
    range = start + "-" + end
    args = "-db " + config.blast_db_file_path + " -entry " + entry + " -range " + range + " -strand " \
           + strand + " -outfmt %f"
    command = "blastdbcmd " + args + " >> " + str(outfile)

    try:
        os.system(command)
    except Exception as err:
        print("\nThe following command has failed:")
        print(command)
        print(err)
        exit()

