from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
import os
import csv
import re
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


def blast_per_region_process(full_path_region_file, blast_region_outfile, blastdbcmd_region_outfile,
                             blastdbcmd_region_outfile_tmp, blast_db_file_path, flanking_length,
                             minimal_flanking_length, minimal_full_length, minimal_identity, minimal_coverage,
                             num_threads):

    # Run blast for each region
    run_blastn(full_path_region_file, blast_region_outfile, blast_db_file_path, minimal_identity, minimal_coverage,
               num_threads)

    # Count the number of hits in the blast output file
    file = open(blast_region_outfile)
    hits_num = len(file.readlines())  # Count the number of lines (=num of hits) in the output file
    hits_by_sample_dict = dict()  # a dict to hold the sample names that have already been found among the blast hits

    # At least two hits -> read the file
    if hits_num >= 2:
        with open(blast_region_outfile) as file:
            reader = csv.reader(file, delimiter='\t')
            for row in reader:
                sample_name = row[0]
                start = row[1]
                end = row[2]
                strand = row[3]

                m = re.search("^(Sample\.\d+)_contig", sample_name)
                if m:
                    sample_name_only = m.group(1)

                # If the current sample is already found in the samples dictionary ->
                # there is more than one hit per sample -> ignore this sample
                if sample_name_only in hits_by_sample_dict:
                    hits_by_sample_dict[sample_name_only]["valid"] = False

                # This is the first hit of this sample
                else:
                    hits_by_sample_dict[sample_name_only] = dict()
                    hits_by_sample_dict[sample_name_only]["valid"] = True
                    hits_by_sample_dict[sample_name_only]["sample_name"] = sample_name
                    hits_by_sample_dict[sample_name_only]["start"] = start
                    hits_by_sample_dict[sample_name_only]["end"] = end
                    hits_by_sample_dict[sample_name_only]["strand"] = strand

        valid_samples = 0
        for sample in hits_by_sample_dict:
            if hits_by_sample_dict[sample]["valid"]:
                valid_samples += 1

        # Continue to run blastdbcmd only if this region has at least two valid samples (with no more than one hit)
        if valid_samples >= 2:
            for sample in hits_by_sample_dict:
                sample_name = hits_by_sample_dict[sample]["sample_name"]
                start = hits_by_sample_dict[sample]["start"]
                end = hits_by_sample_dict[sample]["end"]
                strand = hits_by_sample_dict[sample]["strand"]

                # Plus strand
                if strand == "plus":
                    flank_start = int(start) - int(flanking_length)
                    flank_end = int(end) + int(flanking_length)
                # Minus strand
                else:
                    flank_start = int(end) - int(flanking_length)
                    flank_end = int(start) + int(flanking_length)

                # If the flanking-start is a negative value
                if flank_start <= 0:
                    # If the upstream length meets the minimal flanking length, take this sequence
                    # and change the start position to 1
                    if abs(flank_start) < flanking_length - minimal_flanking_length:
                        flank_start = 1
                    # The overall upstream length is shorter than the minimal flanking length -> ignore the sequence
                    else:
                        continue

                # Run blastdbcmd to get the hit including flanking sequences from the database
                run_blastdbcmd(sample_name, str(flank_start), str(flank_end), strand, minimal_full_length,
                               blastdbcmd_region_outfile, blastdbcmd_region_outfile_tmp, blast_db_file_path)

        else:
            print("\nBLAST output file " + blast_region_outfile + " contains less than two valid samples - "
                                                                "ignoring it...\n")


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


def run_blastdbcmd(entry, start, end, strand, minimal_full_length, outfile, outfile_tmp, blast_db_file_path):
    range = start + "-" + end
    args = "-db " + blast_db_file_path + " -entry " + entry + " -range " + range + " -strand " \
           + strand + " -outfmt %f"
    command = "blastdbcmd " + args + " > " + str(outfile_tmp)

    try:
        os.system(command)
    except Exception as err:
        print("\nThe following command has failed:")
        print(command)
        print(err)
        exit()

    # Read the tmp file to check the sequence's length
    if os.path.isfile(outfile_tmp):
        seq = ""
        header = ""
        with open(outfile_tmp) as file:
            for line in file:

                # Header line
                if re.search(r"^>", line):
                    header = line

                    # Plus strand
                    if re.search(r"^>.+:\d+", line):
                        m = re.search("^>.+:(\d+)-(\d+)\n", line)
                        if m:
                            start = m.group(1)
                            end = m.group(2)

                    # Minus strand
                    elif re.search(r"^>.+:c\d+", line):
                        m = re.search("^>.+:c(\d+)-(\d+)\n", line)
                        if m:
                            end = m.group(1)
                            start = m.group(2)

                # Sequence line
                elif re.search(r"^\w+", line):
                    seq += line

        length = int(end) - int(start) + 1

        # Write this sequence into the blastdbcmd output file only if the length meets the minimal criteria
        if length >= minimal_full_length:
            out_file = open(outfile, "a")
            out_file.write(header)
            out_file.write(seq)
            out_file.close()


