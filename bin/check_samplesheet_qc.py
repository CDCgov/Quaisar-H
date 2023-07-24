#!/usr/bin/env python

# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

import os
import sys
import errno
import argparse
import gzip
import re

def parse_args(args=None):
    Description = "Reformat cdcgov/phoenix samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    sample,fastq_1,fastq_2,fastp_pass_json,fastp_failed_json,spades,mlst,quast,amrfinderplus
    SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz,passed_fastp.json,failed.fastp.json,spades.fasta,mlst.tsv,quast.tsv,amrfinderplus.tsv

    sample             A sample name for the input
    fastq_1            R1 of reads run through Fastp
    fastq_2            R2 of reads run through Fastp
    fastp_pass_json    JSON output from initial Fastp run
    fastp_failed_json  JSON output from rerun of Fastp on failed reads
    spades             Assembly created by SPAdes
    mlst               TSV output from mlst tool
    quast              TSV report generated from quast
    amrfinderplus      TSV report generated from amrfinderplus
    """

    sample_mapping_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 2
        HEADER = ['sample', 'fastq_1', 'fastq_2', 'fastp_pass_json', 'fastp_failed_json', 'spades', 'mlst', 'quast', 'amrfinderplus']

        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, fastq_1, fastq_2, fastp_pass_json, fastp_failed_json, spades, mlst, quast, amrfinderplus = lspl[: len(HEADER)]
            sample = sample.replace(" ", "_")
            if not sample:
                print_error("Sample entry has not been specified!", "Line", line)

            ## Check FastQ file extension
            for fastq in [fastq_1, fastq_2]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"): # If file is not gzipped then gzip it. 
                        fastq_gz = fastq + ".gz"
                        with open(fastq, "rb") as f_in:
                            with gzip.open(fastq_gz, 'wb') as f_out: 
                                f_out.writelines(f_in)
                        print("FastQ file does not have extension '.fastq.gz' or '.fq.gz'! Zipping file.",
                            "Line",
                            line,
                        )

            ## Auto-detect paired-end/single-end
            sample_info = []  ## [single_end, fastq_1, fastq_2]
            if sample and fastq_1 and fastq_2:  ## Paired-end short reads
                sample_info = ["0", fastq_1, fastq_2, fastp_pass_json, fastp_failed_json, spades, mlst, quast, amrfinderplus]
            elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
                sample_info = ["1", fastq_1, fastq_2, fastp_pass_json, fastp_failed_json, spades, mlst, quast, amrfinderplus]
            else:
                print_error("Invalid combination of columns provided!", "Line", line)

            ## Create sample mapping dictionary = { sample: [ single_end, fastq_1, fastq_2 ] }
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[sample].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["sample", "single_end", "fastq_1", "fastq_2", "fastp_pass_json", "fastp_failed_json", "spades", "mlst", "quast", "amrfinderplus"]) + "\n")
            for sample in sorted(sample_mapping_dict.keys()):

                ## Check that multiple runs of the same sample are of the same datatype
                if not all(x[0] == sample_mapping_dict[sample][0][0] for x in sample_mapping_dict[sample]):
                    print_error("Multiple runs of a sample must be of the same datatype!", "Sample: {}".format(sample))

#                for idx, val in enumerate(sample_mapping_dict[sample]):
#                    fout.write(",".join(["{}_T{}".format(sample, idx + 1)] + val) + "\n")
                for idx, val in enumerate(sample_mapping_dict[sample]):
                    if not val[1].endswith(".gz"): # check that forward read is a gzip file
                        val[1] = re.sub(".fastq$", ".fastq.gz", val[1])
                        val[1] = re.sub(".fq$", ".fq.gz", val[1])
                    if not val[2].endswith(".gz"): # check that reverse read is a gzip file
                        val[2] = re.sub(".fastq$", ".fastq.gz", val[2])
                        val[2] = re.sub(".fq$", ".fq.gz", val[2])
                    fout.write(",".join(["{}".format(sample)] + val) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
