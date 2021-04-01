#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Batch process that applies the script "primer_pair_profiling.py"
(which must be in the same directory) to a list of files specified
in this script, collects all outputs. It can also extract all
sequence pair reads into separate file, as in the original script.

Nicolás Pinel
npinelp@eafit.edu.co
Universidad EAFIT
Medellín, Colombia

April 1, 2021

"""

import sys
import argparse
import subprocess
import logging
from pathlib import Path
from collections import defaultdict


####################
# Global variables #
####################

primer_pairs = defaultdict(dict)
logger = logging.getLogger()


###########
# Classes #
###########

class LogFormatter(logging.Formatter):
    """Specifies a custom format for the logging messages."""

    def format(self, record):
        self._style._fmt = "%(asctime)s (%(levelname)s): %(message)s"
        self._style._datefmt = "%Y-%m-%d %H:%M:%S"
        return super().format(record)


#############
# Functions #
#############

def construct_argparser():
    """Create an argparse object to collect command line arguments."""

    parser = argparse.ArgumentParser(prog="primer_pair_profiling_batch.py")

    parser.add_argument(
        "-p",
        "--primers",
        help="fasta-formated file with primer sequences to be searched",
        required=True
    )
    parser.add_argument(
        "-d",
        "--directory",
        help="location for the forward sequences files in format fast[a|q]",
        required=True
    )
    parser.add_argument(
        "-e",
        "--extension",
        help="file extension to be identified as sequences in the directory",
        required=True
    )
    parser.add_argument(
        "-o",
        "--out_file",
        help="output file to print the tallies (defaults to stdout)",
        required=True
    )
    parser.add_argument(
        "-s",
        "--separate",
        action="store_true",
        help="output the primer pair reads into separate sequence files"
    )
    parser.add_argument(
        "-m",
        "--minreads",
        type=int,
        default=1000,
        help="minimum number of reads for the primer pair in orther to \
              separate them into files"
    )

    command_line_arguments = parser.parse_args()
    return command_line_arguments


def get_files(seq_dir, ext):
    """Collect from directory all the forward sequences files."""

    seq_dir = Path(seq_dir)

    seq_files = list(seq_dir.glob(f"*R1*.{ext}"))
    logger.info(f"Found {len(seq_files)} files.")

    return seq_files


def tally_primers(args):
    """Parse each file pair and tally primer pair frequency."""

    files = get_files(args.directory, args.extension)
    file_names = []  # To keep track of all the (base) file names for printing.

    for f in files:
        file_name = f.name
        logger.info(f"Processing file {file_name}")
        file_name = file_name.replace("." + args.extension, "")

        file_names.append(file_name)

        command = ["primer_pair_profiling.py", "-p" + args.primers,
                   "-f" + f.as_posix(), "-m" + str(args.minreads)]

        if args.separate:
            command.append("-s")

        output = subprocess.run(command, text=True, capture_output=True)

        output.check_returncode()  # If non-zero, raise a CalledProcessError.
        output = output.stdout.rstrip().split("\n")  # Reutilize the variable.

        for line in output:
            items = line.split("\t")
            primer_pairs[items[0]][file_name] = items[1]

    logger.info(f"Finished processing all files.")

    return file_names


def print_tallies(file_names, out_file):
    """Output the tallies to file."""

    logger.info(f"Printing results to file {out_file}")

    file_names.sort()

    with open(out_file, "w") as out_handle:

        print("pair\t" + "\t".join(file_names), file=out_handle)

        primer_totals = {}
        for pair in primer_pairs.keys():
            tot = 0
            for file_n in file_names:
                if not file_n in primer_pairs[pair]:
                    primer_pairs[pair][file_n] = 0

                tot += int(primer_pairs[pair][file_n])

            primer_totals[pair] = tot

        for items in sorted(primer_totals.items(), key=lambda x: (-x[1], x[0].lower())):
            freqs = []

            for file_n in file_names:
                freqs.append(primer_pairs[items[0]][file_n])

            freqs = map(str, freqs)
            freqs = "\t".join(freqs)
            print(items[0] + "\t" + freqs, file=out_handle)


def main():
    """Driving code."""

    args = construct_argparser()

    log_handler = logging.StreamHandler(sys.stdout)
    log_handler.setFormatter(LogFormatter())
    logger.setLevel(logging.DEBUG)
    logger.addHandler(log_handler)

    logger.info(
        f"Will collect all .{args.extension} files with forward reads in {args.directory}")

    file_names = tally_primers(args)

    print_tallies(file_names, args.out_file)

    logger.info(f"Finished everything!")


if __name__ == "__main__":
    main()
