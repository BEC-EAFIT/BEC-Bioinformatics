#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Short utility to separate the sequences in a SILVA database into the
three domains. Originally constructed to separate
SILVA_138.1_SSURef_NR99_tax_silva.fasta (obtained from
https://www.arb-silva.de/)

Nicolás Pinel
npinelp@eafit.edu.co
Universidad EAFIT
Medellín, Colombia

March 31, 2021

"""

import re
import argparse
from Bio import SeqIO


####################
# Global variables #
####################

sequences = {}


#############
# Functions #
#############

def construct_argparser():
    """Create an argparse object to collect command line arguments."""

    parser = argparse.ArgumentParser(prog="split_silva.py")

    parser.add_argument(
        "-s",
        "--silva",
        help="fasta-formated file containing the SILVA database",
        required=True,
    )

    command_line_arguments = parser.parse_args()
    return command_line_arguments


def write_sequences(silva_file):
    """Creates the output sequence files, one per domain."""

    for items in sequences.items():
        domain = items[0]
        f_ext = re.search(r"\.\w+$", silva_file).group()
        o_file = silva_file.replace(f_ext, "_" + domain + f_ext)

        SeqIO.write(items[1], o_file, "fasta")


def main():
    """Driving code."""

    args = construct_argparser()

    with open(args.silva, "r") as handle:

        for record in SeqIO.parse(handle, "fasta"):

            record.seq = record.seq.back_transcribe()
            description = record.description.split(" ")

            domain = re.search(r"^(\w+)\;", description[1]).group(1)

            if domain in sequences.keys():
                sequences[domain].append(record)
            else:
                sequences[domain] = [record]

    write_sequences(args.silva)


if __name__ == "__main__":
    main()
