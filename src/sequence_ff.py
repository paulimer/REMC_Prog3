#!/usr/bin/env python3


def get_sequence_from_file(filename):
    """gets a sequence from a fasta file"""
    sequences = []
    with open(filename, "r") as file_in:
        for line in file_in:
            if not line.startswith(">"):
                sequences += [line[:-1]]
    print("Read sequence:")
    print("".join(sequences))
    return "".join(sequences)
