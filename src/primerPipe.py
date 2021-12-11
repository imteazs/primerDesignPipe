import pandas as pd
from Bio import SeqIO
import primer3 as primer
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastafile')
    args = parser.parse_args()
    fasta = args.fastafile

    seq = SeqIO.parse(fasta, 'fasta')

    for item in seq:
        print(item)

