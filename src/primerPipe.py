import pandas as pd
from Bio import SeqIO
import primer3 as primer
import argparse

def primerDesign(seqRec):


if __name__ == "__main__":
    '''
    Creating argument parser so that command line arguments can be passed
    - First for getting path of the fasta file
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastafile')
    args = parser.parse_args()
    fasta = args.fastafile

    # creating SeqIO biopython object
    seq = SeqIO.parse(fasta, 'fasta')

    for record in seq:
        primer.designPrimers(record)

