import pandas as pd
from Bio import SeqIO
import primer3 as primer
import argparse

def primerDesign(seqRec):
    inputDict = {
                    'SEQUENCE_ID': seqRec.id,
                    'SEQUENCE_TEMPLATE': str(seqRec.seq),
                    'SEQUENCE_INCLUDED_REGION': [36,342]
                }
    designSettings = {
                        'PRIMER_OPT_SIZE': 20,
                        'PRIMER_PICK_INTERNAL_OLIGO': 1,
                        'PRIMER_INTERNAL_MAX_SELF_END': 8,
                        'PRIMER_MIN_SIZE': 17,
                        'PRIMER_MAX_SIZE': 35,
                        'PRIMER_OPT_TM': 60.0,
                        'PRIMER_MIN_TM': 57.0,
                        'PRIMER_MAX_TM': 63.0,
                        'PRIMER_MIN_GC': 20.0,
                        'PRIMER_MAX_GC': 80.0,
                        'PRIMER_MAX_POLY_X': 100,
                        'PRIMER_INTERNAL_MAX_POLY_X': 100,
                        'PRIMER_SALT_MONOVALENT': 50.0,
                        'PRIMER_DNA_CONC': 50.0,
                        'PRIMER_MAX_NS_ACCEPTED': 0,
                        'PRIMER_MAX_SELF_ANY': 12,
                        'PRIMER_MAX_SELF_END': 8,
                        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                        'PRIMER_PAIR_MAX_COMPL_END': 8,
                        'PRIMER_PRODUCT_SIZE_RANGE': [[90,100]]
    }
    oligo = primer.designPrimers(inputDict, designSettings)
    for key, value in oligo.items():
        print(key, ':', value)


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
        primerDesign(record)

