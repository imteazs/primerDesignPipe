import pandas as pd
from Bio import SeqIO
import primer3 as primer
import argparse
import csv

def primerDesign(seqRec, start, length):
    inputDict = {
                    'SEQUENCE_ID': seqRec.id,
                    'SEQUENCE_TEMPLATE': str(seqRec.seq),
                    'SEQUENCE_INCLUDED_REGION': [start,length]
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
        print(key, ',', value)
    return oligo

def primerFilter(primerList):
    data = pd.DataFrame(primerList)
    finalData = None
    return finalData


if __name__ == "__main__":
    '''
    Creating argument parser so that command line arguments can be passed
    - First for getting path of the fasta file.
    - adding fields for included region.
    - possibly adding in a field for out directory
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--fastafile', help='Path to fasta file')
    parser.add_argument('--included_region_start', help='start position of subregion', type=int)
    parser.add_argument('--included_region_length', help='length of subregion', type=int)
    args = parser.parse_args()
    fasta = args.fastafile
    include_reg_start = args.included_region_start
    include_reg_len = args.included_region_length

    # creating SeqIO biopython object
    seq = SeqIO.parse(fasta, 'fasta')

    '''
    * This is code is designed to run bulk requests, so it can be used to design primers
      for either qPCR or sequencing applications.
    * Just provide fasta with either multiple genes or single genes. The program should
      be able to handle it.
    '''
    for record in seq:
        design = primerDesign(record, include_reg_start, include_reg_len)
        #resultdf = primerFilter(design)



