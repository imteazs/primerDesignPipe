import pandas as pd
from Bio import SeqIO
import primer3 as primer
import argparse


def primerDesign(seqRec, start, length):
    """
    Take seq record and run through primer3 pipeline along with data to
    determine window where primers will be designed.
    :param seqRec: record from the fasta file
    :param start: base position of interest
    :param length: length of the window of interest
    :return: primer3 data in the form of a dictionary
    """
    inputDict = {
        'SEQUENCE_ID': seqRec.id,
        'SEQUENCE_TEMPLATE': str(seqRec.seq),
        'SEQUENCE_INCLUDED_REGION': [start, length]
    }
    designSettings = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 17,
        'PRIMER_MAX_SIZE': 35,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 65.0,
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
        'PRIMER_PRODUCT_SIZE_RANGE': [[90, 100]]
    }
    oligo = primer.designPrimers(inputDict, designSettings)
    return oligo


def genDF(primerList):
    """
    :param primerList: primer3 dictionary object
    :return: preprocessed dataframe to be ready for further processing
    """
    data = pd.DataFrame(primerList)
    data = data.T
    finalData = pd.DataFrame(columns=['oligo', 'id', 'type', 'oligo_size', 'product_size', 'Tm', 'GC_pct'])
    for index, row in data.iterrows():
        rowname_split = row.name.split('_')

        if 'sequence' in row.name.lower():
            # ['PRIMER', 'LEFT', '0', 'SEQUENCE']
            nrow = {'oligo' : data[0][index], 'id': rowname_split[2], 'type': rowname_split[1]}
            finalData = finalData.append(nrow, ignore_index=True)
        elif 'tm' in row.name.lower():
            idx = finalData.index[
                (finalData['type'] == rowname_split[1]) & (finalData['id'] == rowname_split[2])].tolist()
            finalData['Tm'][idx] = data[0][index]

        elif 'gc_percent' in row.name.lower():
            idx = finalData.index[
                (finalData['type'] == rowname_split[1]) & (finalData['id'] == rowname_split[2])].tolist()
            finalData['GC_pct'][idx] = data[0][index]

        elif 'pair' in row.name.lower():
            idx = finalData.index[finalData['id'] == rowname_split[2]].tolist()
            finalData['product_size'][idx] = data[0][index]

    finalData['oligo_size'] = finalData['oligo'].str.len()
    print(finalData)
    return finalData


def addCalc(primDF):
    restrict_enz = {
        'CviAII': ['CATG', 'GTAC'],
        'FatI': ['CATG', 'GTAC'],
        'Hpy188III' : ['TCNNGA', 'AGNNCT'],
        'NlaIII': ['CATG', 'GTAC'],
        'CviQI': ['GTAC', 'CATG'],
        'RsaI': ['GTAC', 'CATG']
    }
    return None


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
    parser.add_argument('--mv_cations', help='Concentration of monovalent cations in mM', type=float)
    parser.add_argument('--dv_cations', help='concentration of divalent cations in mM', type=float)
    parser.add_argument('--dntp', help='concentration of dntps in uM', type=float)
    parser.add_argument('--DNA', help='concentration of DNA in nM', type=float)
    parser.add_argument('--anneal_temp', help='annealing temp in C', type=float)

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
        resultdf = genDF(design)
