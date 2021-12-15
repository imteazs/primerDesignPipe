import pandas as pd
import primer3
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
    Take the primer3 data and parse through it and transfer select data to a
    :param primerList: primer3 dictionary object
    :return: preprocessed dataframe to be ready for further processing
    """
    data = primerList
    finalData = pd.DataFrame(columns=['id', 'oligo', 'type', 'position', 'primer_len', 'amplicon_size', 'GC_pct'])

    for key, value in data.items():
        rowname_split = key.split('_')

        if 'sequence' in key.lower():
            # ['PRIMER', 'LEFT', '0', 'SEQUENCE']
            nrow = {'oligo' : value, 'id': rowname_split[2], 'type': rowname_split[1]}
            finalData = finalData.append(nrow, ignore_index=True)

        elif 'gc_percent' in key.lower():
            idx = finalData.index[
                (finalData['type'] == rowname_split[1]) & (finalData['id'] == rowname_split[2])].tolist()
            finalData['GC_pct'][idx] = value

        elif 'pair' in key.lower():
            idx = finalData.index[finalData['id'] == rowname_split[2]].tolist()
            finalData['amplicon_size'][idx] = value

        elif isinstance(value, tuple):
            idx = finalData.index[
                (finalData['type'] == rowname_split[1]) & (finalData['id'] == rowname_split[2])].tolist()
            finalData['position'][idx] = value[0]
            finalData['primer_len'][idx] = value[1]

    return finalData


def addCalc(primDF):
    """
    start calculating and compiling all the metrics that will be used for primer picking
    :param primDF: dataframe from the primer pipeline
    :return: primer after calculations
    """
    '''
    Calculate the homodimer of each oligo and get their deltaG and Tm in C
    '''
    primDF['Tm_calc'] = primDF['oligo'].apply(primer3.calcTm, mv_conc=50, dv_conc=4.7, dntp_conc=0.00095, dna_conc=200)
    primDF['Homodimer_thermo'] = primDF['oligo'].apply(primer3.calcHomodimer, mv_conc=50, dv_conc=4.7, dntp_conc=0.00095, dna_conc=200)
    primDF['Homodimer_delG_kcal/mol'] = [i.dg/1000 for i in primDF['Homodimer_thermo'].values]
    primDF['Homodimer_Tm_C'] = [i.tm for i in primDF['Homodimer_thermo'].values]

    '''
    Check the oligos for the whether the restriction sites of these enzymes actually exist
    future plan: drop down menu in shiny app and collapse it all into a list
    '''
    restrict_enz = {
        'CviAII': ['CATG'],
        'FatI': ['CATG'],
        'Hpy188III': ['TCNNGA'],
        'NlaIII': ['CATG'],
        'CviQI': ['GTAC'],
        'RsaI': ['GTAC']
    }
    for key, value in restrict_enz.items():
        if 'Hpy188III' in key:
            primDF[key] = primDF['oligo'].str.contains('TC..GA')
        else:
            for item in value:
                primDF[key] = primDF['oligo'].str.contains(item, regex=False)
    primDF['restric_enz_hit'] = primDF[restrict_enz.keys()].any(axis=1)

    '''
    Count how many Gs and Cs in the last bases of the oligos
    '''
    primDF['last_5_base'] = primDF['oligo'].str.slice(start=-5)
    primDF['3_prime_GC_count_last_5_oligos'] = primDF['last_5_base'].str.count('G') + primDF['last_5_base'].str.count('C')

    '''drop the processing columns'''
    primDF = primDF.drop(columns=['Homodimer_thermo', 'last_5_base'])
    primDF = primDF.drop(columns=restrict_enz.keys())

    '''
    Split the dataframe up into primer type
    '''
    left = primDF[primDF['type'] == 'LEFT']
    left.drop(columns=['type'], inplace=True)

    right = primDF[primDF['type'] == 'RIGHT']
    right.drop(columns=['type'], inplace=True)

    internal = primDF[primDF['type'] == 'INTERNAL']
    internal.drop(columns=['type'], inplace=True)

    '''Merge all the data into one assay id per row'''
    mergedf = left.merge(right, on=['id'], suffixes=['_left', '_right']).merge(internal, on=['id'], suffixes=['_left','_internal'])

    '''
    Calculate the heterodimers ie forward and reverse and probes deltaG and Tm in C
    '''
    #left and right
    mergedf['left_right_heterodimer_thermo'] = [primer3.calcHeterodimer(i, j, mv_conc=50, dv_conc=4.7, dntp_conc=0.00095, dna_conc=200)
                                                for i, j in mergedf[['oligo_left', 'oligo_right']].values]
    mergedf['left_right_heterodimer_kcal/mol'] = [i.dg/1000 for i in mergedf['left_right_heterodimer_thermo'].values]
    mergedf['left_right_heterodimer_Tm_C'] = [i.tm for i in mergedf['left_right_heterodimer_thermo'].values]

    #left and probe
    mergedf['left_internal_heterodimer_thermo'] = [primer3.calcHeterodimer(i, j, mv_conc=50, dv_conc=4.7, dntp_conc=0.00095, dna_conc=200)
                                                   for i, j in mergedf[['oligo_left', 'oligo']].values]
    mergedf['left_internal_heterodimer_kcal/mol'] = [i.dg / 1000 for i in mergedf['left_internal_heterodimer_thermo'].values]
    mergedf['left_internal_heterodimer_Tm_C'] = [i.tm for i in mergedf['left_internal_heterodimer_thermo'].values]

    #right and probe
    mergedf['internal_right_heterodimer_thermo'] = [primer3.calcHeterodimer(i, j, mv_conc=50, dv_conc=4.7, dntp_conc=0.00095, dna_conc=200)
                                                    for i, j in mergedf[['oligo', 'oligo_right']].values]
    mergedf['right_internal_kcal/mol'] = [i.dg / 1000 for i in mergedf['internal_right_heterodimer_thermo'].values]
    mergedf['right_internal_Tm_C'] = [i.tm for i in mergedf['internal_right_heterodimer_thermo'].values]

    mergedf.drop(columns=['left_right_heterodimer_thermo', 'left_internal_heterodimer_thermo',
                          'internal_right_heterodimer_thermo'], inplace=True)

    return mergedf


def pickPrimer(calcdf):
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
    parser.add_argument('--dntp', help='concentration of dntps in mM [according to primer3 python docs]', type=float)
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
        selectdf = addCalc(resultdf)

