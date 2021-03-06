import pandas as pd
from Bio import SeqIO
import primer3 as primer
import argparse


def primerDesign(seqRec, start, length, assay_num):
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
        'PRIMER_NUM_RETURN': assay_num,
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
    print('Assay sequences generated')
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

    print('Primer3 data transferred to dataframe')
    return finalData


def addCalc(primDF):
    """
    start calculating and compiling all the metrics that will be used for primer picking
    :param primDF: dataframe from the primer pipeline
    :return: primer after calculations
    """
    '''Calculate the Tm of all oligos in C'''
    primDF['Tm_calc'] = primDF['oligo'].apply(primer.calcTm, mv_conc=50, dv_conc=4.7, dntp_conc=0.00095, dna_conc=200)

    '''Calculate the homodimer of each oligo and get their deltaG and Tm in C'''
    primDF['Homodimer_thermo'] = primDF['oligo'].apply(primer.calcHomodimer, mv_conc=50, dv_conc=4.7, dntp_conc=0.00095, dna_conc=200)
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

    '''Find more than three repeats'''
    primDF['3plus_repeats'] = primDF['oligo'].str.contains(r'((\w)\2{3,})')

    '''
    Count how many Gs and Cs in the last bases of the oligos
    '''
    primDF['last_5_base'] = primDF['oligo'].str.slice(start=-5)
    primDF['3_prime_GC_count_last_5_oligos'] = primDF['last_5_base'].str.count('G') + primDF['last_5_base'].str.count('C')

    '''drop the processing columns'''
    primDF = primDF.drop(columns=['Homodimer_thermo', 'last_5_base'])
    primDF = primDF.drop(columns=restrict_enz.keys())

    '''Split the dataframe up into primer type'''
    left = primDF[primDF['type'] == 'LEFT']
    left.drop(columns=['type'], inplace=True)

    right = primDF[primDF['type'] == 'RIGHT']
    right.drop(columns=['type'], inplace=True)

    internal = primDF[primDF['type'] == 'INTERNAL']
    internal.drop(columns=['type'], inplace=True)

    '''Merge all the data into one assay id per row'''
    mergedf = left.merge(right, on=['id'], suffixes=['_left', '_right']).merge(internal, on=['id'], suffixes=['_left', '_internal'])

    '''Calculate the heterodimers ie forward and reverse and probes deltaG and Tm in C'''
    #left and right
    mergedf['left_right_heterodimer_thermo'] = [primer.calcHeterodimer(i, j, mv_conc=50, dv_conc=4.7, dntp_conc=0.00095, dna_conc=200)
                                                for i, j in mergedf[['oligo_left', 'oligo_right']].values]
    mergedf['left_right_heterodimer_kcal/mol'] = [i.dg/1000 for i in mergedf['left_right_heterodimer_thermo'].values]
    mergedf['left_right_heterodimer_Tm_C'] = [i.tm for i in mergedf['left_right_heterodimer_thermo'].values]

    #left and probe
    mergedf['left_internal_heterodimer_thermo'] = [primer.calcHeterodimer(i, j, mv_conc=50, dv_conc=4.7, dntp_conc=0.00095, dna_conc=200)
                                                   for i, j in mergedf[['oligo_left', 'oligo']].values]
    mergedf['left_internal_heterodimer_kcal/mol'] = [i.dg / 1000 for i in mergedf['left_internal_heterodimer_thermo'].values]
    mergedf['left_internal_heterodimer_Tm_C'] = [i.tm for i in mergedf['left_internal_heterodimer_thermo'].values]

    #right and probe
    mergedf['internal_right_heterodimer_thermo'] = [primer.calcHeterodimer(i, j, mv_conc=50, dv_conc=4.7, dntp_conc=0.00095, dna_conc=200)
                                                    for i, j in mergedf[['oligo', 'oligo_right']].values]
    mergedf['right_internal_kcal/mol'] = [i.dg / 1000 for i in mergedf['internal_right_heterodimer_thermo'].values]
    mergedf['right_internal_heterodimer_Tm_C'] = [i.tm for i in mergedf['internal_right_heterodimer_thermo'].values]

    mergedf.drop(columns=['left_right_heterodimer_thermo', 'left_internal_heterodimer_thermo',
                          'internal_right_heterodimer_thermo'], inplace=True)

    mergedf['restric_enz_hit_all'] = mergedf[['restric_enz_hit_left', 'restric_enz_hit_right', 'restric_enz_hit']].any(axis=1)
    mergedf['3plus_repeats_all'] = mergedf[['3plus_repeats_left', '3plus_repeats_right', '3plus_repeats']].any(axis=1)

    print('Primer calculations completed')
    return mergedf


def pickPrimer(calcdf):
    '''
    Get the calculations and then perform the filtering steps and get final candidates of assays
    :param calcdf:
    :return: finaldf
    '''
    #remove rows where the enzymes can hit any sites
    finaldf = calcdf[calcdf['restric_enz_hit_all'] == False]

    #remove G/C clamps in the elimination of G or C sequences greater than 3 in the last 5 bases at 3prime end
    finaldf = finaldf[(finaldf['3_prime_GC_count_last_5_oligos_left'] < 4) &
                      (finaldf['3_prime_GC_count_last_5_oligos_right'] < 4) &
                      (finaldf['3_prime_GC_count_last_5_oligos'] < 4)]

    #remove assays with oligos more than 3 repeats
    finaldf = finaldf[finaldf['3plus_repeats_all'] == False]

    #Remove assays where oligos will make heterodimers in reaction temperature
    finaldf = finaldf[(finaldf['left_right_heterodimer_Tm_C'] < 51) &
                      (finaldf['left_internal_heterodimer_Tm_C'] < 51) &
                      (finaldf['right_internal_heterodimer_Tm_C'] < 51)]

    if finaldf.empty:
        print('ERROR: No data in dataframe')

    finaldf = finaldf.sort_values(by=['primer_len_left', 'primer_len_right', 'primer_len', 'Tm_calc_left',
                                      'Tm_calc_right', 'Tm_calc'])

    finaldf = finaldf[['id', 'oligo_left', 'oligo_right', 'oligo',
                       'Tm_calc_left', 'Tm_calc_right', 'Tm_calc',
                       'Homodimer_Tm_C_left', 'Homodimer_Tm_C_right', 'Homodimer_Tm_C']]
    return finaldf


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
    parser.add_argument('--num_of_assays', help='return number of candidate assays from pipeline', type=int)
    parser.add_argument('--mv_cations', help='Concentration of monovalent cations in mM', type=float)
    parser.add_argument('--dv_cations', help='concentration of divalent cations in mM', type=float)
    parser.add_argument('--dntp', help='concentration of dntps in mM [according to primer3 python docs]', type=float)
    parser.add_argument('--DNA', help='concentration of DNA in nM', type=float)
    parser.add_argument('--anneal_temp', help='annealing temp in C', type=float)
    parser.add_argument('--out_path', help='Path to save outputs')

    args = parser.parse_args()
    fasta = args.fastafile
    include_reg_start = args.included_region_start
    include_reg_len = args.included_region_length
    assay_num = args.num_of_assays
    outpath = args.out_path

    # creating SeqIO biopython object
    seq = SeqIO.parse(fasta, 'fasta')

    '''
    * This is code is designed to run bulk requests, so it can be used to design primers
      for either qPCR or sequencing applications.
    * Just provide fasta with either multiple genes or single genes. The program should
      be able to handle it.
    '''
    for record in seq:
        design = primerDesign(record, include_reg_start, include_reg_len, assay_num) #design primers
        resultdf = genDF(design) #transfer designs to a dataframe
        calcprimdf = addCalc(resultdf) #calculate hetero and homo dimers, and search for restriction sites
        rawfile = outpath + '/' + record.id + '_raw' + '.csv'
        calcprimdf.to_csv(rawfile, index=False)
        finalPrim = pickPrimer(calcprimdf) #pick primers
        finalfile = outpath + '/' + record.id + '_final' + '.csv'
        finalPrim.to_csv(finalfile, index=False) #Final output for the record


