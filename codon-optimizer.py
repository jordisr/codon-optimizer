'''
Generate DNA sequence from protein by sampling target codon frequencies.

Codon frequency talbe should be formatted as output from GCG CodonFrequency e.g.
https://downloads.yeastgenome.org/unpublished_data/codon/ysc.gene.cod
'''

import numpy as np
import sys

# amino acid three-letter conversion table
three_to_one  = {
    'Ala':'A',
    'Arg':'R',
    'Asn':'N',
    'Asp':'D',
    'Cys':'C',
    'Glu':'E',
    'Gln':'Q',
    'Gly':'G',
    'His':'H',
    'Ile':'I',
    'Leu':'L',
    'Lys':'K',
    'Met':'M',
    'Phe':'F',
    'Pro':'P',
    'Ser':'S',
    'Thr':'T',
    'Trp':'W',
    'Tyr':'Y',
    'Val':'V',
    'End':'End'
}

# program options -default cutoff of 10%
FREQUENCY_CUTOFF = 0.1

def sample_codon(aa, codon_table):
    return(np.random.choice(codon_table[aa][0],p=codon_table[aa][1]))


def load_codon_table(codon_table_path):
    # main data structure to hold dict of dicts
    codon_table = dict()

    # flag whether data section of file has begun
    READ_TABLE = False

    # parse codon frequency table into data structure
    with open(codon_table_path, 'r') as f:
        for line in f:
            cols = line.split()
            if len(cols) > 0:
                if cols[0] == 'AmAcid':
                    READ_TABLE = True
                elif READ_TABLE:
                    amacid = three_to_one[cols[0]]
                    codon = cols[1]
                    fraction = float(cols[4])

                    if amacid in codon_table:
                        codon_table[amacid][0].append(codon)
                        codon_table[amacid][1].append(fraction)
                    else:
                        codon_table[amacid] = [[codon], [fraction]]

    # remove stop codon
    del three_to_one['End']
                    
    # remove low frequency codons and rescale probabilities
    for amacid in codon_table:

        for i, frequency in enumerate(codon_table[amacid][1]):
            if frequency < FREQUENCY_CUTOFF:
                codon_table[amacid][0].pop(i)
                codon_table[amacid][1].pop(i)
        codon_table[amacid][1] = np.array(codon_table[amacid][1])
        codon_table[amacid][1] /= np.sum(codon_table[amacid][1])

    return(codon_table)
  
# read FASTA file of protein sequences and sample codon-optimzied DNA
# just generate something to test out
# TODO add FASTA reading code and generate DNA sequence for each protein

if __name__ == '__main__':

    codon_table = load_codon_table(sys.argv[1])
    sequence = np.random.choice(list(three_to_one.values()),size=100)

    print('Protein:')
    print(''.join(sequence))
    print('Codon-optimized:')
    print(''.join([sample_codon(aa, codon_table) for aa in sequence]))
