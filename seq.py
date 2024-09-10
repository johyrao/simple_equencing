from Bio import SeqIO # Biopython for fast parsing of FASTA data
import pandas as pd # pandas for dataframes to organize outputs

def complement(seq):
    """
    Takes in a sequence and output the complement of that sequence. Note, this 
    does not reverse the sequence so the output is not the reverse complement.

    seq: String. The sequence to get the complement of.
    Return: String. The complement of the input sequence.
    
    """
    comp = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    comp_seq = "".join(comp.get(base, base) for base in seq)
    return comp_seq

# Load in the text file with the mutation targets
mutations = pd.read_csv('somatic_mutation_data.txt', sep='\t')

sequences = [] # Store the target sequences

for ind in range(len(mutations)):
    # Load information from mutation data table
    chr = 'chr' + mutations['chromosome'][ind] # Chromosome number of target
    loc = mutations['chromosome_start'][ind] - 1 # Location of target
    sense = mutations['chromosome_strand'][ind].strip() # Strand information
    ref = mutations['reference_genome_allele'][ind] # Expected reference base
    
    for sequence in SeqIO.parse(open('hg19.fa'), 'fasta'):
        if sequence.id == chr:
            seq = str(sequence.seq).upper()
            break
    
    # Check if the base at location is the same as the expected reference
    if seq[loc] != ref and sense == '+':
        raise Exception("Error in parsing sequence")

    # Get the 200bp sequence
    snp_seq = seq[loc-100:loc+100]

    # Check to see if the complement sequence needs to be calculated and append
    if sense == '+':
        sequences.append(snp_seq)
    else:
        sequences.append(complement(snp_seq))

# Get the counts of A, C, G, T in each extracted sequences
count_A = [sect.count('A') for sect in sequences]
count_C = [sect.count('C') for sect in sequences]
count_G = [sect.count('G') for sect in sequences]
count_T = [sect.count('T') for sect in sequences]

# Sanity check to make sure all the values add back up to 200 so that there will
#  not be any double or missing count
check_sum = [sum(x) for x in zip(count_A, count_C, count_G, count_T)]
for count, total in enumerate(check_sum):
    if total != 200:
        raise Exception(f"Mutation {count} individual single base total is "
                        f"not 200.")

# Get the 2-gram counts
count_AA = [sect.count('AA') for sect in sequences]
count_AC = [sect.count('AC') for sect in sequences]
count_AG = [sect.count('AG') for sect in sequences]
count_AT = [sect.count('AT') for sect in sequences]

count_CA = [sect.count('CA') for sect in sequences]
count_CC = [sect.count('CC') for sect in sequences]
count_CG = [sect.count('CG') for sect in sequences]
count_CT = [sect.count('CT') for sect in sequences]

count_GA = [sect.count('GA') for sect in sequences]
count_GC = [sect.count('GC') for sect in sequences]
count_GG = [sect.count('GG') for sect in sequences]
count_GT = [sect.count('GT') for sect in sequences]

count_TA = [sect.count('TA') for sect in sequences]
count_TC = [sect.count('TC') for sect in sequences]
count_TG = [sect.count('TG') for sect in sequences]
count_TT = [sect.count('TT') for sect in sequences]

# Put all the results back into data table
mutations['Sequence'] = sequences
mutations['A count'] = count_A
mutations['C count'] = count_C
mutations['G count'] = count_G
mutations['T count'] = count_T
mutations['AA count'] = count_AA
mutations['AC count'] = count_AC
mutations['AG count'] = count_AG
mutations['AT count'] = count_AT
mutations['CA count'] = count_CA
mutations['CC count'] = count_CC
mutations['CG count'] = count_CG
mutations['CT count'] = count_CT
mutations['GA count'] = count_GA
mutations['GC count'] = count_GC
mutations['GG count'] = count_GG
mutations['GT count'] = count_GT
mutations['TA count'] = count_TA
mutations['TC count'] = count_TC
mutations['TG count'] = count_TG
mutations['TT count'] = count_TT

# Display the data table
mutations.head(5)