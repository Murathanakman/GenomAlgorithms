import random
from itertools import product
import time

# File Reading
def read_dna_sequences(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file.readlines()]

# Hamming Distance Calculator
def hamming_distance(s1, s2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

# All k-mers
def all_kmers(k):
    return [''.join(p) for p in product('ACGT', repeat=k)]

# Median String
def median_string(dna, k):
    min_distance = float('inf')
    execution_time = time.time()
    for kmer in all_kmers(k):
        distance = sum(min(hamming_distance(kmer, seq[i:i+k]) for i in range(len(seq) - k + 1)) for seq in dna)
        if distance < min_distance:
            min_distance = distance
            median_kmer = kmer
    
    execution_time = time.time() - execution_time
    return median_kmer, execution_time

def main():
    file_path = 'Data/DNA_String_Motif.txt'
    k_values = [9, 10, 11]
    dna = read_dna_sequences(file_path)
    
    for k in k_values:
        median, execution_time = median_string(dna, k)
        print(f'Median String for k={k}: {median}')
        print(f'Execution Time: {execution_time:.10f} seconds\n')

if __name__ == '__main__':
    main()