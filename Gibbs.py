import random
import time

# File Reading
def read_dna_from_file(file_path):
    with open(file_path, 'r') as file:
        return [line.strip() for line in file.readlines()]

# Randomly select k-mers from each sequence
def random_motifs(dna, k):
    return [random.choice([seq[i:i+k] for i in range(len(seq) - k + 1)]) for seq in dna]

# Create profile matrix
def profile_with_pseudocounts(motifs):
    profile = []
    for i in range(len(motifs[0])):
        col = ''.join([motif[i] for motif in motifs])
        profile.append({base: (col.count(base) + 1) / (len(col) + 4) for base in 'ACGT'})
    return profile

# Randomly select a k-mer from a sequence based on the profile matrix
def weighted_random_kmer(seq, k, profile):
    weights = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        weight = 1
        for j, base in enumerate(kmer):
            weight *= profile[j][base]
        weights.append(weight)
    total_weight = sum(weights)
    r = random.uniform(0, total_weight)
    for i, weight in enumerate(weights):
        r -= weight
        if r <= 0:
            return seq[i:i+k]

# Calculate the score of the motifs
def score(motifs):
    score = 0
    for i in range(len(motifs[0])):
        col = [motif[i] for motif in motifs]
        max_freq = max(col.count(base) for base in 'ACGT')
        score += len(motifs) - max_freq
    return score

# Find the consensus string
def find_consensus(motifs):
    consensus = ''
    for i in range(len(motifs[0])):
        col = [motif[i] for motif in motifs]
        consensus += max('ACGT', key=col.count)
    return consensus

# Gibbs Sampler
def gibbs_sampler(dna, k):
    motifs = random_motifs(dna, k)
    best_motifs = motifs[:]
    best_score = score(motifs)
    scores = [best_score]  # Keep track of the scores

    while True:
        start_time = time.time()

        i = random.randint(0, len(dna) - 1)
        motifs_except_i = motifs[:i] + motifs[i+1:]
        profile = profile_with_pseudocounts(motifs_except_i)
        motifs[i] = weighted_random_kmer(dna[i], k, profile)
        current_score = score(motifs)

        if current_score < best_score:
            best_motifs = motifs[:]
            best_score = current_score
            scores.append(best_score)  # Update
        else:
            scores.append(current_score) # Keep
         
        if len(scores) > 50 and scores[-1] == scores[-51]:
            break  # Stuck

    end_time = time.time() 
    execution_time = end_time - start_time  

    return best_motifs, best_score, execution_time

# Path to the file
file_path = "Data/DNA_String_Motif.txt"
dna = read_dna_from_file(file_path)

# k value and Iteration count
k_values = [9, 10, 11]

# Run the Gibbs Sampler
for k in k_values:
    best_motifs, best_score, execution_time = gibbs_sampler(dna, k)
    consensus = find_consensus(best_motifs)
    print(f"k= {k}")
    print("----------")
    print("Best Score: ", {best_score})
    print("Best Motifs: ", best_motifs)
    print("Consensus: ", consensus)
    print("Execution Time: {:.10f} seconds".format(execution_time))
    print("\n")
