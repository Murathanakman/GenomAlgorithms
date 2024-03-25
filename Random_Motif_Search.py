import random
import time

def create_random_motifs(dna, k):
    return [random.choice([d[i:i+k] for i in range(len(d) - k + 1)]) for d in dna]

def score_motifs(motifs, k):
    score = 0
    for i in range(k):
        motif_column = [motif[i] for motif in motifs]
        most_common = max(set(motif_column), key=motif_column.count)
        score += sum(1 for base in motif_column if base != most_common)
    return score

def profile_with_pseudocounts(motifs):
    profile = []
    k = len(motifs[0])
    for i in range(k):
        column = [motif[i] for motif in motifs]
        profile_column = {base: (column.count(base) + 1) / (len(column) + 4) for base in 'ACGT'}
        profile.append(profile_column)
    return profile

def find_consensus(motifs):
    profile = profile_with_pseudocounts(motifs)
    consensus = ''
    for column in profile:
        consensus += max(column, key=column.get)
    return consensus

def randomized_motif_search(dna, k):
    best_motifs = create_random_motifs(dna, k)
    best_score = score_motifs(best_motifs, k)
    score_history = [best_score]
    
    while True:
        start_time = time.time()

        i = random.randint(0, len(dna[0]) - k)
        motifs = [d[i:i+k] for d in dna]
        current_score = score_motifs(motifs, k)
        
        if current_score < best_score:
            best_motifs, best_score = motifs, current_score
            score_history = [best_score]  # Reset since improvement occurred
        else:
            score_history.append(current_score)
        
        if len(score_history) > 50 and all(s >= best_score for s in score_history[-50:]): # Stuck
            break

    end_time = time.time() 
    execution_time = end_time - start_time

    return best_motifs, best_score, execution_time

# DNA dizilerini dosyadan okuma
with open("Data/DNA_String_Motif.txt", "r") as file:
    dna = [line.strip() for line in file.readlines()]

# k değerleri için döngü
for k in [9, 10, 11]:
    best_motifs, best_score, execution_time = randomized_motif_search(dna, k)
    consensus = find_consensus(best_motifs)
    print(f"k= {k}")
    print("----------")
    print("Best Score: ", {best_score})
    print("Best Motifs: ", best_motifs)
    print("Consensus: ", consensus)
    print("Execution Time: ", execution_time)
    print("\n")
