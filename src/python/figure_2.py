import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import Counter
from typing import List
from tqdm import tqdm
import numpy as np

# 1. loading DNA sequences from a FASTA file
def load_sequences(file_path: str) -> List[str]:
    """Loading DNA sequences from a FASTA file"""
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    return sequences

# 2. calculating k-mer frequencies in a given sequence
def count_kmers(sequence: str, k: int) -> Counter:
    """Calculate k-mer frequencies in a given sequence"""
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    return Counter(kmers)

# 3. calculating the number of k-mer matches between two sequences
def kmer_matches(seq1: str, seq2: str, k: int) -> int:
    """Calculate the number of k-mer matches between two sequences"""
    kmer_count1 = count_kmers(seq1, k)
    kmer_count2 = count_kmers(seq2, k)
    matches = 0
    for kmer in kmer_count1:
        if kmer in kmer_count2:
            matches += min(kmer_count1[kmer], kmer_count2[kmer])
    return matches

# 4. calculating F(k) for different k values
def compute_F_k(seq1: str, seq2: str, k_values: List[int]) -> List[float]:
    """Calculate F(k) for different k values"""
    F_k = []
    L_1 = len(seq1)
    L_2 = len(seq2)
    L_min = min(L_1, L_2)  # minimum sequence length
    q = 0.25  # assume equal frequency of each nucleotide, background match probability

    for k in tqdm(k_values, desc="Calculating F(k) for different k values"):
        matches = kmer_matches(seq1, seq2, k)

        # calculate the number of background matches
        background_matches = 2 * L_1 * L_2 * (q ** k)

        # generate F(k) value
        if matches > background_matches:
            F_k_value = np.log(matches - background_matches)
            F_k.append(F_k_value)
        else:
            F_k.append(0)  # if matches <= background_matches, F(k) = 0

    return F_k

# 5. plotting F(k) vs k
def plot_F_k(k_values: List[int], F_k: List[int]):
    """Plot F(k) vs k"""

    plt.figure(figsize=(10, 6))
    plt.plot(k_values[:len(F_k)], F_k, marker='o', linestyle='-', color='b')
    plt.xlabel('Word Length (k)')
    plt.ylabel('Number of Matches F(k)')
    plt.title('Relationship between F(k) and Word Length (k)')
    plt.grid(True)
    plt.savefig('F_k_vs_k.png') 
    plt.show()

def main():
    # take fasta file as input
    fasta_file = "/home/chenye/project/slope-spam/chenye-implement/dataset/test.fasta"
    sequences = load_sequences(fasta_file)

    # make sure there are at least two sequences
    if len(sequences) < 2:
        print("Error: you need 2 seqs to calculate F(k)")
        return

    seq1 = sequences[0]
    seq2 = sequences[1]

    # define k values range
    k_values = list(range(0, 25))  # k_min = 10, k_max = 24 in paper

    # calculate F(k) for different k values
    F_k = compute_F_k(seq1, seq2, k_values)

    # plot F(k) vs k
    plot_F_k(k_values, F_k)

if __name__ == "__main__":
    main()

