from collections import Counter

# calculating k-mer frequencies in a given sequence
def count_kmers(sequence: str, k: int) -> Counter:
    """Calculate k-mer frequencies in a given sequence"""
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    return Counter(kmers)

# calculating the number of k-mer matches between two sequences
def basic_kmer_matches(seq1: str, seq2: str, k: int) -> int:
    """Calculate the number of k-mer matches between two sequences"""
    kmer_count1 = count_kmers(seq1, k)
    kmer_count2 = count_kmers(seq2, k)
    matches = 0
    for kmer in kmer_count1:
        if kmer in kmer_count2:
            # matches += min(kmer_count1[kmer], kmer_count2[kmer])
            matches += kmer_count1[kmer] * kmer_count2[kmer]
    return matches