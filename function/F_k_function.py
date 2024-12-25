import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(project_root)

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

from utils.file_system import load_sequences
from methods.kmer_matches import *

class FKFuntion():
    """
    calculating F(k) for different k values

    """
    def __init__(self, args):
        self.k_mers_method_valid_list = args.get("k_mers_method_valid_list", ["basic_kmer_matches"])
        self.show_k_range_min = args.get("k_range_min", 1)
        self.show_k_range_max = args.get("k_range_max", 25)
        self.log_epsilon = args.get("log_epsilon", 1e-10)
        self.data_path = args.get("data_path", None)
        self.seq_list = load_sequences(self.data_path)
        self.figure_root_path = args.get("figure_root_path", None)


    def compute_F_k(self, used_k_mers_method):
        F_k = []
        seq_1 = self.seq_list[0]
        seq_2 = self.seq_list[1]
        L_1 = len(seq_1)
        L_2 = len(seq_2)

        q = 0.25  # assume equal frequency of each nucleotide, background match probability

        k_values = list(range(self.show_k_range_min, self.show_k_range_max))
        self.k_values = k_values

        for k in tqdm(k_values, desc="Calculating F(k) for different k values"):
            if used_k_mers_method in self.k_mers_method_valid_list:
                matches = eval(f"{used_k_mers_method}(seq_1, seq_2, k)")
            else:
                raise ValueError(f"Invalid alignment-free method: {used_k_mers_method}")

            # calculate the number of background matches
            background_matches = (L_1 - k + 1) * (L_2 - k + 1) * (q ** k)

            # generate F(k) value
            # if matches > background_matches:
            #     F_k_value = np.log(matches - background_matches)
            #     F_k.append(F_k_value)
            # else:
            #     F_k.append(0)  # if matches <= background_matches, F(k) = 0
            F_k_value = np.log(max(matches - background_matches, self.log_epsilon))
            F_k.append(F_k_value)

        return F_k
    
    def plot_F_k(self, F_k):
        """
        Plot F(k) vs k
        """
        plt.figure(figsize=(10, 6))
        plt.plot(self.k_values[:len(F_k)], F_k, marker='o', linestyle='-', color='b')
        plt.xlabel('Word Length (k)')
        plt.ylabel('Number of Matches F(k)')
        plt.title('Relationship between F(k) and Word Length (k)')
        plt.grid(True)
        plt.savefig(self.figure_root_path + 'F_k_vs_k.png') 
        plt.show()
