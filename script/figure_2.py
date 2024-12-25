import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.abspath(os.path.join(current_dir, os.pardir))
sys.path.append(project_root)
from utils.file_system import load_sequences

import yaml
from easydict import EasyDict

from function.F_k_function import FKFuntion

def main():
    # take fasta file as input
    config_path = "../config/run.yaml"
    with open(config_path) as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    args = EasyDict(config)

    f_k_function_object = FKFuntion(args)

    F_k = f_k_function_object.compute_F_k("basic_kmer_matches")

    # plot F(k) vs k
    f_k_function_object.plot_F_k(F_k)

if __name__ == "__main__":
    main()

