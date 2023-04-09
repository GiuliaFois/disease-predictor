import numpy as np
from typing import Dict, List


class Gene:
    gene_id: str
    protein_id: str

    def __init__(self, gene_id: str, protein_id: str = None):
        self.gene_id = gene_id
        self.protein_id = protein_id


class Network:
    nodes: List[Gene]
    adj_matrix: np.array
    normalized_adj_matrix: np.array
    genes_idx = Dict[str, int]
