import numpy as np
from sklearn.preprocessing import normalize
from typing import Dict, List, Tuple
from disease_predictor.entities import Gene, Network
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime


class DiseasePredictorConfig:

    def __init__(self, network_filename: str,
                 network_separator: str,
                 gene_mapping_filename: str,
                 gene_mapping_separator: str,
                 combined_score_threshold: float = 0.0,
                 random_walk_steps: int = 10,
                 restart_probability: float = 0.2):
        self.network_filename = network_filename
        self.network_separator = network_separator
        self.gene_mapping_filename = gene_mapping_filename
        self.gene_mapping_separator = gene_mapping_separator
        self.combined_score_threshold = combined_score_threshold
        self.random_walk_steps = random_walk_steps
        self.restart_probability = restart_probability

    network_filename: str
    network_separator: str
    gene_mapping_filename: str
    gene_mapping_separator: str
    combined_score_threshold: float
    random_walk_steps: int
    restart_probability: float


class DiseasePredictor:
    config: DiseasePredictorConfig
    scores: List[Tuple[str, np.array]] = None

    def __init__(self, config: DiseasePredictorConfig, network: Network = None):
        self.config = config
        if network is not None:
            self.__network = network
        else:
            print("Building the network...")
            self._build_network()

    '''
    Given a list of seed genes and a list of candidate genes,
    it computes a random-walk-based score for each candidate gene.
    Such score is computed as the average between the values
    obtained by computing a random walk with restart starting from
    each seed gene.
    This method returns a list of tuples (gene,score), sorted by
    decreasing score.
    '''
    def rank_genes(self,
                   seed_genes: List[str],
                   candidate_genes: List[str] = None):
        print("Executing random walk...")
        scores = self._random_walk(seed_genes)
        print("Ranking candidate genes...")
        self.scores = self._rank_candidate_genes(scores, candidate_genes)
        return self.scores

    '''
    Returns the internal representation of the PPI network.
    '''
    def get_network(self):
        return self.__network

    '''
    Plots the scores on a strip plot.
    Genes can be grouped by category (e.g. to compare scores 
    of genes associated to different diseases).
    '''
    def plot_scores(self, categories: List[str]):
        if self.scores is None:
            print("No scores detected. rank_genes has to be called to compute them")
            return
        sns.set_palette("pastel")
        data = []
        for i, s in enumerate(self.scores):
            data.append([categories[i], s[1]])
        df = pd.DataFrame(data, columns=['Disease', 'Score'])
        sns.stripplot(data=df, x="Disease", y="Score", hue="Disease", jitter=0.1)
        plt.savefig(f"plot_{str(datetime.now().timestamp()).split('.')[0]}")

    def _build_network(self):
        network = Network()
        network.nodes = []
        protein_genes = {}

        with open(self.config.gene_mapping_filename, 'r') as f:
            for l in f.read().splitlines()[1:]:
                protein_info = l.split(self.config.gene_mapping_separator)
                protein_genes[protein_info[0]] = protein_info[1]

        # Given N = #genes
        # Build a NxN matrix for keeping track of scores between all genes
        adj_matrix = np.zeros((len(protein_genes.keys()), len(protein_genes.keys())))

        # This map will contain, for each gene in the network, the related
        # position in the adjacency matrix. I adopted this approach because
        # querying for the position of a gene is often performed, and by using
        # a map as an index I reduce this operation's complexity to O(1) each time.
        network.genes_idx = {}

        with open(self.config.network_filename, 'r') as f:
            pos = 0
            for l in f.read().splitlines()[1:]:
                edge = l.split(self.config.network_separator)
                score = float(edge[2])
                try:
                    src_gene = protein_genes[edge[0]]
                    dst_gene = protein_genes[edge[1]]
                except KeyError:
                    print(f"Skipping edge between {edge[0]} and {edge[1]}: one"
                          f"or both proteins not in nodes list")
                    continue
                if not (src_gene in network.genes_idx.keys()):
                    network.genes_idx[src_gene] = pos
                    network.nodes.append(Gene(gene_id=src_gene, protein_id=edge[0]))
                    pos = pos + 1
                if not (dst_gene in network.genes_idx.keys()):
                    network.genes_idx[dst_gene] = pos
                    network.nodes.append(Gene(gene_id=dst_gene, protein_id=edge[1]))
                    pos = pos + 1
                if score >= self.config.combined_score_threshold:
                    src_gene_pos = network.genes_idx[src_gene]
                    dst_gene_pos = network.genes_idx[dst_gene]

                    # Undirected network: M[i,j] = M[j,i]
                    adj_matrix[src_gene_pos, dst_gene_pos] = score
                    adj_matrix[dst_gene_pos, src_gene_pos] = score
        network.adj_matrix = adj_matrix
        # Normalize the matrix by row
        network.normalized_adj_matrix = normalize(adj_matrix, axis=1, norm='l1')
        self.__network = network
        return network

    def _random_walk(self, seed_genes: List[str]):
        if len(seed_genes) == 0:
            raise Exception("The seed genes list cannot be empty")

        gene_scores = {}
        for seed_gene in seed_genes:
            try:
                curr_gene_idx = self.__network.genes_idx[seed_gene]
            except KeyError:
                print(f"Seed gene {seed_gene} not present in network")
                continue

            adj_matrix_dim = self.__network.normalized_adj_matrix.shape[0]
            probabilities = np.zeros(adj_matrix_dim)
            # The random walk starts from the seed gene
            probabilities[curr_gene_idx] = 1
            start_probabilities = probabilities.copy()
            for i in range(self.config.random_walk_steps):
                # Random walk with restart: i-th iteration
                step_probabilities = self.__network.normalized_adj_matrix.dot(probabilities)
                probabilities = (1 - self.config.restart_probability) * step_probabilities + \
                                self.config.restart_probability * start_probabilities
            gene_scores[seed_gene] = probabilities.transpose()
        return gene_scores

    def _rank_candidate_genes(self, scores: Dict[str, np.array], candidate_genes: List[str] = None):
        # Remove duplicates, if any
        if candidate_genes is not None:
            candidate_genes = list(dict.fromkeys(candidate_genes))
        if candidate_genes is None or len(candidate_genes) == 0:
            candidate_genes = [g.gene_id for g in self.__network.nodes if g.gene_id not in list(scores.keys())]
        else:  # Remove seed genes (if any) from candidate genes
            candidate_genes = [cg for cg in candidate_genes if cg not in list(scores.keys())]
        candidate_genes_scores = {}
        score_arrays = [v for k, v in scores.items()]

        for cg in candidate_genes:
            try:
                cg_position = self.__network.genes_idx[cg]
            except KeyError:
                print(f"Candidate gene {cg} not present in network")
                continue
            cg_scores = np.zeros(len(list(scores.keys())))
            for idx, s in enumerate(score_arrays):
                cg_scores[idx] = s[cg_position]
            cg_avg_score = np.mean(cg_scores)
            candidate_genes_scores[cg] = cg_avg_score

        final_scores = [(k, v) for k, v in candidate_genes_scores.items()]
        final_scores.sort(key=lambda x: x[1], reverse=True)
        return final_scores
