from disease_predictor.rank_genes import DiseasePredictorConfig, DiseasePredictor
from disease_predictor.entities import Network, Gene
import numpy as np
import pytest


@pytest.fixture
def dp_config():
    return DiseasePredictorConfig(network_filename="./test/resources/protein_protein_interactions.txt",
                                  network_separator=',',
                                  gene_mapping_filename="./test/resources/protein_genes.txt",
                                  gene_mapping_separator=',')


def test_build_network(dp_config):
    dp = DiseasePredictor(config=dp_config)
    dp._build_network()
    network = dp.get_network()
    assert len(network.nodes) == 4
    assert network.adj_matrix.shape == (4, 4)
    assert network.normalized_adj_matrix.shape == (4, 4)
    expected_matrix = np.array([[0, 5, 4.67, 0],
                                [5, 0, 0, 1],
                                [4.67, 0, 0, 0],
                                [0, 1, 0, 0]])
    for i in range(network.adj_matrix.shape[0]):
        matrix_row = network.adj_matrix[i]
        normalized_matrix_row = network.normalized_adj_matrix[i]
        assert np.allclose(matrix_row, expected_matrix[i])
        expected_normalized_row = expected_matrix[i] / np.linalg.norm(expected_matrix[i], ord=1)
        assert np.allclose(normalized_matrix_row, expected_normalized_row)


@pytest.mark.parametrize(
    "normalized_adj_matrix,expected_scores",
    [
        pytest.param(np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]]), np.array([1, 0, 0])),
        pytest.param(np.array([[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]), np.array([0.375, 0.3125, 0.3125]))
    ]
)
def test_random_walk(dp_config, normalized_adj_matrix, expected_scores):
    dp_config.random_walk_steps = 4
    # To prevent the random element from obtaining a testable result
    dp_config.restart_probability = 0
    network = Network()
    network.nodes = [Gene(gene_id='gene1'), Gene(gene_id='gene2'), Gene(gene_id='gene3')]
    network.genes_idx = {"gene1": 0, "gene2": 1, "gene3": 2}
    network.normalized_adj_matrix = normalized_adj_matrix
    dp = DiseasePredictor(config=dp_config, network=network)
    result = dp._random_walk(seed_genes=['gene1'])
    assert np.allclose(result["gene1"], expected_scores)


def test_random_walk_skips_unknown_seed(dp_config, capsys):
    network = Network()
    network.nodes = [Gene(gene_id='gene1'), Gene(gene_id='gene2'), Gene(gene_id='gene3')]
    network.genes_idx = {"gene1": 0, "gene2": 1, "gene3": 2}
    network.normalized_adj_matrix = np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]]), np.array([1, 0, 0])
    dp = DiseasePredictor(config=dp_config, network=network)
    dp._random_walk(seed_genes=['gene4'])
    captured = capsys.readouterr()
    assert captured.out == "Seed gene gene4 not present in network\n"


def test_random_walk_throws_for_empty_seed(dp_config):
    network = Network()
    dp = DiseasePredictor(config=dp_config, network=network)
    expected_err_msg = "The seed genes list cannot be empty"
    with pytest.raises(Exception, match=expected_err_msg):
        dp._random_walk(seed_genes=[])


def test_rank_candidate_genes(dp_config):
    network = Network()
    network.nodes = [Gene(gene_id=f"gene{i}") for i in range(1, 6)]
    network.genes_idx = {f"gene{i + 1}": i for i in range(5)}
    dp = DiseasePredictor(config=dp_config, network=network)
    scores = {'gene1': np.array([0.1, 0.1, 0.4, 0.2, 0.2]),
              'gene2': np.array([0.025, 0.1, 0.025, 0.15, 0.7])}
    res = dp._rank_candidate_genes(scores=scores)
    assert [(r[0], float("{:.4f}".format(r[1]))) for r in res] == [('gene5', 0.4500), ('gene3', 0.2125), ('gene4', 0.1750)]


def test_rank_candidate_genes_skips_unknown_candidate(dp_config, capsys):
    network = Network()
    network.nodes = [Gene(gene_id=f"gene{i}") for i in range(1, 6)]
    network.genes_idx = {f"gene{i + 1}": i for i in range(5)}
    dp = DiseasePredictor(config=dp_config, network=network)
    scores = {'gene1': np.array([0.1, 0.1, 0.4, 0.2, 0.2]),
              'gene2': np.array([0.025, 0.1, 0.025, 0.15, 0.7])}
    dp._rank_candidate_genes(scores=scores, candidate_genes=['gene10'])
    captured = capsys.readouterr()
    assert captured.out == "Candidate gene gene10 not present in network\n"