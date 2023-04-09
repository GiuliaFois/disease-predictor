from disease_predictor import rank_genes

# Read genes from local files. Genes belong to two categories:
# they can be associated to diabetes or depression, with possible overlap
diabetes_genes = []
depression_genes = []
with open('resources/genes_diabetes.txt') as f:
    for l in f.readlines()[1:]:
        diabetes_genes.append(l.split('\t')[5])

with open('resources/genes_depression.txt') as f:
    for l in f.readlines()[1:]:
        depression_genes.append(l.split('\t')[5])

# Remove the overlap
diabetes_genes, depression_genes = [i for i in diabetes_genes if i not in depression_genes], \
                                        [j for j in depression_genes if j not in diabetes_genes]

# Build the DiseasePredictor object
config = rank_genes.DiseasePredictorConfig(
network_filename="resources/9606.protein.physical.links.v11.5.txt",
                                    network_separator=" ",
                                    gene_mapping_filename="resources/9606.protein.info.v11.5.txt",
                                    gene_mapping_separator="\t",
                                    combined_score_threshold=50,
                                    random_walk_steps=5)
dp = rank_genes.DiseasePredictor(config=config)

# Rank diabetes/depression candidate genes based on diabetes seed genes
seed_genes = diabetes_genes[:100]
candidate_genes = diabetes_genes[100:201] + depression_genes[:101]
scores = dp.rank_genes(seed_genes=seed_genes, candidate_genes=candidate_genes)
count = 0

# Build a plot to show if scores are different based on the disease
positions = []
genes = []
categories = []
for i, s in enumerate(scores):
    if s[0] in diabetes_genes:
        genes.append(s[0])
        count = count + 1
        categories.append('diabetes')
    elif s[0] in depression_genes:
        categories.append('depression')

dp.plot_scores(categories=categories)
