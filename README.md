# Protein Interaction Analysis and Gene Set Enrichment

This repo contains python tools for predicting disease-associated genes from protein interaction networks and evaluating gene sets using gene set enrichment analysis (GSEA). Developed as a course project using endometriosis expression data.

## Methods

### Network-based gene prediction

`ppi.py` builds an undirected protein interaction graph and uses Node2Vec to represent each gene as a 64-dimensional vector. Genes within a chosen cosine-distance threshold of any known disease gene are selected. The returned set also includes the original disease genes.

Embeddings capture network context: a direct interaction does not necessarily imply a small embedding distance.

### Gene set enrichment analysis

`gsea.py` ranks genes by mean patient expression minus mean control expression. Because the input values are already log2-normalized, this difference represents log2 fold-change.

Enrichment is calculated using a Brownian-bridge running sum. Significance is estimated using 100 sample-label permutations, followed by Bonferroni correction.

Together, these methods can prioritize candidate genes and pathways for follow-up. They identify associations, not causal mechanisms.

## Installation

Use Python 3.12:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install numpy pandas scipy networkx node2vec gensim scikit-learn "setuptools<81"
```

The setuptools constraint supports older Node2Vec versions that depend on `pkg_resources`. Matplotlib is optional for plotting.

## Input formats

| File | Format |
| --- | --- |
| Disease genes | One gene symbol per line |
| Interaction network | Whitespace-delimited interacting gene pairs; the third weight column is ignored |
| Expression data | Tab-delimited table with a `SYMBOL` column and one column per sample |
| Sample assignments | Two tab-delimited columns without a header: sample ID and group label |
| Gene sets | GMT format: pathway name, description, and gene symbols separated by tabs |

The current implementation treats `0` as control and `1` as patient. Sample IDs must match the expression column names.

## Usage

The scripts accept:

```bash
python ppi.py disease_gene_list.txt interaction_network_filtered.txt
python gsea.py GSE25628_filtered_expression.txt GSE25628_samples.txt c2.cp.kegg.v6.2.symbols.filtered.gmt
```

Outputs depend on the analysis configured in each script's `main()`. To run a complete analysis through the class interface:

```python
from ppi import PPI
from gsea import GSEA

# Predict genes using a cosine-distance threshold of 0.2.
ppi = PPI()
ppi.load_data("disease_gene_list.txt", "interaction_network_filtered.txt")
nodes, embeddings = ppi.calculate_embedding()
predicted = ppi.get_close_genes(nodes, embeddings, 0.2)

# Save the prediction as a GMT-style gene set.
with open("predicted_genes.txt", "w") as f:
    f.write(
        "predicted_genes\tdescription\t"
        + "\t".join(sorted(predicted))
        + "\n"
    )

# Calculate enrichment of the predicted gene set.
gsea = GSEA()
gsea.load_data(
    "GSE25628_filtered_expression.txt",
    "GSE25628_samples.txt",
    "predicted_genes.txt",
)
print("Predicted genes:", len(predicted))
print("Enrichment score:", gsea.get_enrichment_score("predicted_genes"))

# Identify significant KEGG pathways.
gsea.load_data(
    "GSE25628_filtered_expression.txt",
    "GSE25628_samples.txt",
    "c2.cp.kegg.v6.2.symbols.filtered.gmt",
)
print("Significant pathways:", gsea.get_sig_sets(0.05))
```

For a separate KEGG analysis, create a fresh `GSEA()` object before loading the KEGG file to avoid retaining previously loaded gene sets.

## Reproducibility and limitations

- Node2Vec uses 64 dimensions, walk length 30, 100 walks per node, one worker, and seed 42.
- Word2Vec training uses `window=3`, `min_count=1`, and `batch_words=4`.
- If `node2vec_pretrained` exists, it is loaded instead of retraining. Remove or rename it when changing the network or training settings.
- Permutation results can vary between runs. With only 100 permutations, empirical p-values have limited resolution.
- Results depend on network coverage, predefined gene sets, sample labels, and expression quality.

## Data

The example expression data come from GEO study **GSE25628** (Crispi et al., 2013). The KEGG collection was supplied with the course.

Data files are not bundled by default. Obtain them from authorized sources and check redistribution permissions before sharing them.
