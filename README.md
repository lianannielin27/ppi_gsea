# Protein Interaction and Gene Set Enrichment Analyses

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

Both scripts accept paths to your own input files; filenames do not need to match the examples.

### Predict similar genes

```bash
python ppi.py <disease_gene_file> <interaction_network_file>
```

Replace the placeholders with your file paths:

```bash
python ppi.py data/my_disease_genes.txt data/my_interactions.txt
```

The disease gene file must contain one gene symbol per line. The interaction file must contain whitespace-separated gene pairs; a third interaction-weight column is ignored.

Set the cosine-distance threshold in `main()` when calling `get_close_genes()`.

**When using a different interaction network, remove or rename `node2vec_pretrained` before running the script** so embeddings are trained on your new network.

### Analyze gene set enrichment

```bash
python gsea.py <expression_file> <sample_assignment_file> <gene_set_file>
```

For example:

```bash
python gsea.py data/my_expression.txt data/my_samples.txt data/my_pathways.gmt
```

The expression file must be tab-delimited, with a `SYMBOL` column containing gene names and one column per sample. Expression values must already be log2-normalized.

The sample assignment file must have two tab-delimited columns without a header: sample ID and group label. The current code treats `0` as control and `1` as patient. Sample IDs must match the expression column names.

The gene set file must use GMT format: one pathway per line, with the pathway name, description, and gene symbols separated by tabs.

Gene identifiers must be consistent across the input files. Quote any file paths containing spaces. Results printed or saved depend on the analysis configured in each script's `main()`.

## Outputs

### PPI gene prediction

- `calculate_embedding()` returns gene names and their corresponding 64-dimensional embedding vectors.
- `get_close_genes(nodes, embeddings, threshold)` returns a set containing the known disease genes and predicted similar genes at the chosen cosine-distance threshold.
- Training saves a `node2vec_pretrained` model cache for reuse.

### Gene set enrichment

- `get_gene_rank_order()` returns gene names ranked from highest to lowest patient-versus-control log2 fold-change.
- `get_enrichment_score(geneset)` returns the enrichment score for a gene set named in your input file, rounded to two decimal places.
- `get_sig_sets(p)` returns the names of gene sets significant at the specified Bonferroni-corrected significance level, or an empty list if none qualify.

Printing results or exporting them to files is controlled by each script's `main()` function. Output filenames, prediction thresholds, and gene set names can be changed there to suit your analysis.

## Reproducibility and limitations

- Node2Vec uses 64 dimensions, walk length 30, 100 walks per node, one worker, and seed 42.
- Word2Vec training uses `window=3`, `min_count=1`, and `batch_words=4`.
- If `node2vec_pretrained` exists, it is loaded instead of retraining. Remove or rename it when changing the network or training settings.
- Permutation results can vary between runs. With only 100 permutations, empirical p-values have limited resolution.
- Results depend on network coverage, predefined gene sets, sample labels, and expression quality.

## Data

The example expression data come from GEO study **GSE25628** (Crispi et al., 2013). The KEGG collection was supplied with the course.

Data files are not bundled by default. Obtain them from authorized sources and check redistribution permissions before sharing them.
