import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


class GSEA(object):
    def __init__(self):        
        self.expressions = pd.DataFrame()
        self.sample = {}
        self.genesets = {}
        
    def load_data(self, expfile, sampfile, keggfile):
        #function to load data from files into GSEA object
        #inputs: expfile, sampfile, genesets
        #stores data in respective data structure
        self.expressions = pd.read_csv(expfile, sep ='\t')
        sample_pd = pd.read_csv(sampfile, sep='\t', header = None)
        self.sample = dict(zip(sample_pd[0], sample_pd[1]))
        with open(keggfile, "r") as f:
            for line in f:
                columns = line.strip().split('\t')
                self.genesets[columns[0]] = set(columns[2:])

    def get_gene_rank_order(self):
        #function to rank genes by their logFC (difference between patient and control) and order from greatest to smallest
        #input: 
        #returns: a list of all genes (strings) ranked by logFC in desc order
        sample_names = np.array(list(self.sample.keys()))
        sample_labels = np.array(list(self.sample.values()))

        healthy_samples = sample_names[sample_labels ==0]
        patient_samples = sample_names[sample_labels ==1]
        healthy_exp_mean = self.expressions[healthy_samples].mean(axis= 1)
        patient_exp_mean = self.expressions[patient_samples].mean(axis= 1)

        logFC = patient_exp_mean - healthy_exp_mean
        #onnect logFC values of genes back to gene name in self.expression
        logFC.index = self.expressions['SYMBOL']
        logFC = logFC.sort_values(ascending = False)
        ranked_genes = list(logFC.index)
        return ranked_genes

    def get_enrichment_score(self, geneset):
        ranked_genes = self.get_gene_rank_order()
        return self.calculate_es(geneset, ranked_genes)

    def calculate_es(self, geneset, ranked_genes):
        pathway_genes = self.genesets[geneset]
        ranked_geneset = set(ranked_genes) & pathway_genes

        num_ranked_genes = len(ranked_genes)
        num_ranked_geneset = len(ranked_geneset)
        num_not_ranked_geneset = num_ranked_genes - num_ranked_geneset

        if num_ranked_geneset == 0 or num_not_ranked_geneset == 0:
            return 0.0

        up_score = np.sqrt((num_not_ranked_geneset) / num_ranked_geneset)
        down_score = -np.sqrt(num_ranked_geneset / num_not_ranked_geneset)

        ranked_genes_array = np.array(ranked_genes)
        geneset_in_ranked = np.isin(
            ranked_genes_array,
            list(ranked_geneset)
        )

        score_array = np.where(geneset_in_ranked, up_score, down_score)

        running_score = np.cumsum(score_array)

        max_es = round(float(np.max(running_score)), 2)
        return max_es
        
    def get_sig_sets(self, p):
        #get actual enrichment scores
        ranked_genes = self.get_gene_rank_order()

        actual_scores = {geneset: self.calculate_es(geneset, ranked_genes) for geneset in self.genesets}

        # count how many permuted scores are >= actual scores
        permuted_counts = {geneset: 0 for geneset in self.genesets}

        #save original sample labels
        original_sample = self.sample.copy()

        # permute sample labels 100 times
        sample_names = np.array(list(self.sample.keys()))
        original_labels = np.array(list(self.sample.values()))

        for _ in range(100):
            permuted_labels = np.random.permutation(original_labels)

            self.sample = dict(zip(sample_names, permuted_labels))

            # calculate ranking once for this permutation
            permuted_ranked_genes = self.get_gene_rank_order()

            # calculate ES for each gene set using this ranking
            for geneset in self.genesets:
                permuted_es = self.calculate_es(geneset, permuted_ranked_genes)

                if permuted_es >= actual_scores[geneset]:
                    permuted_counts[geneset] += 1

        #restore original sample assignments
        self.sample = original_sample

        #calculate p-values
        p_values = {geneset: permuted_counts[geneset] / 100 for geneset in self.genesets}

        # Bonferroni correction
        num_tests = len(self.genesets)
        corrected_threshold = p / num_tests

        # Return significant gene sets
        significant_sets = [geneset for geneset in self.genesets
        if p_values[geneset] <= corrected_threshold]

        return significant_sets

def main():
    if len(sys.argv) != 4:
        return

    expfile = sys.argv[1]
    sampfile = sys.argv[2]
    genesets = sys.argv[3]

    gsea = GSEA()
    gsea.load_data(expfile, sampfile, genesets)

    all_genes = []

    for geneset in gsea.genesets.values():
        all_genes.extend(geneset)
    genes, counts = np.unique(all_genes, return_counts=True)

    max_count = np.max(counts)
    most_common_genes = genes[counts == max_count]

    print("Gene(s):", most_common_genes)
    print("Number of occurrences:", max_count)

    #plt.hist(counts, bins=np.arange(1, counts.max() + 2) - 0.5)
    #plt.xlabel("Number of occurrences")
    #plt.ylabel("Number of genes")
    #plt.show()

    #score = gsea.get_enrichment_score("pathway_name")
    #print(score)

    # Calculate enrichment scores for all KEGG gene sets
    #enrichment_scores = {}

    #for geneset in gsea.genesets:
        #enrichment_scores[geneset] = gsea.get_enrichment_score(geneset)

    # Sort from highest ES to lowest ES
    #sorted_scores = sorted(
        #enrichment_scores.items(),
        #key=lambda x: x[1],
        #reverse=True
    #)

    #print(sorted_scores[0])

    # Write tab-delimited debugging file
    #with open("kegg_enrichment_scores.txt", "w") as f:
        #for geneset, score in sorted_scores:
            #f.write(f"{geneset}\t{score:.2f}\n")

if __name__ == "__main__":
    main()