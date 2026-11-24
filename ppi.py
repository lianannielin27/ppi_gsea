import sys
import os
import numpy as np
from scipy.spatial.distance import cosine
import networkx as nx
from sklearn.metrics.pairwise import cosine_distances
from node2vec import Node2Vec
from gensim.models import Word2Vec

class PPI(object):
    def __init__(self):
        #what info does class PPI need?
        #initialise known disease genes as emtpy set
        self.disease_genes = set()
        #initialise graph as empty networkx graph for gene interactions
        self.graph = nx.Graph()

    def load_data(self, diseaseGeneFile, interactionNetworkFile):
        # read the dsiease genes and interaction network file
        # inputs: 
        #   diseasegenefile, path to a file containing oen disease gene per line
        #   interactionnetworkfile, path to a file containing gene interactions
        # returns none
        # store disease genes as set
        with open(diseaseGeneFile, 'r') as f:
            self.disease_genes = set(f.read().splitlines())
        # build/store interaction graph in self.graph
        self.graph = nx.read_edgelist(interactionNetworkFile, data = False)

    def calculate_embedding(self):
        #function to calculate node2vec embeddings for genes in interaction network
        #inputs: none
        #returns: 
        #    nodes: list of gene names from graph
        #    vector_embeddings_list: list of embedding vectors
        #load pretrained model if exists
        if os.path.exists("node2vec_pretrained"):
            model = Word2Vec.load("node2vec_pretrained")
            print("loaded pretrained model successfuly")
        # train node2vec
        else:
            #initialize parametres
            node2vec = Node2Vec(
                self.graph, 
                dimensions= 64,     #embedding dimesions
                walk_length= 30,    # number of nodes in each walk
                num_walks= 100,     # number of walks per node         
                workers=1,          # non-parellel execution
                seed = 42
            )
            #train model to get embeddings
            model = node2vec.fit(window=3, min_count=1, batch_words=4)
            model.save("node2vec_pretrained")

        #extract nodes and vector embeddings
        nodes = list(model.wv.index_to_key)
        vector_embeddings_list = model.wv.vectors.tolist()
        #return both
        return nodes, vector_embeddings_list

    def get_close_genes(self, gene_nodes, gene_embeddings, threshold):
        #function to find genes similar to known disease genes
        #inputs:
        #   gene_nodes: list of gene names corresponding to embeddings
        #   gene_embeddings: list og embedding vectors corresponding to names
        #   threshold: max cosine distance for genes to be considered similar
        #returns: set containing known disease genes and predicted similar genes
        #compare disease genes with embeddings
        gene_nodes_array = np.array(gene_nodes)
        disease_mask = np.isin(
            gene_nodes_array,
            list(self.disease_genes)
        )
        gene_embeddings_array = np.array(gene_embeddings)
        disease_embeddings = gene_embeddings_array[disease_mask]
        distances = cosine_distances(disease_embeddings, gene_embeddings_array)
        #find genes within threshold
        close_mask = np.any(distances <= threshold, axis = 0)
        close_nodes = set(gene_nodes_array[close_mask])
        #return set
        close_genes = self.disease_genes.copy()
        return close_genes.union(close_nodes)

def main():
    # check if correct arg # in command
    if (len(sys.argv) !=3):
        return

    diseaseGeneFile = sys.argv[1]
    interactionNetworkFile = sys.argv[2]

    ppi = PPI()
    ppi.load_data(diseaseGeneFile, interactionNetworkFile)

    gene_nodes, gene_embeddings = ppi.calculate_embedding()

    results = ppi.get_close_genes(gene_nodes, gene_embeddings, 0.2)

    with open("predicted_genes.txt", "w") as f:
        f.write( "pathway_name\tdescription\t"
        + "\t".join(sorted(results))
        + "\n"
    )

    #close_genes = test_ppi.get_close_genes(nodes, embeddings, t)
    #for t in [0.05, 0.1, 0.2, 0.3]:
        #predicted_genes = ppi.get_close_genes(gene_nodes, gene_embeddings,t)
        #print(t, len(predicted_genes))
        #result = test_ppi.get_close_genes(nodes, embeddings, t)
        #print("Threshold:", t, "Number of genes:", len(result))


if __name__ == "__main__":
    main()

