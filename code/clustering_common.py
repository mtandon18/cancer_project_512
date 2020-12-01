from configparser import SafeConfigParser
import itertools
import os
import re

# Read all clusters under clusteringDir; if includePrior=False exclude prior genes
def read_clusters(clusteringDir, includePrior=True):
    if includePrior:
        clusterFiles = [name for name in os.listdir(clusteringDir) if re.match('cluster_[0-9]+$', name)]
    else:
        clusterFiles = [name for name in os.listdir(clusteringDir) if re.match('wo_prior_cluster_[0-9]+$', name)]
    clusters = []
    for name in clusterFiles:
        with open(clusteringDir + '/' + name, 'r') as f:
            genes = f.readlines()
        genes = set(map(str.strip, genes))
        clusters.append(genes)
    return clusters


def read_clusters_with_edges(clusteringDir, network):
    clusters = read_clusters(clusteringDir)
    clusterEdges = [get_edges_in_gene_set(cluster, network) for cluster in clusters]
    return clusters, clusterEdges


def get_edges_in_gene_set(geneSet, network):
    ppis = []
    for u,v in itertools.combinations(geneSet, 2):
        if network.has_edge(u,v):
            ppis.append((u,v))
    return ppis


def get_ppi_filename(paramsFile):
    parser = SafeConfigParser(allow_no_value=True)
    parser.read(paramsFile)
    return parser.get('input', 'ppiFilename')


# Return the set of genes in the given 2-column gene-score prior file
def read_genes_in_prior(filename):
    with open(filename) as f:
        genes = [line.strip().split('\t')[0] for line in f]
    return set(genes)
