import os
import sys
import clustering_common
from networks_common import network_io
import run_over_multiple_diseases
from scipy.stats.stats import pearsonr

def count_clusters(diseaseGeneFile, paramsFile, subdir, baseDir):
    diseaseGenes = run_over_multiple_diseases.parse_disease_gene_file(diseaseGeneFile)
    network = network_io.read_ppi(clustering_common.get_ppi_filename(paramsFile), removeSelfLoops=True)

    count = 0
    numClusters = 0
    clusterSizes, priorSizes = [], []

    hotnetDiseases = []

    for disease in diseaseGenes:
        diseaseDir = baseDir + '/' + disease
        if not os.path.exists(diseaseDir):
            continue

        clusteringDir = diseaseDir + '/' + subdir
        clusters, clusterEdges = clustering_common.read_clusters_with_edges(clusteringDir, network)
        print(len(clusters), disease)

        if len(clusters) > 0:
            count += 1
            if subdir == 'hotnet':
                hotnetDiseases.append(disease)

        numClusters += len(clusters)

        prior = clustering_common.read_genes_in_prior(diseaseDir + '/prior_in_network.txt')
        clusterSizes.append(len(clusters))
        priorSizes.append(len(prior))

        if len(clusters) == 0:
           print('No clusters for disease', disease)
        else:
            clustersUnion = set.union(*clusters)
            predictions = set.difference(clustersUnion, set(prior))
            if len(predictions) == 0:
                print('No predictions for disease', disease)

    print('\nCount =', count, '; Total=', len(diseaseGenes))
    print('Number of clusters =', numClusters)
    print(pearsonr(clusterSizes, priorSizes))

    if subdir == 'hotnet':
        with open(home + '/validation/DISEASES/hotnet_diseases', 'w') as f:
            for d in hotnetDiseases:
                f.write(d.replace('_', ' ') + '\t' + str(len(diseaseGenes[d])) + '\t' + '|'.join(diseaseGenes[d]) + '\n')


home = '/home/bnet/arnonmazza/dis_comp'
count_clusters(home + '/input/disease_genes_priors.txt', home + '/input/general.ini', sys.argv[1], home + '/automatic_diseases')
