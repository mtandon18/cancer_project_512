import os
import re
import sys
import pandas as pd

from network_prop import propagation
from clustering_ilp import find_clusters
import find_prop_threshold as fpt


def run_framework(priorFile, paramsFile, execDir, func=find_clusters.run, pvalue=None, expressionFile=None):
    # Run propagation
    propFile = execDir + '/prop.txt' 

    ''' Change above line of code here '''

    propagation.run(priorFile, paramsFile, propFile, expressionFile)

    if pvalue != None:
        # Get the actual size of the prior (i.e. exclude genes not in the ppi network) and run propagations over random priors of the same size
        params = propagation.read_params(paramsFile)
        priorSize = fpt.prior_size(priorFile, params['ppiFilename'])
        fpt.generate_random_data(paramsFile, [priorSize], override=False)

        # Get the pvalues for each gene
        pvaluesFile = execDir + '/ranking_by_pvalue.txt'
        fpt.calc_pvalues(propFile, fpt.get_exec_dir(paramsFile), priorSize, pvaluesFile)

        # Extract the nodes above the threshold rank into thresholded_prop.txt
        ranking = fpt.read_prop_file(pvaluesFile, index=False)
        ranking = ranking[ranking[1] < pvalue]
        thresholdedPropFile = execDir + '/thresholded_prop.txt'
        with open(thresholdedPropFile, 'w') as tpf:
            for i in range(len(ranking)):
                tpf.write(ranking.iloc[i][0] + '\t' + str(ranking.iloc[i][2]) + os.linesep)

    # Execute the ILP clustering
    os.mkdir(execDir + '/output')
    func(execDir, thresholdedPropFile if pvalue != None else propFile, paramsFile)

    # Create files that are required for convenient calculation of enrichments
    create_derived_files(execDir + '/output', priorFile)


def create_derived_files(outputDir, priorFile):
    # Read the disease prior
    with open(priorFile, 'r') as f:
        lines = f.readlines()
    prior = [line.split('\t')[0] for line in lines]

    # Write a prior node attributes file for visualization
    with open(outputDir + '/prior.na', 'w') as f:
        f.writelines(['Prior (class=java.lang.String)\n'] + [gene + ' = 1\n' for gene in prior])

    clusterFiles = [name for name in os.listdir(outputDir) if re.match('cluster_[0-9]+$', name)]
    clusteringGenes = []
    for clusterFile in clusterFiles:
        with open(outputDir + '/' + clusterFile, 'r') as fi, open(outputDir + '/wo_prior_' + clusterFile, 'w') as fo:
            genes = fi.readlines()
            genes = list(map(str.strip, genes))
            fo.write('\n'.join(set.difference(set(genes), set(prior))) + '\n')
        clusteringGenes = clusteringGenes + genes

    # Write the union of all clustering genes, w/ and w/o the prior
    with open(outputDir + '/clusters_union', 'w') as f1, open(outputDir + '/clusters_union_wo_prior', 'w') as f2:
        f1.write('\n'.join(clusteringGenes) + '\n')
        f2.write('\n'.join(set.difference(set(clusteringGenes), set(prior))) + '\n')
    

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Usage: python run_ilp_clustering.py <priorFile> <paramsFile> <execDir>')
        print('Example: python run_ilp_clustering.py $BNET_DIR/schizophrenia/schiz_genes_prior.txt $BNET_DIR/schizophrenia/params.ini $BNET_DIR/dis_comp/diseases/schiz')
        exit(1)

    run_framework(sys.argv[1], sys.argv[2], sys.argv[3], pvalue=0.01)
