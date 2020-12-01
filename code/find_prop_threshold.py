from configparser import SafeConfigParser
import glob
import itertools
import math
import multiprocessing as mp
import os
import random
import sys
import pandas as pd
from networks_common import network_io
from network_prop import propagation

REPEATS = list(range(1000))


def read_prop_file(propFile, index):
    return pd.read_csv(propFile, sep='\t', index_col=index, lineterminator='\n', header=None, na_values=['X'*20], keep_default_na=False)


def get_exec_dir(paramsFile):
    parser = SafeConfigParser(allow_no_value=True)
    parser.read(paramsFile)
    return '../find_prop_th/' + parser.get('framework', 'networkName')


# Get the number of genes in the prior file that are also in the network file.
# The priorFile format can be genes only or genes with prior score.
def prior_size(priorFile, networkFile):
    data = read_prop_file(priorFile, index=False)
    priorNodes = data[0]
    networkNodes = get_network_nodes(networkFile)
    return len(set.intersection(set(priorNodes), set(networkNodes)))


def get_network_nodes(networkFile):
    # Return the nodes in the given network, excluding nodes whose only neighbors are themselves
    network = network_io.read_ppi(networkFile, removeSelfLoops=True)
    networkNodes = [node for node in network.nodes() if len(network.neighbors(node)) > 0]
    return networkNodes


def run_random_propagation(xxx_todo_changeme):
    # Randomize prior
    (paramsFile, networkNodes, priorSize, repeat) = xxx_todo_changeme
    randNodes = random.sample(networkNodes, priorSize)

    execDir = get_exec_dir(paramsFile)

    # Set output filenames
    suffix = 'S' + str(priorSize) + '_' + str(repeat+1)
    priorFile = execDir + '/prior_' + suffix
    propFile = execDir + '/prop_' + suffix

    # Save prior to disk    
    with open(priorFile, 'w') as f:
        for node in randNodes:
            f.write(node + '\t' + '1' + os.linesep)

    # Run propagation
    propagation.run(priorFile, paramsFile, propFile)


# Create two summary files:
# prior summary: per gene, 1 if it appears in the random prior for test (column) #i, 0 otherwise
# prop summary: per gene, its score in test (column) #i
#
def aggregate_data_by_prior_size(priorSizes, paramsFile, execDir):
    data = []
    for priorSize in priorSizes:
        # Set prefixes of file names
        propFilePrefix = execDir + '/prop_S' + str(priorSize) + '_'
        priorFilePrefix = execDir + '/prior_S' + str(priorSize) + '_'

        # Aggregate data from all prior and prop files
        propData, priorData = [], []
        for repeat in REPEATS:
            propData.append(read_prop_file(propFilePrefix + str(repeat+1), 0))
            priorData.append(read_prop_file(priorFilePrefix + str(repeat+1), 0))

        allGenes = propData[0].index

        # Set output file names
        propSumFile, priorSumFile = execDir + '/prop_summary_S' + str(priorSize), execDir + '/prior_summary_S' + str(priorSize)

        with open(propSumFile, 'w') as f1, open(priorSumFile, 'w') as f2:
            for gene in allGenes:
                f1.write(str(gene) + '\t')
                f2.write(str(gene) + '\t')
                for repeat in REPEATS:
                    currScore = propData[repeat].loc[gene][1]
                    f1.write(format(currScore, '.5f') + ('\t' if repeat < REPEATS[-1] else os.linesep))
                    f2.write('1' if gene in priorData[repeat].index else '0')
                    f2.write('\t' if repeat < REPEATS[-1] else os.linesep)

        # Remove the singleton files
        for f in glob.glob(propFilePrefix + '*'):
            os.remove(f)
        
        for f in glob.glob(priorFilePrefix + '*'):
            os.remove(f)


def generate_random_data(paramsFile, priorSizes, override=False):
    execDir = get_exec_dir(paramsFile)

    done = True
    for prSz in priorSizes:
        if not os.path.exists(execDir + '/prop_summary_S' + str(prSz)) or override:
            done = False
    if done:
        return

    params = propagation.read_params(paramsFile)
    networkFile = params['ppiFilename']

    pool = mp.Pool(mp.cpu_count() - 1)

    # Run random propagations for all <prior sizes> * <number of repeats> combinations
    networkNodes = get_network_nodes(networkFile)
    jobParams = [(paramsFile, networkNodes, prSz, rpt) for prSz,rpt in itertools.product(priorSizes, REPEATS)]
    pool.map(run_random_propagation, jobParams)

    # Aggregate into a single prior file and a single prop scores file
    aggregate_data_by_prior_size(priorSizes, paramsFile, execDir)

    pool.close()
    pool.join()


def calc_pvalues(propFile, randomPropsDir, priorSize, outFile):
    propData = read_prop_file(propFile, 0)
    randomProps = read_prop_file(randomPropsDir + '/prop_summary_S' + str(priorSize), 0)
    randomPriors = read_prop_file(randomPropsDir + '/prior_summary_S' + str(priorSize), 0)

    d = {}
    for gene in propData.index:
        numLower = sum((propData.loc[gene][1] > randomProps.loc[gene]) & (randomPriors.loc[gene] == 0))
        denom = float(sum(randomPriors.loc[gene] == 0))
        d[gene] = 1 - numLower/(denom+1)

    sorted_pvalues = sorted(list(d.items()), key=lambda x:x[1])
    with open(outFile, 'w') as f:
        for gene,pvalue in sorted_pvalues:
            f.write(str(gene) + '\t' + str(pvalue) + '\t' + str(propData.loc[gene][1]) + os.linesep)


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print('Usage: python find_prop_threshold.py <propFile> <priorFile> <paramsFile> <outFile>')
        print('Example: python find_prop_threshold.py $BNET_DIR/schizophrenia/propagation_ranking_with_expression.txt $BNET_DIR/schizophrenia/schiz_genes.txt $BNET_DIR/dis_comp/find_prop_th/schiz/params.ini $BNET_DIR/dis_comp/find_prop_th/schiz/ranking_by_pvalue')
        exit(1)

    propFile = sys.argv[1]
    priorFile = sys.argv[2] # only for getting the size
    paramsFile = sys.argv[3]
    outFile = sys.argv[4]

    params = propagation.read_params(paramsFile)
    priorSize = prior_size(priorFile, params['ppiFilename'])

    #PRIOR_SIZES = [5, 10, 25, 50, 100, 150, 200, 250]
    generate_random_data(paramsFile, [priorSize], override=True)

    randomPropsDir = get_exec_dir(paramsFile)
    calc_pvalues(propFile, randomPropsDir, priorSize, outFile)
