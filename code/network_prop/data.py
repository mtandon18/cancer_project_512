import numpy as np
import scipy.sparse as sps
from .network import Network

def read_ppi(networkFile, expressionFile):
    # decay factor for edges over non-expressed genes
    ALPHA = 0.1

    # Read the list of expressed genes
    if expressionFile != None:
        with open(expressionFile, 'r') as f:
            expressedGenes = f.readlines()
            expressedGenes = [line.strip() for line in expressedGenes]

    # Read the PPI network and perform edge reweighting by gene expression
    fromArr = []
    toArr = []
    confidenceArr = []
    with open(networkFile, 'r') as f:
        for line in f:
            lineArr = line.strip().split('\t')
            fromNode, toNode = lineArr[0], lineArr[2]
            confidence = float(lineArr[1])

            # Perform edge reweight by tissue-specific prince paradigm
            if expressionFile != None:
                if (fromNode in expressedGenes) and (toNode in expressedGenes):
                    pass # Both endpoints are expressed - do nothing
                elif (fromNode not in expressedGenes) and (toNode not in expressedGenes):
                    confidence *= ALPHA**2
                else:
                    confidence *= ALPHA

            fromArr.append(fromNode)
            toArr.append(toNode)
            confidenceArr.append(confidence)

    geneToRow = {}
    allGenes = list(set(fromArr + toArr))
    i = 0
    for gene in allGenes:
        geneToRow[gene] = i
        i += 1

    for i in range(len(fromArr)):
        fromArr[i] = geneToRow[fromArr[i]]
        toArr[i] = geneToRow[toArr[i]]

    n = len(allGenes)
    d = sps.coo_matrix((confidenceArr, (fromArr, toArr)), shape=(n,n))
    sym = d + d.T
    return Network(sym, allGenes, geneToRow)

def read_prior(filename, network):
    prior = np.zeros((network.size,1))
    with open(filename, 'r') as f:
        for line in f:
            name, score = line.strip().split('\t')
            if name in network.geneToRow:
                prior[network.geneToRow[name]] = float(score)
            else:
                print('Warning: prior gene', name, 'does not appear in the background network, ignoring')

    return prior