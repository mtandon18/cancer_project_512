from configparser import SafeConfigParser
import os
import sys
import numpy as np
from . import data

def read_params(paramsFile):
    parser = SafeConfigParser(allow_no_value=True, defaults={'iterations':'1000', 'epsilon':'1e-5'})
    parser.read(paramsFile)
    
    params = {}
    
    propSection = 'prop'
    
    # Read "input" section parameters
    params['ppiFilename'] = parser.get('input', 'ppiFilename')
    params['alpha'] = parser.getfloat(propSection, 'alpha') # relative weight for the network
    params['iterations'] = parser.getint(propSection, 'iterations') # maximum number of iterations to execute
    params['epsilon'] = parser.getfloat(propSection, 'epsilon') # convergence threshold
    
    return params


def is_valid_args(priorFile, expressionFile, ppiFilename, alpha, iterations, epsilon):
    failed = False
    
    if not os.path.exists(priorFile):
        print("Prior file does not exist")
        failed = True
        
    if expressionFile is not None and not os.path.exists(expressionFile):
        print("Gene expression file does not exist")
        failed = True

    if ppiFilename is None or not os.path.exists(ppiFilename):
        print("Network file does not exist")
        failed = True
        
    if not 0 < alpha < 1:
        print("alpha should be in (0,1)")
        failed = True

    if epsilon > 0.5:
        print("Epsilon is too large")
        failed = True

    return not failed


def validate(result, validationFile, network):
    validation = data.read_prior(validationFile, network)
    
    # Print the scores of the genes in the validation set
    with open(validationFile, 'r') as f:
        for line in f:
            gene, _ = line.strip().split('\t')
            geneRow = network.geneToRow[gene]
            print(gene, result[geneRow])
            
    # Calculate an ROC curve
    import sklearn.metrics as metrics
    fpr, tpr, _ = metrics.roc_curve(validation > 0, result)
    auc = metrics.auc(fpr, tpr)
    
    print("fpr:", fpr, '; tpr:', tpr, '; auc:', auc)


def run_propagation(priorFile, expressionFile, outFile, ppiFilename, alpha, iterations, epsilon):
    network = data.read_ppi(ppiFilename, expressionFile)
    prior = data.read_prior(priorFile, network)
    
    if np.any(prior > 1.0):
        sys.stderr.write("Scaling prior, values outside [0,1]\n")
        prior /= np.max(prior)

    result = network.normalize().smooth(
        prior, alpha=alpha, eps=epsilon, max_iter=iterations)
    
    with open(outFile, 'w') as f:
        for item in sorted(zip(network.allGenes, result.flat), key=lambda x: x[1], reverse=True):
            f.write("\t".join(map(str, item)) + '\n')


def run(priorFile, paramsFile, outFile, expressionFile=None):
    params = read_params(paramsFile)
    succeeded = is_valid_args(priorFile, None, **params)
    
    if succeeded:
        run_propagation(priorFile, expressionFile, outFile, **params)
    

if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2], sys.argv[3])
    