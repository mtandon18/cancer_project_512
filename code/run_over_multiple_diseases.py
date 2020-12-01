import os
import shutil
import sys
import run_ilp_clustering
import find_prop_threshold as fpt
from network_prop import propagation

def parse_disease_gene_file(diseaseGeneFile):
    d = {}
    with open(diseaseGeneFile) as f:
        for line in f:
            line = line.strip().split('\t')
            diseaseName = line[0].replace(' ', '_')
            genes = line[2].split('|')
            d[diseaseName] = genes
    return d


def prior_in_network(priorGenes, networkFile):
    networkNodes = fpt.get_network_nodes(networkFile)
    return set.intersection(set(priorGenes), set(networkNodes))


def create_prior_file(diseaseDir, genes, filename):
    priorFile = diseaseDir + '/' + filename
    with open(priorFile, 'w') as f:
        for gene in genes:
            f.write(gene + '\t1' + os.linesep)
    return priorFile


def run(diseaseGeneFile, baseDir, paramsFile, override=False):
    logFile = open(baseDir + '/exec.log', 'w')

    # Parse disease-gene association file
    diseaseGenesMap = parse_disease_gene_file(diseaseGeneFile)

    for disease, genes in list(diseaseGenesMap.items()):
        logFile.write('Working on disease: ' + disease + '\n')
        diseaseDir = baseDir + '/' + disease

        # If disease dir already exists, then override/skip according to the value of the parameter override
        if os.path.exists(diseaseDir):
            if override:
                logFile.write('Disease folder ' + diseaseDir + ' already exists, OVERRIDING\n')
                shutil.rmtree(diseaseDir)
            else:
                logFile.write('Disease folder ' + diseaseDir + ' already exists, SKIPPING\n')
                continue

        # Find the actual prior size after intersecting with the network genes
        params = propagation.read_params(paramsFile)
        priorInNetwork = prior_in_network(genes, params['ppiFilename'])

        if len(priorInNetwork) > 0:
            os.mkdir(diseaseDir)
            priorFile = create_prior_file(diseaseDir, genes, 'prior.txt')
            create_prior_file(diseaseDir, priorInNetwork, 'prior_in_network.txt')
            run_ilp_clustering.run_framework(priorFile, paramsFile, diseaseDir, pvalue=0.01)
        else:
            logFile.write('Disease ' + disease + ' has no prior genes in the network\n')

    logFile.close()


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('Usage: python run_over_multiple_diseases.py <diseaseGeneFile> <baseDir> <paramsFile>')
        print('Example: python run_over_multiple_diseases.py $BNET_DIR/disease_complexes/disease_genes/prior/disease_genes_FINAL.txt $BNET_DIR/disease_complexes/out $BNET_DIR/disease_complexes/params.ini')
        exit(1)

    run(sys.argv[1], sys.argv[2], sys.argv[3])
