from configparser import SafeConfigParser
import os
import sys
import time
from gurobipy import *
from . import clustering_viz
from networks_common import network_io, network_utils
from network_prop import prop_file_reader

''' Find dense clusters in propagation output.
    Usage example: find_clusters C:/On work/Schizophrenia/params.ini
'''

def edge_var_name(u, v):
    return 'e_{' + u + ',' + v + '}'


def get_edge_var(m, u, v):
    var = m.getVarByName(edge_var_name(u, v))
    if var == None:
        var = m.getVarByName(edge_var_name(v, u))
    return var


def find_max_cluster(nodeScores, network, logFile, outDir, numClusterNodes, nodeWeightCoeff, mipGap, timeLimit, forceStrongConnectivity):
    numNodes = len(nodeScores)

    m = Model("Clustering")

    # Add the node variables
    nodeVars = []
    for node,score in list(nodeScores.items()):
        var = m.addVar(vtype=GRB.BINARY, name=node, obj=score*nodeWeightCoeff)
        nodeVars.append(var)

    # Add the edge variables
    edgeVars = []
    allWeights = []
    for nodeI,nodeJ,props in network.edges(data=True):
        weight = props['weight']
        allWeights.append(weight)
        var = m.addVar(vtype=GRB.CONTINUOUS, lb=0, ub=1, name=edge_var_name(nodeI,nodeJ), obj=weight)
        edgeVars.append(var)

    minWeight, maxWeight, avgWeight = min(allWeights), max(allWeights), sum(allWeights) / len(allWeights)
    logFile.write('Input network stats: min weight = ' + str(minWeight)  + '; max weight = ' + str(maxWeight) + '; average weight = ' + str(avgWeight) + os.linesep)

    # Save the model variables
    m.update()

    # Define the optimization direction
    m.modelSense = GRB.MAXIMIZE

    # Add the following family of constraints:
    # For every two nodes u,v, if they are in the cluster then set e_uv = 1; otherwise set e_uv=0.
    for nodeI,nodeJ in network.edges():
        nodeIVar = m.getVarByName(nodeI)
        nodeJVar = m.getVarByName(nodeJ)
        edgeVar = m.getVarByName(edge_var_name(nodeI, nodeJ))
        m.addConstr(edgeVar <= nodeIVar)
        m.addConstr(edgeVar <= nodeJVar)
        m.addConstr(edgeVar >= nodeIVar + nodeJVar - 1)

    numNodesInClusterExpr = quicksum(nodeVars)

    # Force each cluster node to be connected to at least half of the cluster nodes (round down).
    # IMPORTANT: if this constraint is removed, the constraint (*) must remain !
    for node in list(nodeScores.keys()):
        # If a node had degree 1 and its neighbor was selected to a cluster, then a variable will exist for that node, but it will not exist in the updated PPI network
        nodeNeighbors = set(network.neighbors(node)) if node in network else []

        # Express the sum of edge variables between node and its neighbors
        nodeEdgeVars = []
        for neighbor in nodeNeighbors:
            nodeEdgeVars.append(get_edge_var(m, node, neighbor))
        sumNodeEdgeVars = quicksum(nodeEdgeVars)

        nodeVar = m.getVarByName(node)
        
        if forceStrongConnectivity:
            m.addConstr(2*sumNodeEdgeVars + 1 >= numNodesInClusterExpr - 2*numNodes*(1-nodeVar))

        # Prevent isolated turned on nodes (*)
        m.addConstr(sumNodeEdgeVars >= nodeVar)

    # Restrict the cluster size
    m.addConstr(numNodesInClusterExpr == numClusterNodes)
    
    # Set gurobi's log file and avoid output to standard output
    m.setParam(GRB.Param.LogFile, outDir + '/gurobi.log')
    m.setParam(GRB.Param.LogToConsole, 0)

    # Set optimization parameters
    m.setParam(GRB.Param.MIPGap, mipGap)
    m.setParam(GRB.Param.TimeLimit, timeLimit)
    m.setParam(GRB.Param.Threads, 20)
    #m.setParam(GRB.Param.MIPFocus, 3)

    # Add to the objective function a penalty for the cluster non-edges
    numEdgesNotInClusterExpr = LinExpr((numClusterNodes**2 - numNodesInClusterExpr) / 2 - quicksum(edgeVars))
    m.setObjective(LinExpr(m.getObjective() - 2.3 * numEdgesNotInClusterExpr))

    # Solve the ILP
    m.update()
    m.write(outDir + '/model.lp')
    m.optimize()

    if m.status == GRB.Status.INFEASIBLE or m.status == GRB.Status.INF_OR_UNBD:
        logFile.write('ILP is infeasible' + os.linesep)
        return [], 0

    # Find the output cluster and print a report
    clusterNodes = get_ilp_cluster(nodeVars, nodeScores, logFile)

    clusterSize = len(clusterNodes)
    numPotentialEdges = clusterSize * (clusterSize-1) / 2
    numEdgesInCluster = 0

    for i in range(clusterSize):
        u = clusterNodes[i]
        for j in range(i+1, clusterSize):
            v = clusterNodes[j]
            uvEdge = network.get_edge_data(u,v)
            isInPPI = (uvEdge != None)
            if isInPPI:
                edgeExists = get_edge_var(m, u, v).x # must be 1 - printed for verification only
                edgeWeight = uvEdge['weight']
                numEdgesInCluster += 1
            else:
                edgeExists = 0
                edgeWeight = -2.3
            logFile.write('(' + u + ',' + v + ')' + '\t' + str(edgeWeight) + '\t' + str(edgeExists) + os.linesep)

    logFile.write('Number of nodes in cluster: ' + str(clusterSize) + os.linesep)
    logFile.write('Number of edges: ' + str(numEdgesInCluster) + ' out of potential ' + str(numPotentialEdges) + os.linesep)
    logFile.write('Final score: ' + str(m.objVal) + os.linesep + os.linesep)

    return clusterNodes, m.objVal


def get_ilp_cluster(nodeVars, nodeScores, logFile):
    clusterNodes = []
    for var in nodeVars:
        if var.x > 0.99: # in solution cluster, account for numerical errors
            name = var.getAttr('VarName')
            clusterNodes.append(name)
            logFile.write(name + '\t' + str(nodeScores[name]) + os.linesep)
    return clusterNodes


def perform_clustering(nodeScoresFile, outDir, logFile, minPropagationScore, maxPropagationNodes, ignoredNodesFilename, ppiFilename, numClusters, minClusterNodes, maxClusterNodes, **extras):
    startTime = time.time()
    nodeScores = prop_file_reader.read_node_scores(nodeScoresFile, minPropagationScore, maxPropagationNodes, ignoredNodesFilename)
    ppiNetwork = network_io.read_ppi(ppiFilename, removeSelfLoops=True)

    # Reweight edges by their probability in a random degree preserving graph
    reweightedNetwork = network_utils.reweight_by_prob_in_rand_graph(ppiNetwork)
    network_io.save_network(reweightedNetwork, 'weight', outDir + '/reweighted_network.txt')

    # Retain only edges between high scoring genes
    networkForClustering = reweightedNetwork.subgraph(list(nodeScores.keys()))
    network_io.save_network(networkForClustering, 'weight', outDir + '/network_for_clustering.txt')

    clusters = []

    for i in range(numClusters):
        logFile.write('Computing cluster no. ' + str(i+1) + os.linesep)

        # Find the cluster with the maximal score
        bestCluster, objVal = None, float('-inf')
        for sz in range(minClusterNodes, maxClusterNodes + 1):
            logFile.write('Invoking find_max_cluster for cluster size ' + str(sz) + os.linesep)
            bestClusterForSz, objValForSz = find_max_cluster(nodeScores, networkForClustering, logFile, outDir, sz, **extras)
            if objValForSz > objVal:
                bestCluster = bestClusterForSz
                objVal = objValForSz
                
        logFile.write('Best score for cluster no. ' + str(i+1) + ': ' + str(objVal) + os.linesep)

        # If the ILP was infeasible then stop looking for clusters
        if len(bestCluster) == 0:
            break

        clusters.append(bestCluster)

        # Write the nodes of the current cluster to a file
        clusterFile = open(outDir + '/cluster_' + str(i+1), 'w')
        for node in bestCluster:
            clusterFile.write(node + os.linesep)
        clusterFile.close()

        # Remove the cluster nodes from the propagation results and from the PPI network
        networkForClustering.remove_nodes_from(bestCluster)
        for node in bestCluster:
            nodeScores.pop(node)

    endTime = time.time()
    logFile.write('Total computation time = ' + str(endTime - startTime) + os.linesep)
    
    return clusters


def read_params(paramsFile):
    parser = SafeConfigParser(allow_no_value=True, defaults={'minPropagationScore':'0', 'maxPropagationNodes':'999999', 'forceStrongConnectivity':'1'})
    parser.read(paramsFile)
    
    params = {}
    
    # Read "input" section parameters
    params['minPropagationScore'] = parser.getfloat('input', 'minPropagationScore')
    params['maxPropagationNodes'] = parser.getint('input', 'maxPropagationNodes')
    params['ignoredNodesFilename'] = parser.get('input', 'ignoredNodesFilename')
    params['ppiFilename'] = parser.get('input', 'ppiFilename')
    
    # Read "clustering" section parameters
    params['maxClusterNodes'] = parser.getint('clustering', 'maxClusterNodes')
    params['nodeWeightCoeff'] = parser.getfloat('clustering', 'nodeWeightCoeff')
    params['numClusters'] = parser.getint('clustering', 'numClusters')
    params['forceStrongConnectivity'] = parser.getboolean('clustering', 'forceStrongConnectivity')
    
    # Read "ilp" section parameters
    params['mipGap'] = parser.getfloat('ilp', 'mipGap')
    params['timeLimit'] = parser.getint('ilp', 'timeLimit')
    
    params['minClusterNodes'] = 4
    
    return params
    

def run(execDir, nodeScoresFile, paramsFile):
    outDir = execDir + '/output'
    if not os.path.isdir(outDir):
        raise ValueError('Output dir should be created manually: '+ outDir)
    
    logFile = open(outDir + '/clustering.log', 'w')
    try:
        params = read_params(paramsFile)
        
        logFile.write('Starting fixed_size_find_clusters.py with parameters:' + os.linesep)
        logFile.write(str(params) + os.linesep)
        
        # Run the clustering algorithm
        clusters = perform_clustering(nodeScoresFile, outDir, logFile, **params)
        
        clustering_viz.prepare_cluster_visualization( \
            clusters, outDir + '/reweighted_network.txt', outDir + '/clustering.na', \
            outDir + '/clustering.sif', outDir + '/clustering_is_internal_edge.eda', outDir + '/clustering_edge_weights.eda')
    finally:
        logFile.close()


if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2], sys.argv[3])
    