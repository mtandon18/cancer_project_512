import itertools
import os
from networks_common import network_io

def prepare_cluster_visualization(clusters, inWeightedNetworkFile, outNaFile, outClusteringFile, outIsInternalEdaFile, outEwEdaFile):
    ''' Save the network induced by the clustering as well as node and edge attribute files for visualization '''
    # Write node attributes file: for each node, "node = <cluster-id>", to be used for coloring
    with open(outNaFile, 'w') as f:
        f.write('Cluster (class=java.lang.String)' + os.linesep)
        clusteringMap = {}
        for i in range(len(clusters)):
            for node in clusters[i]:
                clusteringMap[node] = i
                f.write(node + ' = ' + str(i) + os.linesep)

    # Get the ppi network induced by the clustered nodes
    clusteringNetwork = network_io.read_ppi(inWeightedNetworkFile, removeSelfLoops=True).subgraph(list(itertools.chain(*clusters)))

    # Save the following data:
    # - the clustering network
    # - edge attributes file - for each edge, is it internal or cross clusters ?
    # - edge attributes file - for each edge, its score after the reweighting
    #
    with open(outClusteringFile, 'w') as f1, \
         open(outIsInternalEdaFile, 'w') as f2, \
         open(outEwEdaFile, 'w') as f3:

        f2.write('IsInternalEdge (class=java.lang.String)' + os.linesep)
        f3.write('EdgeWeight (class=java.lang.Double)' + os.linesep)

        for u,v,attrs in clusteringNetwork.edges(data=True):
            f1.write(u + '\t1\t' + v + os.linesep)
            inSameCluster = 1 if clusteringMap[u] == clusteringMap[v] else 0
            f2.write(u + ' (1) ' + v + ' = ' + str(inSameCluster) + os.linesep)
            f3.write(u + ' (1) ' + v + ' = ' + str(attrs['confidence']) + os.linesep)
