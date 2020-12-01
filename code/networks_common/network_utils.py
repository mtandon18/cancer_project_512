import math
import networkx as nx
from networks_common import network_io

BETA = 0.95

'''
We want to inject values for beta that reward interactions for similar expression
'''

def calc_weight(i, j, network, numEdges, edgeExists):
    degI, degJ = network.degree(i), network.degree(j)
    prob = degI * degJ / (2 * float(numEdges))
    if prob >= BETA:
        prob = BETA
    return math.log(BETA / prob) if edgeExists else math.log((1-BETA) / (1-prob))
    

def reweight_by_prob_in_rand_graph(network, addNonEdges=True):
    ''' Return a new network where each edge is weighted according to its (approximated) probability in a random degree preserving graph.
        This weight will be stored as the "weight" attribute.
        If addNonEdges=True then non-edges are also weighted according to the same scheme, and the resulting network is complete.
        In order to discern between edges and non-edges, an additional attribute "edgeExists" is set (1/0).
    '''
    numEdges = network.number_of_edges()
    result = nx.Graph()
    
    for i,j in network.edges():
        weight = calc_weight(i, j, network, numEdges, True)
        result.add_edge(i, j, weight=weight, edgeExists=1)
        
    if addNonEdges:
        for i,j in nx.complement(network).edges():
            weight = calc_weight(i, j, network, numEdges, False)
            result.add_edge(i, j, weight=weight, edgeExists=0)

    return result


def contract_edge(network, edge):
    '''
    a = nx.Graph()
    a.add_edge('a', 'b')
    a.add_edge('a', 'c')
    a.add_edge('c', 'd')
    a.add_edge('c', 'e')
    contract_edge(a, ('a','c'))
    print a.edges()
    '''
    node = edge[0]
    other = edge[1]
    if network.get_edge_data(node, other) == None:
        raise ValueError(str.format('contract_edge: edge {0} not in network', edge))

    newNode = node + '_' + other

    for neighbor in network.neighbors(node):
        conf = network.get_edge_data(node, neighbor).get('confidence')
        network.add_edge(newNode, neighbor, confidence=conf)

    for neighbor in network.neighbors(other):
        conf = network.get_edge_data(other, neighbor).get('confidence')
        network.add_edge(newNode, neighbor, confidence=conf)

    network.remove_node(node)
    network.remove_node(other)
    network.remove_edge(newNode, newNode)
    
    
def merge_networks(networkFiles):
    networks = []
    networkNames = []
    combinedNetwork = nx.Graph()
    
    with open(networkFiles) as f:
        for line in f:
            line = line.strip()
            network = network_io.read_ppi(line)
            networks.append(network)
            networkNames.append(line[line.rfind('/')+1 :])
            combinedNetwork.add_edges_from(network.edges())
            
    return combinedNetwork, networks, networkNames
