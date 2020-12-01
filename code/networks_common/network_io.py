import networkx as nx
import os

def read_ppi(filename, removeSelfLoops=False):
    ''' Read a weighted ppi network in sif format, naming the weight attribute "confidence" '''
    ppiNetwork = nx.Graph()
    with open(filename, 'r') as f:
        for line in f:
            lineArr = line.strip().split('\t')
            ppiNetwork.add_edge(lineArr[0], lineArr[2], confidence=float(lineArr[1]))
    if removeSelfLoops:
        ppiNetwork.remove_edges_from(ppiNetwork.selfloop_edges())
    return ppiNetwork


def save_network(network, attrName, filename):
    ''' Save a network in sif format using the given attribute name as the middle value '''
    with open(filename, 'w') as f:
        for u,v,attrs in network.edges(data=True):
            f.write(u + '\t' + str(attrs[attrName]) + '\t' + v + os.linesep)


def save_network_default_weight(network, filename, weight=0.99):
    ''' Like save_network, but using the same weight for all edges '''
    with open(filename, 'w') as f:
        for u,v in network.edges():
            f.write(u + '\t' + str(weight) + '\t' + v + os.linesep)
            