def read_node_scores(scoresFilename, minScore, maxNodes, ignoredFilename):
    ''' Return a node->score dictionary such that:
        1. only nodes with score > minScore are included
        2. no more than maxNodes nodes are returned
        3. ignored nodes are excluded (e.g. hubs)
    '''
    ignored = set()
    if ignoredFilename != None:
        with open(ignoredFilename, 'r') as f:
            for line in f:
                ignored.add(line.strip())

    nodeScores = {}
    numNodes = 0
    with open(scoresFilename, 'r') as f:
        for line in f:
            lineArr = line.strip().split('\t')
            node = lineArr[0]
            score = float(lineArr[1])
            if (score > minScore) and (numNodes < maxNodes) and (node not in ignored):
                nodeScores[node] = score
                numNodes += 1
    return nodeScores
