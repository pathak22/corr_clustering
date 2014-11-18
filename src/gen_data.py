# ---------------------------------------------------------------------
# Python script to generate data using geometric random graph
# ---------------------------------------------------------------------
import random
import networkx as nx
import numpy as np


# ---------------------------------------------------------------------
# Initalize parameters
# ---------------------------------------------------------------------
random.seed(222)
k = 10			# Number of cluster
th = 0.3    	# GRG threshold
ni = 50      	# Number of nodes per cluster
n = ni*k    	# Total number of nodes
p = 0.2			# Intra-cluster flipping probability
q = 0.4			# Inter-cluster connected flipping probability
epsilon = 0.05 	# Inter-cluster disconnected flipping probability


# ---------------------------------------------------------------------
# Generate geometric random graph : gr
# ---------------------------------------------------------------------
gr = nx.random_geometric_graph(k,th)


# ---------------------------------------------------------------------
# Generate inital connected graph : g0
# ---------------------------------------------------------------------
g0 = nx.Graph()
for i in range(k):
	nodeList = range(ni*i,ni*i+ni)
	edgeList = [(x,y) for x in nodeList for y in nodeList if x<y]
	g0.add_nodes_from(nodeList)
	g0.add_edges_from(edgeList)


# ---------------------------------------------------------------------
# Generate final graph after flipping : gf
# ---------------------------------------------------------------------
nodeList = g0.nodes()
edgeList = [(x,y) for x in nodeList for y in nodeList if x<y]
gf = g0.copy()
for e in edgeList:
	# Both vertices in same cluster
	if int(e[0]/ni)==int(e[1]/ni):
		if np.random.binomial(1,p): 
			gf.remove_edge(e[0],e[1])

	# Both vertices in diff cluster and connected in gr
	elif gr.has_edge(e[0],e[1]):
		if np.random.binomial(1,q):
			gf.add_edge(e[0],e[1])

	# Both vertices in diff cluster and not-connected in gr
	else:
		if np.random.binomial(1,epsilon):
			gf.add_edge(e[0],e[1])


# ---------------------------------------------------------------------
# Write the data
# ---------------------------------------------------------------------
fid = open('../data/data.txt', 'w')
fid.write('{} {} {} {}'.format(n,p,q,epsilon))
for v in gf.nodes_iter():
	fid.write('\n{} {}'.format(int(v/ni),v))
	for i in gf.neighbors(v):
		fid.write(' {}'.format(i))
fid.close()




