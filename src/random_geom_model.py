# ---------------------------------------------------------------------
# Model 1 : Random Geometric Graph
# Python script to generate data using model-1 and perform greedy pivot
# algorithm.
# 
# Note: All Indices, whatsoever, start from 0.
# ---------------------------------------------------------------------
import random
import networkx as nx
import numpy as np


# ---------------------------------------------------------------------
# Initalize parameters
# ---------------------------------------------------------------------
seed = 222
k = 5			# Number of cluster
th = 0.3    	# GRG threshold
ni = 20      	# Number of nodes per cluster
n = ni*k    	# Total number of nodes
p = 0.2			# Intra-cluster flipping probability
q = 0.4			# Inter-cluster connected flipping probability
epsilon = 0.05 	# Inter-cluster disconnected flipping probability


# ---------------------------------------------------------------------
# Generate geometric random graph : gr
# ---------------------------------------------------------------------
random.seed(seed)
np.random.seed(seed)
gr = nx.random_geometric_graph(k,th,seed)


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
fid = open('./data/data_randomGeomModel.txt', 'w')
fid.write('{} {} {} {}'.format(n,p,q,epsilon))
for v in gf.nodes_iter():
	fid.write('\n{} {}'.format(int(v/ni),v))
	for i in gf.neighbors(v):
		fid.write(' {}'.format(i))
fid.close()



# ---------------------------------------------------------------------
def pivot_algorithm (gOriginal):
# ---------------------------------------------------------------------
# This algorithm is the CC-Pivot (Greedy Algorithm) Algorithm suggested 
# Ailon et. al. [STOC 2005] and again mentioned in Elsner et. al. [2009]. 
# For edges with signs +/- and minimization objective, this algorithm 
# is 3OPT. While for weighted edges, with w_ij^+ + w_ij^- = 1, this 
# algorithm is 5OPT.
#
# This function takes a graph g which is of the form networkx.Graph() 
# and returns the clustering labels for each node with attribute key as 
# 'clusterId'.
#
# This is tail recursion implementation of algorithm formatted in the 
# form of loop, till nodes() become empty in temporary graph.
	clusterId = 100
	gNew = gOriginal.copy()

	while gNew.nodes():
		pivot = gNew.nodes()[random.randint(0,len(gNew.nodes())-1)]
		gOriginal.node[pivot]['clusterId'] = clusterId
		for v in gNew.neighbors(pivot):
			gOriginal.node[v]['clusterId'] = clusterId
		gNew.remove_nodes_from(gNew.neighbors(pivot))
		gNew.remove_node(pivot)
		clusterId = clusterId + 1




# ---------------------------------------------------------------------
# Run Greedy and save output to file
# ---------------------------------------------------------------------
pivot_algorithm(gf)
fid = open('./data/pivot_randomGeomModel.txt', 'w')
fid.write('# Note:  GroundTruthClusterID ObtainedClusterID dont correspond, they just denote grouping of nodes.\n')
fid.write('Node GroundTruthClusterID ObtainedClusterID')
for v in gf.nodes_iter():
	fid.write('\n{} {} {}'.format(v,int(v/ni),gf.node[v]['clusterId']))
fid.close()
















