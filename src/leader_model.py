# ---------------------------------------------------------------------
# Model 2 : Leader Generative Graph Model
# Python module to generate data using model-1 and perform greedy pivot
# algorithm.
# 
# Note: All Indices, whatsoever, start from 0.
# ---------------------------------------------------------------------
import random
import networkx as nx
import numpy as np



# ---------------------------------------------------------------------
def gen_data(k,ni,p,epsilon,dataDir):
# ---------------------------------------------------------------------
# This function is to generate data using the leader-based generative
# graph model. The input parameters are :
# k = 5				# Number of cluster
# ni = 50      		# Number of nodes per cluster
# p = 0.2			# Intra-non-leader flipping probability
# epsilon = 0.0 	# Leader-neighbor flipping probability
# It saves the data to a file and returns the pointer to final graph

	# Inital settings 
	n = ni*k    	# Total number of nodes

	# Generate inital connected graph : g0
	g0 = nx.Graph()
	for i in range(k):
		nodeList = range(ni*i,ni*i+ni)
		edgeList = [(x,y) for x in nodeList for y in nodeList if x<y]
		g0.add_nodes_from(nodeList)
		g0.add_edges_from(edgeList)


	# Generate final graph after flipping : gf
	# Note: We assume that the first node (as per the index set) in every
	# cluster is the leader. This is without loss of generality as the graph
	# anyway has fully connected components without noise.
	nodeList = g0.nodes()
	edgeList = [(x,y) for x in nodeList for y in nodeList if x<y]
	gf = g0.copy()
	for e in edgeList:
		# If one of the vertex is a leader
		if (e[0]%ni==0) or (e[1]%ni==0):
			# Both vertices in same cluster
			if int(e[0]/ni)==int(e[1]/ni):
				if np.random.binomial(1,epsilon): 
					gf.remove_edge(e[0],e[1])

			# Both vertices in diff cluster
			else:
				if np.random.binomial(1,epsilon):
					gf.add_edge(e[0],e[1])

		# If none of the vertex is a leader
		else:
			# Both vertices in same cluster
			if int(e[0]/ni)==int(e[1]/ni):
				if np.random.binomial(1,p): 
					gf.remove_edge(e[0],e[1])

			# Both vertices in diff cluster
			else:
				if np.random.binomial(1,p):
					gf.add_edge(e[0],e[1])


	# Write the data
	fid = open(dataDir, 'w')
	fid.write('# First Line : n p epsilon\n')
	fid.write('# Following Lines : leaderOrNot groundTruthClusterID vertexID listOfNeighbours\n')
	fid.write('{} {} {}'.format(n,p,epsilon))
	for v in gf.nodes_iter():
		fid.write('\n{} {} {}'.format(int(v%ni==0),int(v/ni),v))
		for i in gf.neighbors(v):
			fid.write(' {}'.format(i))
	fid.close()

	# Return final graph
	return gf





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
def density_pivot_algorithm (gOriginal):
# ---------------------------------------------------------------------
# This algorithm is motivated from CC-Pivot (Greedy Algorithm) Algorithm 
# suggested  Ailon et. al. [STOC 2005] and again mentioned in Elsner et. al. 
# [2009]. We pick the pivot with highest density chosen by the ordering in
# starting.
#
# This function takes a graph g which is of the form networkx.Graph() 
# and returns the clustering labels for each node with attribute key as 
# 'clusterId'.
#
# This is tail recursion implementation of algorithm formatted in the 
# form of loop, till nodes() become empty in temporary graph.

	# Returns the list of vertices in decreasing order of their density
	def sort_vertices(g):
		densityList = []
		verticesList = []
		for v in g.nodes_iter():
			totalEdges = g.degree(v)*(g.degree(v)-1)/2
			presentEdges = 0
			nodeList = g.neighbors(v)
			edgeList = [(x,y) for x in nodeList for y in nodeList if x<y]
			for e in edgeList:
				if g.has_edge(e[0],e[1]):
					presentEdges = presentEdges + 1
			density = presentEdges/totalEdges if totalEdges else 0
			densityList.append(density)
			verticesList.append(v)
		sortedIndices = sorted(range(len(densityList)), key=lambda k: densityList[k], reverse=True)
		sortedVertices = [verticesList[i] for i in sortedIndices]
		# print('densityList : ',densityList,'\n verticesList:', verticesList, '\n sortedVertices:',sortedVertices)
		return sortedVertices

	'''
	# Returns the vertex with maximum density
	def max_density_vertex(g):
		maxDensity = 0
		maxDensityVertex = 0
		for v in g.nodes_iter():
			totalEdges = g.degree(v)*(g.degree(v)-1)/2
			presentEdges = 0
			nodeList = g.neighbors(v)
			edgeList = [(x,y) for x in nodeList for y in nodeList if x<y]
			for e in edgeList:
				if g.has_edge(e[0],e[1]):
					presentEdges = presentEdges + 1
			density = presentEdges/totalEdges if totalEdges else 0
			if maxDensity <= density:
				maxDensity = density
				maxDensityVertex = v
		# print('density : ',maxDensity,' vertex:', maxDensityVertex, ' v:',v)
		return maxDensityVertex
	'''

	clusterId = 500
	gNew = gOriginal.copy()

	sortedVertices = sort_vertices(gNew)
	pivotIndex = 0
	while gNew.nodes():
		# pivot = max_density_vertex(g)
		if gNew.has_node(sortedVertices[pivotIndex]):
			pivot = sortedVertices[pivotIndex]
			pivotIndex = pivotIndex + 1
		else:
			pivotIndex = pivotIndex + 1
			continue

		gOriginal.node[pivot]['clusterId'] = clusterId
		for v in gNew.neighbors(pivot):
			gOriginal.node[v]['clusterId'] = clusterId
		gNew.remove_nodes_from(gNew.neighbors(pivot))
		gNew.remove_node(pivot)
		clusterId = clusterId + 1





# ---------------------------------------------------------------------
def vote_algorithm (gOriginal,typeG):
# ---------------------------------------------------------------------
# This algorithm is the VOTE/BEST (Greedy Algorithm) Algorithm suggested in 
# Elsner et. al. [2009].
#
# This function takes a graph which is of the form networkx.Graph(), 
# and a type which is either 'best'/'vote' deciding kind of greedy
# and returns the clustering labels for each node with attribute key as 
# 'clusterId'.
	def get_quality(node,gOriginal,c,typeG):
		sumQuality = 0
		maxQuality = -10000
		for v in gOriginal.nodes_iter(data=True):
			if ('clusterId' in v[1]) and (v[1]['clusterId']==c):
				quality = 1 if gOriginal.has_edge(node,v[0]) else -1
				sumQuality = sumQuality + quality
				if quality > maxQuality:
					maxQuality = quality

		if typeG =='best':
			return maxQuality
		else:
			return sumQuality


	clusterIdStart = 1000
	clusterId = []
	nodesList = list(np.random.permutation(gOriginal.nodes()))
	for v in nodesList:
		if not clusterId:
			clusterId.append(clusterIdStart)
			gOriginal.node[v]['clusterId'] = clusterIdStart
			continue

		quality = []
		for c in clusterId:
				quality.append(get_quality(v,gOriginal,c,typeG))
		cHat = quality.index(max(quality))
		if quality[cHat]>0:
			gOriginal.node[v]['clusterId'] = clusterId[cHat]
		else:
			clusterId.append(clusterId[-1]+1)
			gOriginal.node[v]['clusterId'] = clusterId[-1]














