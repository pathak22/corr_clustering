class Greedy;

#ifndef GREEDY_H
#define GREEDY_H

#include "common.h"
#include "WeightMatrix.h"

//supports various greedy algorithms based on searching nodes in a random
//order

class Greedy
{
public:
	Greedy(WeightMatrix& wt);
	virtual ~Greedy();

	//reset the internal permutation to the identity
	//note: makes sense for linearly-ordered datasets like chat,
	//but is an egregious form of cheating on block-ordered datasets
	//like 20 newsgroups
	virtual void setIdentityPerm();
	virtual WeightMatrix& solve();

protected:
	WeightMatrix& _input;
	WeightMatrix _result;
	ints _perm;

	virtual void markNode(int node, int permIdx) = 0;
	virtual void assoc(int node, int other);
};

class Soon : public Greedy
{
	//implements the simple greedy algorithm of Soon et al '01
	//for each node, search all previous nodes:
	//  stop as soon as an edge with p > .5 is found, and link the two
	//if not linked, form a new cluster

public:
	Soon(WeightMatrix& wt);

protected:
	virtual void markNode(int node, int permIdx);
};

class BestLink : public Greedy
{
	//implements the best-link clustering algorithm of Ng + Cardie '02
	//for each node, search all previous nodes:
	//  take the best edge with p > .5 and link the two
	//if not linked, form a new cluster

public:
	BestLink(WeightMatrix& wt);

protected:
	virtual void markNode(int node, int permIdx);
};

class VotedLink : public Greedy
{
	//implements the voted-link clustering algorithm used in 
	//Elsner + Charniak '08
	//for each node, search all previous nodes, adding their edge weight
	//to the score for their cluster
	//  take the best cluster with p > .5 and link the two
	//if not linked, form a new cluster

public:
	VotedLink(WeightMatrix& wt);
	virtual WeightMatrix& solve();

protected:
	virtual void markNode(int node, int permIdx);

	ints _clusters;
	int _nextCluster;
};

class Pivot : public Greedy
{
	//implements the pivoting algorithm CC-Pivot of 
	//Ailon + Charikar + Newman '05 (by Thm 10, this is a 5-approximation)
	//for each node, take all neighboring nodes with p > .5 that remain
	//unclustered

public:
	Pivot(WeightMatrix& wt);

protected:
	virtual void markNode(int node, int permIdx);

	ints _taken;
};

#endif
