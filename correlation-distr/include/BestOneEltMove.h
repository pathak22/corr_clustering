class BestOneEltMove;

#ifndef BEST_ONE_ELT_MOVE_H
#define BEST_ONE_ELT_MOVE_H

#include "common.h"
#include "WeightMatrix.h"

//makes the best one-element move repeatedly until no further progress
//is possible

class BestOneEltMove
{
public:
	BestOneEltMove(WeightMatrix& wt);
	BestOneEltMove(WeightMatrix& wt, ints& initial);
	BestOneEltMove(WeightMatrix& wt, WeightMatrix& initial);

	WeightMatrix& solve();
	WeightMatrix & solveSimAnneal();
protected:
	void initCosts();
	void initClusters();
	void restoreMatrix(); //writes to _result
	//the naive algorithm: quadratic search
	//moveNode and toCluster are out parameters. Returns true if an improving move exists.
	bool findMove(int& moveNode, int& toCluster);
	//pick a random move. Return its score.
	double randomMove(int& moveNode, int& toCluster);
	void makeMove(int moveNode, int toCluster);
	
	WeightMatrix& _input;
	WeightMatrix _result;
	ints _clustering; //maps node number into cluster number ?
	intToIntSet _clusters; //maps cluster number into set of nodes in it
	intToIntToDouble _enter; //maps cluster number, node number, into (minus the cost of placing this node in this cluster)
	int _nextCluster; //equals number of clusters?
};

#endif
