#include "BestOneEltMove.h"

BestOneEltMove::BestOneEltMove(WeightMatrix& wt):
	_input(wt),
	_result(wt.nodes(), W_UNALTERED),
	_clustering(wt.nodes(), 1), //shouldn't the "1" be a "0"?
	_nextCluster(0)
{
}

BestOneEltMove::BestOneEltMove(WeightMatrix& wt, ints& initial):
	_input(wt),
	_result(wt.nodes(), W_UNALTERED),
	_clustering(initial),
	_nextCluster(0)
{
}

BestOneEltMove::BestOneEltMove(WeightMatrix& wt, WeightMatrix& initial):
	_input(wt),
	_result(wt.nodes(), W_UNALTERED),
	_nextCluster(0)
{
	int clusters = initial.nodeLabels(_clustering);
	assert(clusters != -1);
}

WeightMatrix& BestOneEltMove::solve()
{
	//make the cluster-to-nodes table
  initClusters();

	initCosts();

	int moveNode;
	int toCluster;

	int nMoves = 0;

	while(findMove(moveNode, toCluster))
	{
		makeMove(moveNode, toCluster);

		++nMoves;
	}

	cerr<<"Took "<<nMoves<<" element move steps.\n";

	//must transform back into a matrix
	restoreMatrix();

	return _result;
}

//make the cluster-to-nodes table
void BestOneEltMove::initClusters() {
  int n = _input.nodes();
  for(int node = 0; node < n; ++node)
    {
      _clusters[_clustering[node]].insert(node);
      if(_clustering[node] >= _nextCluster)
	{
	  _nextCluster = _clustering[node] + 1;
	}
    }
}

void BestOneEltMove::restoreMatrix() {
    foreach(intToIntSet, cluster, _clusters)
    {
      foreach(intSet, node, cluster->second)
	{
	  foreach(intSet, other, cluster->second)
	    {
	      if(*node != *other)
		{
		  _result(*node, *other) = 1;
		}
	    }
	}
    }
}

WeightMatrix& BestOneEltMove::solveSimAnneal()
{

  initClusters();	//make the cluster-to-nodes table
  initCosts();
  
  int moveNode;
  int toCluster;
  
  int nMoves = 0;
  
  int n = _input.nodes();


  const int nItr = 2000 * n;

  double tpr = 10 * sqrt(n); //was: n. The current temperature
  double tempReduction = pow(tpr / 0.001, -1./nItr); //temperature ends at 0.001 (accepts little but improving moves)
  
  for (int it = 0; it < nItr; it++) {
    double score = randomMove(moveNode, toCluster);
    tpr *= tempReduction;
    if (score >= 0 || rand() <= RAND_MAX * exp(score / tpr)) {
	  makeMove(moveNode, toCluster);
	  ++nMoves;
    }
    double obj = 0;
    for (int node = 0; node < n; node++)
      obj -= _enter[_clustering[node]][node] / 2.;

    int reportPeriod = nItr / 40;
    if (it % reportPeriod == 0)
      cerr << "Tmprtr: " << tpr << " score: " << score << " obj: " << obj << " clusters: " << _clusters.size() << std::endl;
  }
  
  cerr<<"Took "<<nMoves<<" element move steps.\n";
  
  //must transform back into a matrix
  restoreMatrix();
  
  return _result;
}

void BestOneEltMove::initCosts()
{
	int n = _input.nodes();

	foreach(intToIntSet, cluster, _clusters)
	{
		for(int node = 0; node < n; ++node)
		{
			//score obtained by having node in this cluster
			double enter = 0;

			intSet& cl = _clusters[cluster->first];
			foreach(intSet, other, cl)
			{
				if(*other != node)
				{
					enter += _input.wdiff(node, *other);
				}
			}
			//enter term collects edges which would be joined
			//if we moved this node in
			_enter[cluster->first][node] = enter;
		}
	}

	//move has goodness (enter(item, to) - enter(item, from))
}

bool BestOneEltMove::findMove(int& moveNode, int& toCluster)
{
	moveNode = -1;
	toCluster = -1;

	int n = _input.nodes();

	double best = 0; //find a *positive* max

	foreach(intToIntSet, cluster, _clusters)
	{
	  intToDouble& entries = _enter[cluster->first]; //"first" is STL-speak for the key, i.e. the cluster number

		for(int node = 0; node < n; ++node)
		{
			int currCluster = _clustering[node];

			double score = entries[node] - _enter[currCluster][node];

			if(score > best)
			{
				best = score;
				moveNode = node;
				toCluster = cluster->first;
			}
		}
	}

	for(int node = 0; node < n; ++node)
	{
		int currCluster = _clustering[node];

		double score = 0 - _enter[currCluster][node];

		if(score > best)
		{
			best = score;
			moveNode = node;
			toCluster = _nextCluster;
		}
	}

	return (moveNode != -1);
}

double BestOneEltMove::randomMove(int & moveNode, int & toCluster) {
	int n = _input.nodes();
	moveNode = int(double(rand()) * n / (RAND_MAX+1.));

	//choose which cluster to move this to. Either one of the other clusters, or the new one.
	int clustersToSkip = int(double(rand()) * (_clusters.size()) / (RAND_MAX+1.));
	

	toCluster = -1;
	
	if (clustersToSkip == _clusters.size() - 1) {
	  //make a new cluster
	  toCluster = _nextCluster;
	  int currCluster = _clustering[moveNode];
	  
	  double score = 0 - _enter[currCluster][moveNode];
	  return score;
	} else {
	  foreach(intToIntSet, cluster, _clusters)
	    {
	      int clusterNum = cluster->first; 
	      if (clusterNum != _clustering[moveNode]) {
		if (clustersToSkip == 0) {
		  toCluster = cluster->first;
		  int currCluster = _clustering[moveNode];
		  
		  double score = _enter[toCluster][moveNode] - _enter[currCluster][moveNode];
		  return score;
		}
		clustersToSkip--;
	      }
	    }
	}

	assert(false);
}

void BestOneEltMove::makeMove(int moveNode, int toCluster)
{
	int fromCluster = _clustering[moveNode];

	if(toCluster == _nextCluster)
	{
		++_nextCluster;
	}

	_clustering[moveNode] = toCluster;
	_clusters[fromCluster].erase(moveNode);
	if(_clusters[fromCluster].empty())
	{
		_clusters.erase(fromCluster);
		_enter.erase(fromCluster);
	}
	_clusters[toCluster].insert(moveNode);

	intToDouble& oldEntries = _enter[fromCluster];
	intToDouble& newEntries = _enter[toCluster];
	int n = _input.nodes();

	for(int node = 0; node < n; ++node)
	{
		oldEntries[node] -= _input.wdiff(moveNode, node);
		newEntries[node] += _input.wdiff(moveNode, node);
	}
}
