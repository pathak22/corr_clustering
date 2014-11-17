#include "WeightMatrix.h"
#include "BestOneEltMove.h"

void usage()
{
	cerr<<"testBOEM [classifier decisions]\n";
}

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		usage();
		return 1;
	}

	ifstream decisions(argv[1]);
	WeightMatrix input(decisions, W_UNALTERED);

	BestOneEltMove solver(input);

	WeightMatrix& result = solver.solve();

//	result.print(cerr);

	ints matrixLabels;
	int k = result.nodeLabels(matrixLabels);
	if(k == -1)
	{
		cerr<<"error: failed to cluster...\n";
	}

	cerr<<"k-clustered as "<<k<<" clusters\n";
	for(int ii = 0; ii < matrixLabels.size(); ii++)
	{
		cout<<matrixLabels[ii]<<"\n";
	}
}
