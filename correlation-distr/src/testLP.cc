#include "WeightMatrix.h"
#include "LPSolver.h"

void usage()
{
	cerr<<"testLP [-lp|-ilp|-lazy|-lazyilp] [classifier decisions]\n";
}

int main(int argc, char* argv[])
{
	if(argc < 3)
	{
		usage();
		return 1;
	}

	ifstream decisions(argv[2]);
	WeightMatrix input(decisions, W_UNALTERED);

	LPSolver* solver = NULL;

	string mode(argv[1]);
	if(mode == "-lp")
	{
		cerr<<"Using lp (integral solution not guaranteed)\n";
		solver = new LPSolver(input);
	}
	else if(mode == "-ilp")
	{
		cerr<<"Using ilp\n";
		solver = new ILPSolver(input);
	}
	else if(mode == "-lazy")
	{
		cerr<<"Using lazy lp (integral solution not guaranteed)\n";
		solver = new LazyLPSolver(input, 1.0);
	}
	else if(mode == "-lazyilp")
	{
		cerr<<"Using lazy ilp\n";
		solver = new LazyILPSolver(input);
	}
	else
	{
		cerr<<"Unrecognized mode "<<mode<<"\n";
		usage();
		return 1;
	}

	WeightMatrix& result = solver->solve();

	result.print(cout);

	ints matrixLabels;
	int k = result.nodeLabels(matrixLabels);
	if(k != -1)
	{
		cout<<"No rounding required!\n";

		for(int ii = 0; ii < matrixLabels.size(); ii++)
		{
			cout<<matrixLabels[ii]<<"\n";
		}
	}
	else
	{
		cout<<"Solution isn't integral, or some problem occurred.\n";
	}
}
