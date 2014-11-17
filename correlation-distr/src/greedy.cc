#include "WeightMatrix.h"
#include "Greedy.h"

void usage()
{
	cerr<<"greedy [-a algorithm|-l (list)] [classifier decisions]\n";
}

int main(int argc, char* argv[])
{
	srand(time(NULL));

	if(argc < 2)
	{
		usage();
		return 1;
	}

	Greedy* solver = NULL;

	if(string(argv[1]) == "-l")
	{
		cout<<"list:\n";
		cout<<"First-link (Soon): first\n";
		cout<<"Best-link (Ng): best\n";
		cout<<"Voted-link (Elsner?): vote\n";
		cout<<"Pivot (Ailon): pivot\n";
		return 1;
	}
	
	if(string(argv[1]) != "-a")
	{
		usage();
		return 1;
	}

	if(argc < 4)
	{
		usage();
		return 1;
	}

	ifstream decisions(argv[3]);
	WeightMatrix input(decisions, W_UNALTERED);

	string alg(argv[2]);

	if(alg == "first")
	{
		solver = new Soon(input);
	}
	else if(alg == "best")
	{
		solver = new BestLink(input);
	}
	else if(alg == "vote")
	{
		solver = new VotedLink(input);
	}
	else if(alg == "pivot")
	{
		solver = new Pivot(input);
	}
	else
	{
		cerr<<"Unrecognized algorithm: "<<alg<<"\n";
		return 1;
	}

	WeightMatrix& ans = solver->solve();
	ints matrixLabels;
	int k = ans.nodeLabels(matrixLabels);
	if(k == -1)
	{
		cerr<<"error: failed to cluster...\n";
		ans.print(cerr);
	}

	cerr<<"k-clustered as "<<k<<" clusters\n";
	for(int ii = 0; ii < matrixLabels.size(); ii++)
	{
		cout<<matrixLabels[ii]<<"\n";
	}
}
