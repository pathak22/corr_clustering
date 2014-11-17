#include "WeightMatrix.h"
#include "SDPSolver.h"

void usage()
{
	cerr<<"testSDPFile [classifier decisions]\n";
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

//	input.print(cerr);
//	input.write(cerr);

	SDPSolver solver(input);

	solver.writeFile(cout);
}
