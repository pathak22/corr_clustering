#include "WeightMatrix.h"
#include "Greedy.h"
#include "BestOneEltMove.h"

int main(int argc, char* argv[])
{
	cerr<<"reading "<<argv[1]<<"\n";
	ifstream ins(argv[1]);

	WeightMatrix mat(ins, W_LOG);

	WeightMatrix result(mat.nodes(), W_UNALTERED);

	int samples = 100;
	int edges = mat.edges();

	for(int ii = 0; ii < samples; ++ii)
	{
		VotedLink solver(mat);
		WeightMatrix& votedResult = solver.solve();

		ints clustering;
		votedResult.nodeLabels(clustering);
		BestOneEltMove boem(mat, clustering);
		WeightMatrix& boemResult = boem.solve();

		for(int item = 0; item < edges; ++item)
		{
			result[item] += boemResult[item];
		}
	}

	for(int item = 0; item < edges; ++item)
	{
		result[item] /= samples;
	}

	result.print(cout);
	result.write(cout);
}
