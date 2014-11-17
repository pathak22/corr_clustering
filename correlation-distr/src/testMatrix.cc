#include "WeightMatrix.h"

int main(int argc, char* argv[])
{
	{
		WeightMatrix test1(2, W_UNALTERED);
		test1.print(cout);
	}

	{
		if(argc > 1)
		{
			cerr<<"reading "<<argv[1]<<"\n";
			ifstream ins(argv[1]);
			WeightMatrix test2(ins, W_UNALTERED);
			test2.print(cout);
			test2.write(cout);
		}
	}

}
