#include "WeightMatrix.h"

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		cerr<<"truth [original node file]\n";
		return 1;
	}

	ifstream nfile(argv[1]);
	//discard training section
	string key;
	getline(nfile,key);
	while(1)
	{
		if(key == "")
		{
			break;
		}
		getline(nfile,key);
	}

	ints labels;
	int label;
	while(nfile>>label)
	{
		getline(nfile, key);
		labels.push_back(label);
	}

	int n = labels.size();
	WeightMatrix ans(n, W_UNALTERED);
	for(int ii = 0; ii < n; ii++)
	{
		for(int jj = ii + 1; jj < n; jj++)
		{
			if(labels[ii] == labels[jj])
			{
				ans(ii, jj) = 1;
			}
		}
	}

//	ans.write(cout);
	ints matrixLabels;
	int k = ans.nodeLabels(matrixLabels);
	cerr<<"k-clustered as "<<k<<" clusters\n";
	for(int ii = 0; ii < matrixLabels.size(); ii++)
	{
		cout<<matrixLabels[ii]<<"\n";
	}
}
