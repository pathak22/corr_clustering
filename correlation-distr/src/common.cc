#include "common.h"

void randPermutation(ints& iperm, int size)
{
	iperm.clear();
	for(int i=0; i<size; i++)
	{
		iperm.push_back(i);
	}

	random_shuffle(iperm.begin(), iperm.end()); 	//STL random shuffle
}

bool isIdentity(ints& perm)
{
	int size = perm.size();
	for(int i = 0; i < size; i++)
	{
		if(perm[i] != i)
		{
			return false;
		}
	}
	return true;
}

string intToString(int x)
{
	ostringstream cvt;
	cvt<<x;
	return cvt.str();
}
