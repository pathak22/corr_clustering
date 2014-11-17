class WeightMatrix;

#ifndef WEIGHT_MATRIX_H
#define WEIGHT_MATRIX_H

#include "common.h"

enum { W_UNALTERED, W_LOG };

const float EPSILON = 1e-5;

class WeightMatrix
{
public:
	WeightMatrix(int n, int mode);
	WeightMatrix(WeightMatrix & other);
	WeightMatrix(istream& is, int mode);
	~WeightMatrix();
	WeightMatrix & operator=(WeightMatrix & other);

	void read(istream& is);
	void write(ostream& os);
	void print(ostream& os);

	int nodes();
	int edges();

	int index(int row, int col);

	//the semantically sensitive matrix index
	float& operator()(int row, int col);

	//indexes straight onto the data array
	float& operator[](int idx);

	float wplus(int row, int col);
	float wminus(int row, int col);
	//wplus - wminus
	float wdiff(int row, int col);

	//returns the number of clusters, and fills the results vector
	//with a set of node labels
	//or alternately, returns -1 and an empty vector (if the matrix
	//doesn't define a clustering)
	int nodeLabels(ints& result);

private:
	int _n;
	int _mode;
	int _dSize;
	float* _data;

	void init();
};

#endif
