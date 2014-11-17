#include "WeightMatrix.h"

#include <assert.h>
#include <iomanip>

WeightMatrix::WeightMatrix(int n, int mode):
	_n(n),
	_mode(mode)
{
	init();
}

WeightMatrix::WeightMatrix(WeightMatrix & other):
	_n(other._n),
	_mode(other._mode)
{
	(*this) = other;
}

WeightMatrix::WeightMatrix(istream& is, int mode):
	_n(0),
	_mode(mode)
{
	read(is);
}

void WeightMatrix::init()
{
	_dSize = (_n * (_n-1))/2;
	_data = new float[_dSize + 1];

	for(int i = 0; i < _dSize; ++i)
	{
		_data[i] = 0;
	}
	_data[_dSize] = 1;
}

WeightMatrix::~WeightMatrix()
{
	delete[] _data;
}

void WeightMatrix::read(istream& is)
{
	if(_n > 0)
	{
		delete[] _data;
	}

	is>>_n;
	init();

	for(int i = 0; i < _dSize; ++i)
	{
		is>>_data[i];
	}
}

void WeightMatrix::write(ostream& os)
{
	os<<_n<<"\n";

	for(int i = 0; i < _dSize; ++i)
	{
		os<<_data[i]<<"\n";
	}
}

void WeightMatrix::print(ostream& os)
{
	int prevPrec = os.precision(3);
	for(int i = 0; i < _n; ++i)
	{
		for(int j = 0; j < _n; ++j)
		{
			os.width(4);
			os<<(*this)(i, j)<<" ";
		}
		os<<"\n";
	}
	os.precision(prevPrec);
}

int WeightMatrix::nodes()
{
	return _n;
}

int WeightMatrix::edges()
{
	return _dSize;
}

int WeightMatrix::index(int row, int col)
{
	assert(row < _n && col < _n);
	if(row > col)
	{
		return index(col, row);
	}
	else if(row == col)
	{
		return _dSize; //is 1.0
	}
	else
	{
		int idx = (_n - 1) * row - (row * (row - 1))/2 + col - (row + 1);
		assert(idx < _dSize);
		return idx;
	}
}

float& WeightMatrix::operator()(int row, int col)
{
	return _data[index(row, col)];
}

float& WeightMatrix::operator[](int idx)
{
	assert(idx >= 0 && idx < _dSize);
	return _data[idx];
}

float WeightMatrix::wplus(int row, int col)
{
	float wt = (*this)(row, col);
	switch(_mode)
	{
		case W_UNALTERED:
			return wt;
		case W_LOG:
			if(wt < EPSILON)
			{
				return log(EPSILON);
			}
			return log(wt);
	};

	//can't happen
	assert(0);
}

float WeightMatrix::wminus(int row, int col)
{
	float wt = (*this)(row, col);
	switch(_mode)
	{
		case W_UNALTERED:
			return 1 - wt;
		case W_LOG:
			if(wt > (1.0 - EPSILON))
			{
				return log(EPSILON);
			}
			return log(1 - wt);
	};

	//can't happen
	assert(0);
}

float WeightMatrix::wdiff(int row, int col)
{
	return wplus(row, col) - wminus(row, col);
}

int WeightMatrix::nodeLabels(ints& res)
{
	int nextLabel = 1;
	res.resize(nodes());

	for(int i = 0; i < _n; ++i)
	{
		if(res[i] == 0)
		{
			res[i] = nextLabel++;
		}

		for(int j = i + 1; j < _n; ++j)
		{
			float cwt = (*this)(i, j);
			if(cwt == 1) //nodes are connected
			{
				if(res[j] == 0) //so far unlabeled
				{
					res[j] = res[i];
				}
				else if(res[j] != res[i]) //triangle ineq. violated
				{
					res.clear();
					return -1;
				}
			}
			else
			{
				if(cwt != 0) //matrix is not 0-1
				{
					res.clear();
					return -1;
				}
			}
		}
	}

	return nextLabel - 1;
}

WeightMatrix& WeightMatrix::operator=(WeightMatrix & other)
{
	_n = other._n;
	_mode = other._mode;
	init();
	for(int i = 0; i < _dSize; i++)
	{
		_data[i] = other._data[i];
	}
	return *this;
}
