#include "SDPSolver.h"

SDPSolver::SDPSolver(WeightMatrix& wt, bool addNonNegC):
	_input(wt),
	_result(wt.nodes(), W_UNALTERED),
	VERBOSE(1),
	_addNonNegC(addNonNegC)
{
}

SDPSolver::~SDPSolver()
{
}

WeightMatrix& SDPSolver::solve()
{
  bool nonNegC = _addNonNegC;
	DSDP dsdp;
	int error;

	int nodes = _input.nodes();
	int edges = _input.edges();

	error = DSDPCreate(nodes + (nonNegC ? edges : 0), &dsdp);
	assert(!error);

	SDPCone sdpcone;
	//number of blocks: 1 block for the matrix, 1 additional per
	//(edgewise) slack var
	int nblocks = nonNegC ? edges + 1 : 1;
	error = DSDPCreateSDPCone(dsdp, nblocks, &sdpcone);
	assert(!error);

	error = SDPConeSetBlockSize(sdpcone, 0, nodes);
	assert(!error);

	if (nonNegC) {
	  for(int ii = 1; ii < nblocks; ++ii)
	    {
	      error = SDPConeSetBlockSize(sdpcone, ii, 1);
	      assert(!error);
	    }
	}

	mainProgramSetup(dsdp, sdpcone);
	if (nonNegC)
	  addNonNegConstraints(dsdp, sdpcone);
	
	bool worked = runDSDP(dsdp);
	assert(worked);
	readback(dsdp, sdpcone);
	freeAll(dsdp);

	return _result;
}

void SDPSolver::writeFile(ostream& os)
{
	os<<"*File description given at "<<
		"http://infohost.nmt.edu/~sdplib/FORMAT\n";

	int nodes = _input.nodes();
	int edges = _input.edges();

	//one Y_ii constraint per node, one Y_ij constraint per pair
	//of distinct nodes
	int nMatrices = nodes + edges;

	os<<nMatrices<<"\n";

	//one n*n block for the costs, and 1-cell blocks for all else
	os<<(1 + edges)<<"\n";

	//block sizes
	os<<nodes;
	for(int ii = 0; ii < edges; ++ii)
	{
		os<<" "<<1;
	}
	os<<"\n";

	//objective terms: the first n variables are real, and the
	//remaining ones are slack variables
	for(int ii = 0; ii < nodes; ++ii)
	{
		os<<1<<" ";
	}

	for(int ii = 0; ii < edges; ++ii)
	{
		os<<0<<" ";
	}
	os<<"\n";

	//constraint entries
	//matrix, block, i, j, val

	//F0: the cost matrix
	for(int ii = 0; ii < nodes; ++ii)
	{
		for(int jj = ii + 1; jj < nodes; ++jj)
		{
			os<<0<<" "<<1<<" "<<(ii + 1)<<" "<<(jj + 1)<<
				" "<<_input.wdiff(ii, jj)<<"\n";
		}
	}

	//F1-N: Y_ii constraints
	for(int ii = 0; ii < nodes; ++ii)
	{
		os<<(ii + 1)<<" "<<1<<" "<<(ii + 1)<<" "<<(ii + 1)<<" "<<1<<"\n";
	}

	//FN+1-..: Y_ij constraints
	int ct = 0;
	for(int ii = 0; ii < nodes; ++ii)
	{
		for(int jj = ii + 1; jj < nodes; ++jj,++ct)
		{
			os<<(ct + nodes + 1)<<
				" "<<1<<" "<<(ii + 1)<<" "<<(jj + 1)<<" "<<1<<"\n";
			os<<(ct + nodes + 1)<<
				" "<<(ct + 2)<<" "<<1<<" "<<1<<" "<<-1<<"\n";
		}
	}
}

void SDPSolver::mainProgramSetup(DSDP& dsdp, SDPCone& sdpcone)
{
	int nodes = _input.nodes();
	int edges = _input.edges();

	int error;

	//self-edges must equal 1; also sets up the objective
	_diag = new double[nodes];
	_iptr = new int[nodes];

	for(int ii = 0; ii < nodes; ++ii)
	{
		_iptr[ii] = ii * (ii + 1)/2 + ii;
		_diag[ii] = 1.0;
	}

	for(int ii = 0; ii < nodes; ++ii)
	{
		error = DSDPSetDualObjective(dsdp, ii + 1, 1.0);
		assert(!error);

		if(VERBOSE > 10)
		{
			cerr<<"ii "<<ii<<" "<<_iptr[ii]<<" "<<_diag[ii]<<"\n";
		}

		error = SDPConeSetASparseVecMat(sdpcone,
										0, //block
										ii + 1, //var
										nodes, //array size
										1.0, //alpha (multiplies value)
										0, //ishift (scary magic)
										_iptr + ii, //index array
										_diag + ii, //value array
										1 //nonzero elements in arrays
			);
		assert(!error);

		if(VERBOSE > 5)
		{
			cerr<<"Matrix: "<<(ii + 1)<<"\n";
			error = SDPConeViewDataMatrix(sdpcone,0,ii+1);
		}
	}
	if(VERBOSE > 10)
	{
		error = SDPConeView(sdpcone);
	}

	//actual weight matrix
	int nWeights = nodes + edges;
	_weightMat = new double[nWeights];

	int ct = 0;
	for(int ii = 0; ii < nodes; ++ii)
	{
		for(int jj = 0; jj <= ii; ++jj,++ct)
		{
			if(ii == jj)
			{
				_weightMat[ct] = 0;
			}
			else
			{
				//for some reason this needs to be negated!
				_weightMat[ct] = - _input.wdiff(ii, jj);
			}
			//cerr<<ii<<" "<<jj<<" "<<ct<<" "<<weightMat[ct]<<"\n";
		}
	}

	error = SDPConeSetADenseVecMat(sdpcone,
								   0, //block
								   0, //var 0 is the cost matrix?
								   nodes, //array size
								   1.0, //alpha (multiplier)
								   _weightMat, //data
								   nWeights //data array length
		);
	assert(!error);

	if(VERBOSE > 5)
	{
		error = SDPConeViewDataMatrix(sdpcone,0,0);
	}
	if(VERBOSE > 10)
	{
		error = SDPConeView(sdpcone);
	}
}

void SDPSolver::addNonNegConstraints(DSDP& dsdp, SDPCone& sdpcone)
{
	int nodes = _input.nodes();

	int error;

	//edge non-negativity
	_nonNeg = new double[1];
	_nonNeg[0] = 1;
	_neg = new double[1];
	_neg[0] = -1;
	_zeroIndex = new int[1];
	_zeroIndex[0] = 0;

	int varCt = 0;
	int edgeCt = 0;
	for(int ii = 0; ii < nodes; ++ii)
	{
		for(int jj = 0; jj <= ii; ++jj,++edgeCt)
		{
			if(ii == jj)
			{
				continue;
			}

			++varCt;

			int* nonNegInd = new int[1];
			nonNegInd[0] = edgeCt;
			_nonNegInds.push_back(nonNegInd);

			int varNum = varCt + nodes;
			int blockNum = varCt;

			if(VERBOSE > 5)
			{
				cerr<<"setting constraint for variable "<<varNum<<
					" block "<<blockNum<<
					" at matrix-array index "<<edgeCt<<"\n";
			}

			error = SDPConeSetASparseVecMat(sdpcone,
											0, //block
											varNum, //var
											nodes, //block size?
											1.0, //alpha (multiplies value)
											0, //ishift (scary magic)
											nonNegInd, //index array
											_nonNeg, //value array
											1 //nonzero elements in arrays
				);
			assert(!error);
			if(VERBOSE > 5)
			{
				error = SDPConeViewDataMatrix(sdpcone,0,varNum);
			}

			error = SDPConeSetASparseVecMat(sdpcone,
											blockNum, //block
											varNum, //var
											1, //block size?
											1.0, //alpha (multiplies value)
											0, //ishift (scary magic)
											_zeroIndex, //index array
											_neg, //value array
											1 //nonzero elements in arrays
				);
			assert(!error);
			if(VERBOSE > 5)
			{
				error = SDPConeViewDataMatrix(sdpcone,blockNum,varNum);
			}
		}
	}
	if(VERBOSE > 10)
	{
		error = SDPConeView(sdpcone);
	}
}

bool SDPSolver::runDSDP(DSDP& dsdp)
{
	int error;

	//actually solve
	error = DSDPSetup(dsdp);
	assert(!error);
	DSDPSetStandardMonitor(dsdp, 1);
	assert(!error);

	
	////BEGINNING OF SDP TUNING PARAMETERS

	int n = _input.nodes();
	//choose feasible y0 (default is 0)
	//the default y0=0 seems to work as well as anything else I can come up with, so commented out
	//for(int i = 0; i < _input.edges(); i++) {
	//  error=DSDPSetY0(dsdp, i+1+n , sum / _input.edges());
	//  assert(!error);
	//}

	//set r to make dual feasible
	error=DSDPSetR0(dsdp, n);
	assert(!error);

	//set bounds on dual vars y 
	error=DSDPSetYBounds(dsdp, -3*n, 3*n);
	assert(!error);
       
	//gamma=2*n*n , a bound on trace. The node vars have trace n and the slack variables have trace at most n(n-1). 
		error=DSDPSetPenaltyParameter(dsdp, 2*n*n);
	assert(!error);

	//upper-bound zero works
	error=DSDPSetZBar(dsdp, 0);
	assert(!error);
	
	//As noted in documentation 15 is a good value for SDPs like this one with lots of constraints. Experiments confirm that.
	error=DSDPReuseMatrix(dsdp, 15);
	assert(!error);

	//default works fine (this barrier param is aka nu)
	//error=DSDPSetBarrierParameter(dsdp, 1);
	//assert(!error);


	//10 works ~4% faster than default of 4
	error=DSDPSetPotentialParameter(dsdp, 10);
	assert(!error);

	//The default is better.
	//error=DSDPUsePenalty(dsdp, 1);

	/////END OF SDP TUNING PARAMETERS

	DSDPSolve(dsdp);
	assert(!error);

	DSDPSolutionType feasible;
	error = DSDPGetSolutionType(dsdp, &feasible);
	assert(!error);

	if(feasible == DSDP_UNBOUNDED)
	{
		cerr<<"DSDP: Dual Unbounded, Primal Infeasible\n";
		return false;
	}
	else if(feasible == DSDP_INFEASIBLE)
	{
		cerr<<"DSDP: Primal Unbounded, Dual Infeasible\n";
		return false;
	}
	else
	{
		cerr<<"DSDP: normal termination\n";
		return true;
	}
}

void SDPSolver::readback(DSDP& dsdp, SDPCone& sdpcone)
{
	int nodes = _input.nodes();

	int error;

	error = DSDPComputeX(dsdp);
	assert(!error);

	//view the solution and stuff
	double* xmat;
	int nsize;
	error = SDPConeGetXArray(sdpcone, 0, &xmat, &nsize);
	assert(!error);

	if(VERBOSE > 10)
	{
		SDPConeViewX(sdpcone, 0, nodes, xmat, nsize);
	}

	//read back the output variables
	int ct = 0;
	for(int ii = 0; ii < nodes; ++ii)
	{
		for(int jj = 0; jj <= ii; ++jj,++ct)
		{
			if(ii != jj)
			{
				_result(ii, jj) = xmat[ct];
			}
		}
	}
}

void SDPSolver::freeAll(DSDP& dsdp)
{
	//free all the resources
	delete[] _diag;
	delete[] _iptr;
	delete[] _weightMat;
	delete[] _nonNeg;
	delete[] _neg;
	delete[] _zeroIndex;

	foreach(intPtrs, ptr, _nonNegInds)
	{
		delete[] *ptr;
	}
	_nonNegInds.clear();

	DSDPDestroy(dsdp);
}

LazySDPSolver::LazySDPSolver(WeightMatrix& input, float addLazyFraction,
	int iters):
	SDPSolver(input),
	_addLazyFraction(addLazyFraction),
	_iters(iters)
{
}

WeightMatrix& LazySDPSolver::solve()
{
	int error;

	int nodes = _input.nodes();
	int edges = _input.edges();

	int violated = 1; //always do initial iteration

	int iter = 0;

	while(violated && iter != _iters)
	{
		if(VERBOSE)
		{
			cerr<<"iteration "<<iter<<" of "<<_iters<<"\n";
		}

		++iter;

		int nconstraints = _constraints.size();
		assert(nconstraints <= edges);

		DSDP dsdp;

		error = DSDPCreate(nodes + nconstraints, &dsdp);
		assert(!error);

		//number of blocks: 1 block for the matrix, 1 additional per
		//(edgewise) slack var
		int nblocks = nconstraints + 1;
		SDPCone sdpcone;
		error = DSDPCreateSDPCone(dsdp, nblocks, &sdpcone);
		assert(!error);

		error = SDPConeSetBlockSize(sdpcone, 0, nodes);
		assert(!error);

		for(int ii = 1; ii < nblocks; ++ii)
		{
			error = SDPConeSetBlockSize(sdpcone, ii, 1);
			assert(!error);
		}

		mainProgramSetup(dsdp, sdpcone);
		addNonNegConstraints(dsdp, sdpcone);
		bool worked = runDSDP(dsdp);
		assert(worked);
		readback(dsdp, sdpcone);
		violated = createConstraints();

		freeAll(dsdp);

		/*
		  //uncomment to print out intermediate matrices to files
		std::string s("lazysdp_");
		s += ('a' + iter);
		s += ".mat";
		std::ofstream outFile(s.c_str());
		_result.write(outFile);
		*/
	}
	return _result;
}

int LazySDPSolver::matrixInd(int row, int col)
{
	//array format is the lower triangle, including the diagonal
	return (row * (row + 1))/2 + col;
}

void LazySDPSolver::addNonNegConstraints(DSDP& dsdp, SDPCone& sdpcone)
{
	int nodes = _input.nodes();

	int error;

	//edge non-negativity
	_nonNeg = new double[1];
	_nonNeg[0] = 1;
	_neg = new double[1];
	_neg[0] = -1;
	_zeroIndex = new int[1];
	_zeroIndex[0] = 0;

	int varCt = 0;
	foreach(Constraints, constraint, _constraints)
	{
		++varCt;

		int* nonNegInd = new int[1];
		nonNegInd[0] = matrixInd(constraint->row, constraint->col);
		_nonNegInds.push_back(nonNegInd);

		int varNum = varCt + nodes;
		int blockNum = varCt;

		if(VERBOSE > 5)
		{
			cerr<<"setting constraint for variable "<<varNum<<
				" block "<<blockNum<<
				" at matrix-array index "<<nonNegInd[0]<<"\n";
		}

		error = 
			SDPConeSetASparseVecMat(sdpcone,
									0, //block
									varNum, //var
									nodes, //block size?
									1.0, //alpha (multiplies value)
									0, //ishift (scary magic)
									nonNegInd, //index array
									_nonNeg, //value array
									1 //nonzero elements in arrays
				);
		assert(!error);
		if(VERBOSE > 5)
		{
			error = SDPConeViewDataMatrix(sdpcone,0,varNum);
		}

		error = 
			SDPConeSetASparseVecMat(sdpcone,
									blockNum, //block
									varNum, //var
									1, //block size?
									1.0, //alpha (multiplies value)
									0, //ishift (scary magic)
									_zeroIndex, //index array
									_neg, //value array
									1 //nonzero elements in arrays
				);
		assert(!error);
		if(VERBOSE > 5)
		{
			error = SDPConeViewDataMatrix(sdpcone,blockNum,varNum);
		}
	}
}

int LazySDPSolver::createConstraints()
{
	int violated = 0;
	int created = 0;

	int nodes = _input.nodes();

	for(int ii = 0; ii < nodes; ++ii)
	{
		for(int jj = 0; jj <= ii; ++jj)
		{
			if(_result(ii, jj) < 0)
			{
				++violated;

				if(created == 0 || rand() <= _addLazyFraction * RAND_MAX)
				{
					++created;
					_constraints.push_back(Constraint(ii, jj));
				}
			}
		}
	}

	if(VERBOSE)
	{
		cerr<<violated<<" constraints violated, "<<created<<" generated, "<<
			_constraints.size()<<" total.\n";
	}

	return violated;
}

SumLazySDPSolver::SumLazySDPSolver(WeightMatrix& input, int iters):
	LazySDPSolver(input, 1.0, iters)
{
}

WeightMatrix& SumLazySDPSolver::solve()
{
	//really this block should be inherited but there
	//are some stupid little changes that I don't want to factor right now

	int error;

	int nodes = _input.nodes();

	int violated = 1; //always do initial iteration

	int iter = 0;

	while(violated && iter != _iters)
	{
		if(VERBOSE)
		{
			cerr<<"iteration "<<iter<<" of "<<_iters<<"\n";
		}

		++iter;

		int nconstraints = _sumConstraints.size();

		DSDP dsdp;

		error = DSDPCreate(nodes + nconstraints, &dsdp);
		assert(!error);

		//number of blocks: 1 block for the matrix, 1 additional per
		//(edgewise) slack var
		int nblocks = nconstraints + 1;
		SDPCone sdpcone;
		error = DSDPCreateSDPCone(dsdp, nblocks, &sdpcone);
		assert(!error);

		error = SDPConeSetBlockSize(sdpcone, 0, nodes);
		assert(!error);

		for(int ii = 1; ii < nblocks; ++ii)
		{
			error = SDPConeSetBlockSize(sdpcone, ii, 1);
			assert(!error);
		}

		mainProgramSetup(dsdp, sdpcone);
		addNonNegConstraints(dsdp, sdpcone);
		bool worked = runDSDP(dsdp);
		assert(worked);

		//	DSDPEventLogSummary();

		readback(dsdp, sdpcone);
		violated = createConstraints();

		freeAll(dsdp);

		/* //uncomment this to print out intermediate matrices to files sumsdp_<letter>.mat
		std::string s("sumsdp_");
		s += ('a' + iter);
		s += ".mat";
		std::ofstream outFile(s.c_str());
		_result.write(outFile);
		*/
	}

	roundNegatives();

	return _result;
}

void SumLazySDPSolver::addNonNegConstraints(DSDP& dsdp, SDPCone& sdpcone)
{
	int nodes = _input.nodes();

	int error;

	//edge non-negativity
	_nonNeg = new double[1];
	_nonNeg[0] = 1;
	_neg = new double[1];
	_neg[0] = -1;
	_zeroIndex = new int[1];
	_zeroIndex[0] = 0;

	int varCt = 0;
	foreach(SumConstraints, constraint, _sumConstraints)
	{
		++varCt;

		int varNum = varCt + nodes;
		int blockNum = varCt;

		foreach(Constraints, edge, *constraint)
		{
			int* nonNegInd = new int[1];
			nonNegInd[0] = matrixInd(edge->row, edge->col);
			_nonNegInds.push_back(nonNegInd);

			if(VERBOSE > 5)
			{
				cerr<<"\t"<<edge->row<<" "<<edge->col<<"\n";
			}

			error = 
				SDPConeAddASparseVecMat(sdpcone,
										0, //block
										varNum, //var
										nodes, //block size?
										1.0, //alpha (multiplies value)
										0, //ishift (scary magic)
										nonNegInd, //index array
										_nonNeg, //value array
										1 //nonzero elements in arrays
					);
			assert(!error);
		}

		if(VERBOSE > 5)
		{
			cerr<<"\n";
		}

		error = 
			SDPConeSetASparseVecMat(sdpcone,
									blockNum, //block
									varNum, //var
									1, //block size?
									1.0, //alpha (multiplies value)
									0, //ishift (scary magic)
									_zeroIndex, //index array
									_neg, //value array
									1 //nonzero elements in arrays
				);
		assert(!error);
		if(VERBOSE > 5)
		{
			error = SDPConeViewDataMatrix(sdpcone,blockNum,varNum);
		}
	}
}

int SumLazySDPSolver::createConstraints()
{
	int violated = 0;
	int created = 0;

	int nodes = _input.nodes();

	Constraints newConstraint;

	for(int ii = 0; ii < nodes; ++ii)
	{
		for(int jj = 0; jj <= ii; ++jj)
		{
			if(_result(ii, jj) < 0)
			{
				++violated;
				++created;
				newConstraint.push_back(Constraint(ii, jj));
			}
		}
	}

	if(VERBOSE)
	{
		cerr<<violated<<" constraints violated, "<<created<<" generated, "<<
			_sumConstraints.size()<<" total.\n";
	}

	if(violated)
	{
		_sumConstraints.push_back(newConstraint);
	}

	return violated;
}


void SumLazySDPSolver::roundNegatives()
{
	int nodes = _input.nodes();

	for(int ii = 0; ii < nodes; ++ii)
	{
		for(int jj = 0; jj <= ii; ++jj)
		{
			if(_result(ii, jj) < 0)
			{
				_result(ii, jj) = 0;
			}
		}
	}
}

