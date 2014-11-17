#include "LPSolver.h"

LPSolver::LPSolver(WeightMatrix& wt):
	_input(wt),
	_result(wt.nodes(), W_UNALTERED)
{
}

LPSolver::~LPSolver()
{
}

WeightMatrix& LPSolver::solve()
{
	IloEnv env;

	try
	{
		IloModel model(env);
		//one variable per edge
		IloNumVarArray xij(env);

		modelSetup(env, model, xij);

		IloCplex cplex(model);

		solverSetup(env, cplex, xij);

		int status = cplex.solve();

		if(!status)
		{
			cerr<<"Failed to optimize LP.\n";
			cerr<<"CPLEX status "<<cplex.getCplexStatus()<<"\n";
			cerr<<"Status "<<cplex.getStatus()<<"\n";

			return _result;
		}

		if(cplex.getCplexStatus() != CPX_STAT_OPTIMAL)
		{
			cerr<<"Failed to find an optimal solution.\n";
			cerr<<"CPLEX status "<<cplex.getCplexStatus()<<"\n";
			cerr<<"Status "<<cplex.getStatus()<<"\n";

			return _result;
		}

		cerr<<"Status "<<cplex.getStatus()<<"\n";

		IloNumArray vals(env);
		cplex.getValues(vals, xij);

		recordResult(vals);
	}
	catch(IloException& e)
	{
		cerr<<"Concert exception caught: "<<e<<"\n";
	}
	catch(...) 
	{
		cerr<<"Unknown exception caught.\n";
	}

	env.end();

	return _result;
}

void LPSolver::modelSetup(IloEnv& env, IloModel& model, IloNumVarArray& xij)
{
 	IloExpr objective(env);

	int ct = 0;
	int n = _input.nodes();
	for(int ii = 0; ii < n; ++ii)
	{
		for(int jj = ii + 1; jj < n; ++jj,++ct)
		{
			newVar(env, xij);

			objective += xij[ct] * _input.wdiff(ii, jj);
			assert(_input.index(ii, jj) == ct);

			//cerr<<"cost "<<ii<<","<<jj<<" = "<<_input.wdiff(ii, jj)<<"\n";
		}
	}

	model.add(IloMaximize(env, objective));

	int nct = 0;
	for(int ii = 0; ii < n; ++ii)
	{
		for(int jj = 0; jj < n; ++jj)
		{
			if(ii == jj)
			{
				continue;
			}

			for(int kk = ii + 1; kk < n; ++kk)
			{
				if(kk == jj)
				{
					continue;
				}

				++nct;

// 				cerr<<": "<<ii<<","<<jj<<" and "<<jj<<","<<kk<<" ==> "<<
// 					ii<<","<<kk<<"\n";

				IloExpr triangle(env);
				//constrain edge ik
				triangle += xij[_input.index(ii,jj)]
					+ xij[_input.index(jj, kk)]
					- xij[_input.index(ii, kk)];
				model.add(triangle <= 1);
			}
		}
	}

	cerr<<nct<<" constraints generated\n";
}

void LPSolver::recordResult(IloNumArray& vals)
{
//	cerr<<"values "<<vals<<"\n";

	int n = _result.edges();
	for(int ii = 0; ii < n; ++ii)
	{
		_result[ii] = vals[ii];
	}
}

void LPSolver::newVar(IloEnv& env, IloNumVarArray& addTo)
{
	addTo.add(IloNumVar(env, 0.0, 1.0));
}

void LPSolver::solverSetup(IloEnv& env, IloCplex& cplex,
						   IloNumVarArray& xij)
{
}

ILPSolver::ILPSolver(WeightMatrix& input):
	LPSolver(input)
{
}

void ILPSolver::newVar(IloEnv& env, IloNumVarArray& addTo)
{
	addTo.add(IloNumVar(env, 0.0, 1.0, ILOINT));
}

void ILPSolver::recordResult(IloNumArray& vals)
{
//	cerr<<"values "<<vals<<"\n";

	int n = _result.edges();
	for(int ii = 0; ii < n; ++ii)
	{
		if(approxEq(vals[ii], 1))
		{
			_result[ii] = 1;
		}
		else if(approxEq(vals[ii], 0))
		{
			_result[ii] = 0;
		}
		else
		{
			cerr<<"Non-integral solution to an ILP: "<<vals[ii]<<"\n";
			assert(0);
		}
	}
}

LazyLPSolver::LazyLPSolver(WeightMatrix& wt, float addLazyFraction):
	LPSolver(wt),
	_addLazyFraction(addLazyFraction)
{
}

WeightMatrix& LazyLPSolver::solve()
{
	IloEnv env;

	try
	{
		IloModel model(env);
		//one variable per edge
		IloNumVarArray xij(env);

		modelSetup(env, model, xij);

		IloCplex cplex(model);

		solverSetup(env, cplex, xij);

		bool done = false;

		while(!done)
		{
			int status = cplex.solve();

			if(!status)
			{
				cerr<<"Failed to optimize LP.\n";
				cerr<<"CPLEX status "<<cplex.getCplexStatus()<<"\n";
				cerr<<"Status "<<cplex.getStatus()<<"\n";

				return _result;
			}

			if(cplex.getCplexStatus() != CPX_STAT_OPTIMAL)
			{
				cerr<<"Failed to find an optimal solution.\n";
				cerr<<"CPLEX status "<<cplex.getCplexStatus()<<"\n";
				cerr<<"Status "<<cplex.getStatus()<<"\n";

				return _result;
			}

			cerr<<"Status "<<cplex.getStatus()<<"\n";

			IloNumArray vals(env);
			cplex.getValues(vals, xij);

			done = addConstraints(env, model, xij, vals);
		}

		IloNumArray vals(env);
		cplex.getValues(vals, xij);
		recordResult(vals);
	}
	catch(IloException& e)
	{
		cerr<<"Concert exception caught: "<<e<<"\n";
	}
	catch(...) 
	{
		cerr<<"Unknown exception caught.\n";
	}

	env.end();

	return _result;
}

void LazyLPSolver::modelSetup(IloEnv& env, IloModel& model, 
							  IloNumVarArray& xij)
{
 	IloExpr objective(env);

	int ct = 0;
	int n = _input.nodes();
	for(int ii = 0; ii < n; ++ii)
	{
		for(int jj = ii + 1; jj < n; ++jj,++ct)
		{
			newVar(env, xij);

			objective += xij[ct] * _input.wdiff(ii, jj);
			assert(_input.index(ii, jj) == ct);

			//cerr<<"cost "<<ii<<","<<jj<<" = "<<_input.wdiff(ii, jj)<<"\n";
		}
	}

	model.add(IloMaximize(env, objective));
}

bool LazyLPSolver::addConstraints(IloEnv& env, IloModel& model,
								  IloNumVarArray& xij, IloNumArray& vals)
{
	int n = _input.nodes();

	int nct = 0;
	int violated = 0;
	int added = 0;
	//first count the number of violated constraints
	for(int ii = 0; ii < n; ++ii)
	{
		for(int jj = 0; jj < n; ++jj)
		{
			if(ii == jj)
			{
				continue;
			}

			for(int kk = ii + 1; kk < n; ++kk)
			{
				if(kk == jj)
				{
					continue;
				}

				++nct;

				double triangle = vals[_input.index(ii,jj)]
					+ vals[_input.index(jj, kk)]
					- vals[_input.index(ii, kk)];

				if(triangle - 1.0 > 1e-5)
				{
					++violated;
				}
			}
		}
	}
	//expected number of constraints to add: _addLazyFraction * n * n
	//If _addLazyFraction = 0.5, the number of constraints added is equal to the number of variables (size of basis), which makes sense.
	double sampleFraction = fmin(1., _addLazyFraction * double(n * n) / violated);
	//now add violated constraints
	for(int ii = 0; ii < n; ++ii)
	{
		for(int jj = 0; jj < n; ++jj)
		{
			if(ii == jj)
			{
				continue;
			}

			for(int kk = ii + 1; kk < n; ++kk)
			{
				if(kk == jj)
				{
					continue;
				}

				double triangle = vals[_input.index(ii,jj)]
					+ vals[_input.index(jj, kk)]
					- vals[_input.index(ii, kk)];

				if(triangle - 1.0 > 1e-5)
				{

// 					cerr<<": "<<ii<<","<<jj<<" and "<<jj<<","<<kk<<" ==> "<<
// 						ii<<","<<kk<<", currently "<<triangle<<"\n";
// 					cerr<<"check "<<(triangle > 1)<<"\n";

					if(violated == 1 || rand() <= sampleFraction * RAND_MAX)
					{
						IloExpr triangleConstraint(env);
						//constrain edge ik
						triangleConstraint += xij[_input.index(ii,jj)]
							+ xij[_input.index(jj, kk)]
							- xij[_input.index(ii, kk)];
						model.add(triangleConstraint <= 1);

						++added;
					}
				}
			}
		}
	}
	cerr<<nct<<" lazy constraints checked, "<<violated<<" violated, "<<
	  added<<" generated, i.e. about n^2*" << _addLazyFraction << "="<<_addLazyFraction * n*n << std::endl;
	return (violated == 0);
}

LazyILPSolver::LazyILPSolver(WeightMatrix& wt):
	ILPSolver(wt)
{
}

void LazyILPSolver::modelSetup(IloEnv& env, IloModel& model, 
							  IloNumVarArray& xij)
{
 	IloExpr objective(env);

	int ct = 0;
	int n = _input.nodes();
	for(int ii = 0; ii < n; ++ii)
	{
		for(int jj = ii + 1; jj < n; ++jj,++ct)
		{
			newVar(env, xij);

			objective += xij[ct] * _input.wdiff(ii, jj);
			assert(_input.index(ii, jj) == ct);

			//cerr<<"cost "<<ii<<","<<jj<<" = "<<_input.wdiff(ii, jj)<<"\n";
		}
	}

	model.add(IloMaximize(env, objective));
}

void LazyILPSolver::solverSetup(IloEnv& env, IloCplex& cplex, 
							   IloNumVarArray& xij)
{
	int n = _input.nodes();

	int nct = 0;
	for(int ii = 0; ii < n; ++ii)
	{
		for(int jj = 0; jj < n; ++jj)
		{
			if(ii == jj)
			{
				continue;
			}

			for(int kk = ii + 1; kk < n; ++kk)
			{
				if(kk == jj)
				{
					continue;
				}

				++nct;

// 				cerr<<": "<<ii<<","<<jj<<" and "<<jj<<","<<kk<<" ==> "<<
// 					ii<<","<<kk<<"\n";

				IloExpr triangle(env);
				//constrain edge ik
				triangle += xij[_input.index(ii,jj)]
					+ xij[_input.index(jj, kk)]
					- xij[_input.index(ii, kk)];
				cplex.addLazyConstraint(triangle <= 1);
			}
		}
	}

	cerr<<nct<<" lazy constraints generated\n";
}
