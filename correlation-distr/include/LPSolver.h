/* This class interfaces with CPLEX to solve the LP relaxation of correlation clustering. Some lines were adapted from the sample files in that package.

There are actually 4 related solvers here. The simplest, LPSolver, simply solves the LP by passing it to CPLEX, including all O(n^3) constraints. This scales to around 200 nodes. The ILPSolver does the same but with the variables restricted to be integral. The LazyLPSolver does not add the triangle inequalities initially, but rather repeatedly solves the LP and adds violated triangle inequalities in batches and resolves until none are violated. It ends up needing to tell CPLEX about O(n^2) of the O(n^3) constraints, dramatically improving scalability to around 500 nodes. The LazyILPSolver handles the lazy constraints in a completely unrelated way, using CPLEX's addLazyConstraint function.

If the clusters are large relative to the noise level the LP is often exact and hence the ILP scales as well (and as poorly) as the LP. If the clusters are small relative to the noise level then the LP is frequently a very poor relaxation and sets all X_ij greedily to either 1/2 or 0. For these instances the ILP does horribly.

Scalability of LazyILP??
*/

class LPSolver;
class ILPSolver;

#ifndef LP_SOLVER_H
#define LP_SOLVER_H

#include "WeightMatrix.h"

//limits.h is required for some of the CPLEX includes to work right apparently
#include <limits.h>
#include <ilcplex/ilocplex.h>
#include <ilconcert/ilomodel.h>

class LPSolver
{
  //the model implemented here is equivalent to the ACN '05 LP (figure 1)
  //it uses a single xij variable instead of splitting it into xij+ and xij-
  //and it solves max (wplus - wminus)xij instead of
  //min (wminus - wplus)xij

 public:
  LPSolver(WeightMatrix& wt);
  virtual ~LPSolver();

  virtual WeightMatrix& solve();

 protected:
  virtual void modelSetup(IloEnv& env, IloModel& model, IloNumVarArray& xij);
  virtual void solverSetup(IloEnv& env, IloCplex& cplex, 
			   IloNumVarArray& xij);
  virtual void recordResult(IloNumArray& vals);
  virtual void newVar(IloEnv& env, IloNumVarArray& addTo);

  WeightMatrix& _input;
  WeightMatrix _result;
};

class LazyLPSolver : public LPSolver
{
  //this model is the same as the above but adds the triangle constraints
  //lazily (by solving the model, then modifying it, until optimality)
  //the parameter addLazyFraction controls the number of violated constraints added at a time;
  //it adds about addLazyFraction * n * n violated constraints at a time.
  //A value of around 0.5-0.6 for addLazyFraction seems to work well.
 public:
  LazyLPSolver(WeightMatrix& wt, float addLazyFraction);

  virtual WeightMatrix& solve();

 protected:
  virtual void modelSetup(IloEnv& env, IloModel& model, IloNumVarArray& xij);
  virtual bool addConstraints(IloEnv& env, IloModel& model,
			      IloNumVarArray& xij, IloNumArray& vals);

  float _addLazyFraction;
};

class ILPSolver : public LPSolver
{
  //this model is the same as the above but with integer constraints

 public:
  ILPSolver(WeightMatrix& wt);

 protected:
  virtual void recordResult(IloNumArray& vals);
  virtual void newVar(IloEnv& env, IloNumVarArray& xij);
};

class LazyILPSolver : public ILPSolver
{
  //this model is the same as the above but adds the triangle constraints
  //lazily (using Cplex's addLazyConstraint function)

 public:
  LazyILPSolver(WeightMatrix& wt);

 protected:
  virtual void modelSetup(IloEnv& env, IloModel& model, IloNumVarArray& xij);
  virtual void solverSetup(IloEnv& env, IloCplex& cplex, 
			   IloNumVarArray& xij);
};

#endif
