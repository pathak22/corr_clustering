/* This class interfaces with the SDP Solver based on Christoph Helmberg's Conic Bundle package. Some lines were adapted from the sample files in that package.

This solver iteratively finds better dual feasible solutions. Each of these provides a valid lower bound despite the fact that it hasn't converged and the primal X is often very infeasible.
*/

class SDPSolverTwo;

#ifndef SDP_SOLVER_TWO_H
#define SDP_SOLVER_TWO_H

#include "WeightMatrix.h"

//Load various Conic Bundle headers
#include "MatFCBSolver.hxx"
#include "MatSDPfun.hxx"
#include "cmsymden.hxx"
#include "cmsingle.hxx"

//structure for storing triples of integers identifying a triangle constraint
struct Triple {
  int _i, _j, _k;
  Triple(int i, int j, int k) : _i(i), _j(j), _k(k) {
    assert(i<j && i != k && j != k);
  }
  //Lexicographic comparison for std::set data structure
  bool operator<(const Triple & other) const {
    if(_i < other._i) 
      return true;
    else if (_i > other._i)
      return false;
    else {
      
      if(_j < other._j) 
	return true;
      else if (_j > other._j)
	return false;
      else {
	if(_k < other._k) 
	  return true;
	else if (_k > other._k)
	  return false;
	else 
	  return false; //equal
      }
    }
  }
};

class SDPSolverTwo
{
 public:
  //If the doTriangles parameter is true, it will add triangle inequalities that are violated in attempt to improve the relaxation. This is quite helpful for instances with at most a few hundred nodes, but of questionable value for more nodes since even the basic relaxation takes so long to converge.
  //If the lazy parameter is true, it adds the X_ij>=0 constraints lazily. Setting imrpoves the relaxation value found after a few minutes a lot, but delays eventual convergence slightly.
  SDPSolverTwo(WeightMatrix& wt, bool doTriangles=false, bool lazy=false);
  virtual ~SDPSolverTwo();
  
  virtual WeightMatrix& solve();
  
 protected:
  void AddConstraint(int i, int j); //queue a X_ij>=0 or X_ii = 1 constraint for addition to the SDP soon
  void AddTriangleConstraint(int i, int j, int k); //queue a triangle constraint
  void ReallyAddConstraints(); //Tell the SDP solver about the queued constraints

  WeightMatrix& _input; //the classifier outcomes.
  WeightMatrix _result; //solve() returns a reference to this 
  int VERBOSE;
  int _nodes;
  ConicBundle::MatrixSDPfunction _sdp;
  ConicBundle::MatrixCBSolver _cbsolver;
  CH_Matrix_Classes::Symmatrix _boundVars; //which variables have we added the X_ij>=0 or X_ii = 1 constraints for?
  std::set<Triple> _triangleCons; //which triangle inequalities have we added?

  int _nBoundVars; //Number of non-zeros in lower triangle of_boundVars (includes diagonal)

  //next 5 variables are for constraints that are in the process of being added
  int _numNewCon;
  std::map<CH_Matrix_Classes::Integer,ConicBundle::SparseSDPCoeffVector> _newOpA; //map from constraint number, block number to constraint matrix
  CH_Matrix_Classes::Matrix _newRHSides;
  CH_Matrix_Classes::Matrix _newLBs;
  CH_Matrix_Classes::Matrix _newUBs;

  bool _doTriangles; //add triangle constraints?
  bool _lazy; //add the X_ij >= 0 constraints lazily?
};

#endif
