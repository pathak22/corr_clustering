/* This class interfaces with the DSDP SDP solver. Some lines were adapted from the sample files in that package.

This solver does not produces valid lower bounds until it has solved the SDP to high accuracy. It is reasonably fast for up to around 100 nodes but scales extremely poorly to larger data sets. 
*/

class SDPSolver;
class LazySDPSolver;

#ifndef SDP_SOLVER_H
#define SDP_SOLVER_H

#include "WeightMatrix.h"

#include "dsdp5.h"

class SDPSolver
{
 public:
  //Setting parameter addNonNegC to false eliminates the non-negativity constraints. This improves scalability but at the expense of making the relaxation too loose to be useful.
  SDPSolver(WeightMatrix& wt, bool addNonNegC = true);
  virtual ~SDPSolver();

  //Write SDP description to a file for solving by the command-line DSDP solver. Useful for debugging.
  virtual void writeFile(ostream& os);
  virtual WeightMatrix& solve();

 protected:
  virtual void mainProgramSetup(DSDP& dsdp, SDPCone& sdpcone);
  virtual void addNonNegConstraints(DSDP& dsdp, SDPCone& sdpcone);
  virtual bool runDSDP(DSDP& dsdp);
  virtual void readback(DSDP& dsdp, SDPCone& sdpcone);
  virtual void freeAll(DSDP& dsdp);

  WeightMatrix& _input;
  WeightMatrix _result;
  int VERBOSE;
  double* _diag;
  int* _iptr;
  double* _weightMat;
  double* _neg;
  double* _nonNeg;
  int* _zeroIndex;
  intPtrs _nonNegInds;
 private:
  bool _addNonNegC; //add non-negative constraints?
};

struct Constraint
{
  Constraint(int r, int c):row(r),col(c){}
  int row;
  int col;
};

typedef std::vector<Constraint> Constraints;

//most basic lazy algorithm: iteratively add violated edge constraints
//(like the LP)
class LazySDPSolver : public SDPSolver
{
 public:
  LazySDPSolver(WeightMatrix& wt, float addLazyFraction, int iters);

  virtual WeightMatrix& solve();

 protected:
  virtual void addNonNegConstraints(DSDP& dsdp, SDPCone& sdpcone);
  virtual int createConstraints();
  virtual int matrixInd(int row, int col);

  float _addLazyFraction;
  Constraints _constraints;
  int _iters;
};

typedef std::vector<Constraints> SumConstraints;

//WS suggested lazy algorithm: iteratively add a sum constraint
//for all negative edges
class SumLazySDPSolver : public LazySDPSolver
{
 public:
  //use iteration count -1 for unbounded
  SumLazySDPSolver(WeightMatrix& wt, int iters);

  virtual WeightMatrix& solve();

 protected:
  virtual void addNonNegConstraints(DSDP& dsdp, SDPCone& sdpcone);
  virtual int createConstraints();
  virtual void roundNegatives();

  float _addLazyFraction;
  SumConstraints _sumConstraints;
};

#endif
