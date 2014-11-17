#include "SDPSolverTwo.h"
#include <queue>
#include "cmsymspa.hxx"
#include "clock.hxx"

using namespace CH_Matrix_Classes;
using namespace ConicBundle;

//structure for noting that a triple of vertices (i,j,k) violates either a triangle inequality or non-negativity constraint. _k=-1 denotes non-negativity constraint X_{ij} >= 0 and _k >=0 denotes a triangle constraint.
struct ConstraintVio {
  int _i;
  int _j;
  int _k;
  double _v; 
  ConstraintVio(int i, int j, double v) : _i(i), _j(j), _k(-1), _v(v) { assert(i<j && i >= 0); } //for bound constraints
  ConstraintVio(int i, int j, int k, double v) : _i(i), _j(j), _k(k), _v(v) {assert(i<j && i != k && j != k);  } //for triangle constraints
  bool operator<(const ConstraintVio & other) const {
    return _v < other._v;
  }
  bool IsTriangle() {
    return _k >= 0;
  }
};

SDPSolverTwo::SDPSolverTwo(WeightMatrix& wt, bool doTriangles, bool lazy):
  _input(wt), 
  _result(wt.nodes(), W_UNALTERED),
  VERBOSE(1),
  _nodes(_input.nodes()),
  _sdp(),
  _cbsolver(),
  _boundVars(_nodes, 0.),
  _triangleCons(),
  _nBoundVars(0),
  _numNewCon(0),
  _newOpA(),
  _newRHSides(0,0,0.),
  _newLBs(0,0,0.),
  _newUBs(0,0,0.),
  _doTriangles(doTriangles),
  _lazy(lazy)
{
}

SDPSolverTwo::~SDPSolverTwo()
{
}

void SDPSolverTwo::ReallyAddConstraints() {
  assert(_numNewCon == _newLBs.dim());
  assert(_numNewCon == _newUBs.dim());
  assert(_numNewCon == _newRHSides.dim());
  assert(_numNewCon == _newOpA.size());

  if (_numNewCon > 0) {
    _sdp.append_constraints(_newOpA, _newRHSides);
    _cbsolver.append_variables(_numNewCon, &_newLBs, &_newUBs, 0); //null costs since those are done automatically by MatrixSDPfunction

    _numNewCon = 0;
    _newOpA.clear();
    _newRHSides.init(0,0,0.);
    _newLBs.init(0,0,0.);
    _newUBs.init(0,0,0.);
  }
}

void SDPSolverTwo::AddConstraint(int i, int j) {
  assert(_boundVars(i, j) == 0);

  _boundVars(i,j) = 1;
  _nBoundVars++;
  
  double scale = 1; //The idea of this was to choose scale > 1 to emphasize the diagonal constraints, but it hurts convergence too much to be worth it.
  
  _newOpA[_numNewCon][0] = new CMsingleton(_nodes, i, j, scale);
  if (i == j) {
    _newRHSides.concat_below(scale);
    _newLBs.concat_below(CB_minus_infinity);
    _newUBs.concat_below(CB_plus_infinity);
  } else {
    _newRHSides.concat_below(0.);
    _newLBs.concat_below(CB_minus_infinity);
    _newUBs.concat_below(0);
  }
  _numNewCon++;
}

void SDPSolverTwo::AddTriangleConstraint(int i, int j, int k) {
  Triple t(i,j,k);
  assert(_triangleCons.find(t) == _triangleCons.end());
  assert(i != j && i != k && j != k);
  
  _triangleCons.insert(t);  
  
  //Make matrix for constraint X_ik + X_jk - X_ij <= 1
  Indexmatrix is(3,1,-1);
  Indexmatrix js(3,1,-1);
  Matrix vals(3,1,0.);
  is[0] = i; js[0] = k; vals[0] = 0.5;
  is[1] = j; js[1] = k; vals[1] = 0.5;
  is[2] = i; js[2] = j; vals[2] = -0.5;
  Sparsesym A(_nodes, 3, is, js, vals);


  (_newOpA[_numNewCon])[0] = new CMsymsparse(A);

  _newRHSides.concat_below(1);
  _newLBs.concat_below(0);
  _newUBs.concat_below(CB_plus_infinity);

  _numNewCon++;
}

WeightMatrix& SDPSolverTwo::solve()
{
  CH_Tools::Clock clock; //defined in the Conic Bundle somewhere
  int error;

  int nodes = _nodes;

  Symmatrix C(nodes, 0.); //Objective function matrix

  Symmatrix & bndCons = _boundVars; //which X_ij >=0 or X_ii = 1 constraints added so far
  int & nCon = _nBoundVars; //number of such constraints added so far

  double totalPlus(0.);

  for (int i = 0; i < nodes; i++)
    for (int j = 0; j < i; j++) 
      {
	totalPlus += _input.wplus(i,j);
	C(i,j) = _input.wdiff(i,j) / 2.;
      }
  //Set diagonal of cost matrix so that the constant term comes out right
  Matrix constCorrection(nodes, 1, totalPlus / Real(nodes));
  C-=sparseDiag(constCorrection);

  MatrixSDPfunction & mc = _sdp;
  mc.set_trace(Real(nodes), SDPtrace_fixed); //Add redundant trace X = n constraint, which allows the max cut SDP to be expressed as an eigenvalue maximization problem, enabling the use of spectral bundle method.
  SparseSDPCoeffVector cp;
  cp[0]=new CMsymdense(C); //first block has objective matrix C
  std::map<Integer,SparseSDPCoeffVector> opA; //opA[constraintNumber][blockNumber], or is it the other way around?
  Indexmatrix Xdim(1,1,nodes); //one n by n block
  mc.append_variables(Xdim,cp,opA); //add the primal variables, a n by n matrix, and corresponding objective matrix C. We pass a blank opA because we have no constraints yet
  cp.clear();
  opA.clear();
	
  mc.set_generating_primal(new DenseSDPPrimal(C));
  mc.set_out(&cout,0);

  MatrixCBSolver & cbsolver = _cbsolver;

  error=cbsolver.init_problem(0); //create dual problem with no constraints so far
  assert(!error);
  error=cbsolver.add_function(mc); //dual objective consists of the mc function created earlier
  assert(!error);


  for (Integer i=0;i<nodes;i++){ //add the constraints X_ii = 1 and X_ij >= 0
    for (Integer j=0;j<=i;j++){ 
      if (i == j || !_lazy || double(rand()) / RAND_MAX < 0.01) {
	AddConstraint(i,j);
      }
    }
  }

       
  cbsolver.set_min_weight(0.5); //was 1. IIRC this is a step-size parameter.
  cbsolver.set_max_weight(2); //was 1. 
  cbsolver.set_out(&cerr,1);
  cbsolver.set_inner_update_limit(3); //was 3
	
  SDPBundleParameters bp;
  bp.n_bundle_size=10; //was 15 . Meaning (from BaseSDPOracle.hxx) : number of old ritz vectors
  bp.n_new_subgradients=5; //was 5 Meaning: number of new ritz vectors
  bp.n_keep=7; //was 10 Meaning: min number of old ritz vectors to keep
  bp.n_aggregates=1; //was 1
  bp.aggregation_tolerance=0.001; 
  cbsolver.set_bundle_parameters(mc,bp);
	
  int cnt=0; //number of SDP solving iterations

  DenseSDPPrimal X; //current primal approx
  bool done; //for terminating do-while loop
  do { //Solve the SDP
	  
    //add the constraints we've queued up
    ReallyAddConstraints();

    //Call
    error=cbsolver.do_descent_step();
    assert(!error);
	  
    cnt++;

    error=cbsolver.get_approximate_primal(mc,X);
    assert(!error);
	  
    //Check to see what constraints are violated
    std::priority_queue<ConstraintVio> vios;
    std::priority_queue<ConstraintVio> triangleVios;

	  
    double maxDiagV(0.), maxCOffV(0.), maxUOffV(0.), maxUTV(0.), maxCTV(0), //diagonal, checked off-diag, unchecked off-diag, triangle inequality violations
      totDiagV(0.), totCOffV(0.), totUOffV(0.), totUTV(0.), totCTV(0); 
    int nBigTV(0); //number of big triangle violations
    int nUselessCOff(0); //number of checked constraints that are far from tight
    int nUselessCT(0); //ditto for triangle
    for(int ii = 0; ii < nodes; ++ii) {
      for(int jj = ii; jj < nodes; ++jj) {
	if(ii == jj) {
	  double viol = fabs(X(ii,jj) - 1.);
	  if (viol > maxDiagV) {
	    maxDiagV = viol;
	  }
	  totDiagV += viol;
	} else {
	  double viol = X(ii,jj) < 0. ? -X(ii,jj) : 0.;
	  if(bndCons(ii,jj) != 0) {
	    //checked
	    if (viol > maxCOffV) {
	      maxCOffV = viol;
	    }
	    totCOffV += viol;
	    if (X(ii,jj) > 0.01)
	      nUselessCOff++;
	  } else {
	    //no constraint for this one
	    if (viol > maxUOffV) {
	      maxUOffV = viol;
	    }
	    if (viol > 0.0001) {
	      vios.push(ConstraintVio(ii,jj,viol));
	    }
	    totUOffV += viol;
	  }
	}
      }
    }

    bool trianglesChecked = false;
    int numLargeXij = 0; //number of Xij >= 0.5, which determines runtime of triangle checking
    //Compute the violated triangle inequalities, but only after everything else has stabilized since this step, which is
    // cubic in the number of vertices, is expensive.
    if (_doTriangles && (totDiagV < 0.1 * nodes || !_triangleCons.empty())) {
      trianglesChecked = true;

      //examine the "checked" triangles
      for(std::set<Triple>::iterator t = _triangleCons.begin(); t!=_triangleCons.end(); t++) {
	int i = t->_i;
	int j = t->_j;
	int k = t->_k;
	double preViol = X(i, k) + X(k, j) - X(i, j) - 1. ;
	double violT = fmax(preViol, 0.);
	assert(_triangleCons.find(Triple(i,j,k)) != _triangleCons.end());
	//unchecked
	totCTV += violT;
	maxCTV = fmax(maxCTV, violT);
	if (-preViol >= 0.05) {
	  nUselessCT++;
	}
      } //end loop over "checked" triangles

      for(int ii = 0; ii < nodes; ++ii) {
	for(int kk = 0; kk < nodes; ++kk) {
	  if (ii == kk || X(ii,kk) <= 0.5) //If vio at least one is > 0.5
	    continue;
	  numLargeXij++;
	  for(int jj = 0; jj < nodes; jj++) {
	    if ( jj == ii || jj == kk)
	      continue; //triangle ineq. trivially true if not distinct
	    int i = min(ii,jj); //ensure i < j so find works ok.
	    int j = max(ii,jj);
	    assert(i != j);
	    double preViol = X(i, kk) + X(kk, j) - X(i, j) - 1. ;
	    double violT = fmax(preViol, 0.);
		
	    if (violT > 0 && _triangleCons.find(Triple(i,j,kk)) == _triangleCons.end()) {
	      //not checked
	      totUTV += violT;		    
	      maxUTV = fmax(maxUTV, violT);
		  
	      if (violT > fmax(maxCTV, maxCOffV)) {
		nBigTV++;
		triangleVios.push(ConstraintVio(i,j,kk,violT));
	      }
	    }
	  }//for jj
	} //for kk
      } //for ii
    } //if we should check triangles

	  
    double aveDiagV = totDiagV / nodes;
    //Now add new constraints!
    int added(0), tAdded(0);
    if (cnt > 5 && maxDiagV < 2) { //give it a few iterations to get in the right ballpark before adding constraints
      double cutoff = (aveDiagV + maxCOffV)*(aveDiagV + maxCOffV) / 2;
      double addProb = 1; // maxDiagV+maxCOffV < 1 ? 1 : 0.01;
	    

      if(vios.size() > 0 && vios.top()._v > cutoff) {
	int maxToAdd = nCon;
	while(vios.size() > 0 && added < maxToAdd) {
	  const ConstraintVio & vio = vios.top(); 
	  if(double(rand()) / RAND_MAX < addProb) {
	    AddConstraint(vio._i,vio._j);
	    added++;
	  }
	  vios.pop();
	}
      }

      double tCutoff = fmax(fmax(fmax(3*aveDiagV, maxCOffV*2), maxUOffV*2), maxCTV*1.3);
      if (_doTriangles && trianglesChecked && tCutoff < 0.2) {
	      
	if (!triangleVios.empty() && triangleVios.top()._v > tCutoff) {
	  //It's reasonably well converged; consider adding triangle constraints
	  int maxToAdd = min((5*nCon - _triangleCons.size())/2, nodes * nodes / 40) ;
	  double probToAdd = fmin(1, double(maxToAdd) / triangleVios.size()); 
	  while(!triangleVios.empty()) { // &&  triangleVios.top()._v > maxCTV
	    if (double(rand()) / RAND_MAX < probToAdd) { //was 1
	      ConstraintVio cv = triangleVios.top();
	      if(_triangleCons.find(Triple(cv._i,cv._j,cv._k)) == _triangleCons.end()) //need to check this explicitly since with random sampling we might enqueue the same ijk twice.
		{ 
		  AddTriangleConstraint(cv._i, cv._j,cv._k);
		  tAdded++;
		}
	    }
	    triangleVios.pop();	      
	  } //while adding triangles
	} //if add some triangles
      } //if OK to add triangles

    } //if enough iterations gone by to add constraints

	  
    cerr.precision(4);
    cerr 
      << " diagV: " << maxDiagV << "/" << aveDiagV  
      << " CoffV: "  << nCon - nodes << "/" << added << "/" << nUselessCOff << "/" << maxCOffV << "/" << totCOffV  
      << " UoffV: " << maxUOffV << "/" << totUOffV;
    if(trianglesChecked) {
      //	  ConstraintVio tv1(-7,-6,-5,0);
      //const ConstraintVio & tv = triangleVios.size() > 0 ? triangleVios.top() : tv1;
      cerr
	<< " CTri: " << _triangleCons.size() << "/" << tAdded <<  "/" << nUselessCT << "/" << maxCTV << "/" << totCTV 
	<< " UTri: "  << nBigTV << "/" << maxUTV << "/" << totUTV;  
    //  "/" << tv._i << "," << tv._j << "," << tv._k <<
    }
    cerr << std::endl;

	    
    /*
      if(cnt==20) { //print out dual
      WeightMatrix dual(nodes, W_UNALTERED);
      Matrix y;
      cbsolver.get_center(y);
	    
      constraintNum = 0;
      for (Integer i=0;i<nodes;i++){ //set bounds for dual vars for constraints X_ii = 1 and X_ij >=0
      for (Integer j=0;j<=i;j++){ 
      dual(i,j) = y(constraintNum);
      constraintNum++;
      }
      }
      std::ofstream file("dual.mat");
      dual.write(file);
      }
    */

    bool earlyFinish = false; //maxDiagV < 0.005 && maxUOffV < 0.005 && maxCOffV < 0.005;
    const double maxHours = 120;
    done = earlyFinish || cbsolver.termination_code() || clock.time().roundsecs() > maxHours*3600;
  } while (!done);
	
  cerr.precision(8);
  cerr << "SDP lower-bound: " << (-cbsolver.get_objval()) << std::endl;

  cbsolver.print_termination_code(cerr);

  for(int ii = 0; ii < nodes; ++ii)
    {
      for(int jj = 0; jj <= ii; ++jj)
	{
	  if(ii != jj)
	    {
	      _result(ii, jj) = X(ii,jj);
	    }
	}
    }
	
  return _result;
}

