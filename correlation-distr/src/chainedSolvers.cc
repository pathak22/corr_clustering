#include "WeightMatrix.h"
#include "Greedy.h"

#ifdef CPLEX
#include "LPSolver.h"
#endif

#include "BestOneEltMove.h"

#ifdef DSDP_SDP
#include "SDPSolver.h"
#endif

#ifdef CBUNDLE_SDP
#include "SDPSolverTwo.h"
#endif

void usage()
{
	cerr<<"chainedSolvers [log] [ordered list of solvers] "<<
		"[classifier decisions]\n";
}

int main(int argc, char* argv[])
{
	srand(time(NULL));

	if(argc < 3)
	{
		usage();
		return 1;
	}

	bool log = false;

	strings solvers;
	for(int ii = 1; ii < argc - 1; ++ii)
	{
		if(string(argv[ii]) == "log")
		{
			cerr<<"LOG MODE ON\n";
			log = true;
			continue;
		}
		solvers.push_back(argv[ii]);
	}

	cerr<<"Reading "<<argv[argc - 1]<<"\n";
	ifstream decisions(argv[argc - 1]);
	if(!decisions.good()) {
	  cerr <<"Error: cannot open weight file: " << argv[argc - 1] << std::endl;
	  exit(1);
	}
	WeightMatrix input(decisions, log?W_LOG:W_UNALTERED);

	WeightMatrix process = input;

	foreach(strings, solver, solvers)
	{
		cerr<<*solver<<"\n";

		if(*solver == "first")
		{
			Soon* ga = new Soon(process);
			process = ga->solve();
		}
		else if(*solver == "best")
		{
			BestLink* ga = new BestLink(process);
			process = ga->solve();
		}
		else if(*solver == "vote")
		{
			VotedLink* ga = new VotedLink(process);
			process = ga->solve();
		}
		else if(*solver == "pivot")
		{
			Pivot* ga = new Pivot(process);
			process = ga->solve();
		}
		else if(*solver == "first-id")
		{
			Soon* ga = new Soon(process);
			ga->setIdentityPerm();
			process = ga->solve();
		}
		else if(*solver == "best-id")
		{
			BestLink* ga = new BestLink(process);
			ga->setIdentityPerm();
			process = ga->solve();
		}
		else if(*solver == "vote-id")
		{
			VotedLink* ga = new VotedLink(process);
			ga->setIdentityPerm();
			process = ga->solve();
		}
		else if(*solver == "pivot-id")
		{
			Pivot* ga = new Pivot(process);
			ga->setIdentityPerm();
			process = ga->solve();
		}
#ifdef CPLEX
		else if(*solver == "lp")
		{
			LPSolver* la = new LPSolver(process);
			process = la->solve();
		}
		else if(*solver == "lazylp")
		{
		  LazyLPSolver* la = new LazyLPSolver(process, 0.5);
			process = la->solve();
		}
		else if(*solver == "ilp")
		{
			ILPSolver* la = new ILPSolver(process);
			process = la->solve();
		}
		else if(*solver == "lazyilp")
		{
			LazyILPSolver* la = new LazyILPSolver(process);
			process = la->solve();
		}
#endif
#ifdef DSDP_SDP
		else if(*solver == "sdp")
		{
			SDPSolver* sa = new SDPSolver(process);
			process = sa->solve();
		}
		else if(*solver == "sdp_nnc") //no non-negativity constraints
		{
		  SDPSolver* sa = new SDPSolver(process, false);
		  process = sa->solve();
		}
#endif
#ifdef CBUNDLE_SDP
		else if(*solver == "sdp2")
		{
		  SDPSolverTwo* sa = new SDPSolverTwo(process, false, false);
		  process = sa->solve();
		}
		else if(*solver == "sdp2l") //lazy
		{
		  SDPSolverTwo* sa = new SDPSolverTwo(process, false, true);
		  process = sa->solve();
		}

		else if(*solver == "sdp2t") //t for triangle inequalities
		{
		  SDPSolverTwo* sa = new SDPSolverTwo(process, true, false);
		  process = sa->solve();
		}
		else if(*solver == "sdp2lt") //triangle ineq. and lazy
		{
		  SDPSolverTwo* sa = new SDPSolverTwo(process, true, true);
		  process = sa->solve();
		}
#endif
#ifdef DSDP_SDP
		else if(*solver == "lazysdp")
		{
			SDPSolver* sa = new LazySDPSolver(process, 0.5, 7);
			process = sa->solve();
		}
		else if(*solver == "sumsdp")
		{
			SumLazySDPSolver* sa = new SumLazySDPSolver(process, 7);
			process = sa->solve();
		}
#endif
		else if(*solver == "boem")
		{
			ints clustering;
			int clusters = process.nodeLabels(clustering);
			if(clusters != -1)
			{
				BestOneEltMove* boem = new BestOneEltMove(input, clustering);
				process = boem->solve();
			}
			else
			{
				BestOneEltMove* boem = new BestOneEltMove(input);
				process = boem->solve();
			}
		} 	
		else if(*solver == "sa")
		{
			ints clustering;
			int clusters = process.nodeLabels(clustering);
			if(clusters != -1)
			{
				BestOneEltMove* boem = new BestOneEltMove(input, clustering);
				process = boem->solveSimAnneal();
			}
			else
			{
				BestOneEltMove* boem = new BestOneEltMove(input);
				process = boem->solveSimAnneal();
			}
		} 
		else if(*solver == "stats") {

		  //A trivial lower-bound: the sum over pairs of objects of the smaller of the plus and minus weights
		  double lb = 0;
		  double edges = 0;
		  
		  //print out the current objective value. For ILP this is exact, for LP, lazylp, SDP and sumsdp this is lower-bound, for sdp2 this is neither, for others upper-bound.
		  double accum = 0;
		  
		  double minErrors = 1e200; //min and max errors incident to one vertex
		  double maxErrors = 0; //todo: replace with better infinity
		  
		  double minSavings = 1e200; //min and max savings relative to singletons
		  double maxSavings = 0; //todo: replace with better infinity
		  for(int i = 0; i < input.nodes(); i++) {
		    double thisError = 0; //error incident to vertex i
		    double thisSavings = 0; //savings relative to all singletons
		    for(int j = 0; j < input.nodes(); j++) {
		      if (i == j)
			continue;
		      double error = (input.wplus(i,j) * (1 - process(i,j)) + input.wminus(i,j) * process(i,j));
		      if (i < j) {
			accum += error;
			lb += fmin(input.wplus(i,j), input.wminus(i,j));
			edges += input.wplus(i,j);
		      }
		      thisError += error;
		      thisSavings += input.wplus(i,j) - error;
		    }
		    if (thisError < minErrors)
		      minErrors = thisError;
		    if (thisError > maxErrors)
		      maxErrors = thisError;
		    
		    if (thisSavings < minSavings)
		      minSavings = thisSavings;
		    if (thisSavings > maxSavings)
		      maxSavings = thisSavings;
		    
		  }
		  cerr  << "Objective value: " << accum << " (depending on solver may be lower-bound, upper-bound or neither)" << std::endl;
		  cerr << "Min/max/average error incident to one vertex: " << minErrors << " " << maxErrors << " " << 2. * accum / input.nodes() << std::endl;
		  cerr << "Min/max/avg savings incident to one vertex relative to singletons: " << minSavings << " " << maxSavings << " " << 2. * (edges - accum) / input.nodes() << std::endl;
		  cerr << "Trivial lower bound: " << lb << std::endl;
		  cerr << "LP feasible solution: " << (lb + edges) / 2 << std::endl;
		  cerr << "Total positive weight: " << edges << std::endl;
		}
		else if(*solver == "print")
		{
		  std::ofstream file("clustering.mat");
		  process.write(file);
		}
		else if(*solver == "read")
		{
		  //load process matrix from file. "chainedSolvers sdp2 print" then "chainedSolvers read pivot boem" should be equivalent to "chainedSolvers sdp2 pivot boem", allowing one sdp solution to be rounded multiple times.
		  std::ifstream file("clustering.mat");
		  process.read(file);
		}
		else
		{
			cerr<<"Unrecognized solver: "<<*solver<<"\n";
			cerr<<"If this message is unexpected, check the Makefile to ensure the appropriate LP/SDP solvers are installed." << std::endl;
			return 1;
		}
		//print out the current objective value. For ILP this is exact, for LP, lazylp, SDP and sumsdp this is lower-bound, for sdp2 this is neither, for others upper-bound.

		double accum = 0;

		for(int i = 0; i < input.nodes(); i++) 
		{
			for(int j = 0; j < i; j++) 
			{
				double error = (input.wplus(i,j) * (1 - process(i,j)) + 
								input.wminus(i,j) * process(i,j));
				accum += error;
			}
		}

		cerr  << "Objective value after " << *solver<< ": " << accum << " (depending on solver may be lower-bound, upper-bound or neither)" << std::endl;
	}

	ints matrixLabels;
	int k = process.nodeLabels(matrixLabels);
	if(k == -1)
	{
		cerr<<"error: failed to cluster...\n";
		return 1;
	}

	cerr<<"k-clustered as "<<k<<" clusters\n";
	for(int ii = 0; ii < matrixLabels.size(); ii++)
	{
		cout<<matrixLabels[ii]<<"\n";
	}
}
