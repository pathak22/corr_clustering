=========================================
Correlation Clustering
=========================================

.. To convert this document to HTML, use command line tool rst2html.

.. contents::

Setup
=======

This is the software for [ElsnerSchudy09]_. It contains software for
creating, analyzing and solving correlation clustering problems.

License
----------

Our code is available under the terms of the GPL version 2.0 or later. In addition, as a special exception, we also give permission to link our code with linear programming and semi-definite programming solvers regardless of the licenses of such solvers. When exercising your distribution rights under the GPL you need not include the source code for such solvers. Note that people who make modified versions of our code are not obligated to grant this special exception for their modified versions; it is their choice whether to do so. The GPL gives permission to release a modified version without this exception; this exception also makes it possible to release a modified version which carries forward this exception.

Before distributing binaries please ensure that you do not violate the licenses of any LP or SDP solvers that you link against; two examples follow. The CPLEX license presumably does not allow distribution of binaries containing CPLEX. The ConicBundle is licensed under the GPL without the above exception, therefore only distribute binaries containing ConicBundle if all solvers linked with are GPL compatible.

We are not lawyers and this is not legal advice.

Downloading and compiling
--------------------------

First download the clustering.tgz tarball from *cs.brown.edu/~melsner*
and "tar xzf" it, yielding the *correlation* directory.

The executable that actually does the clustering, *chainedSolvers*, is
written in C++. See the C++ building section of this manual for
building instructions, or pray and then run the following:

::

	cd correlation
	mkdir bin32
	mkdir lib32
	make chainedSolvers

Our support code that does everything else, such as evaluating a
clustering, generating artificial data, and visualization, is written
in Python. Please set your python path to include all necessary
packages (notably the waterworks utility library and pylab). Also edit
the script *script/megam.py* to point to the correct *megam* binary.

External Software
-----------------

Some features of this software require external packages. You can do
without most of these packages if you are willing to forgo various
features of the software:

   [waterworks]_ : David McClosky and others (python utility package,
   including ClusterMetrics library for evaluating clusterings)

   pylab plotting library, which in turn needs [matplotlib]_ (plotting)   

   [megam]_ : Hal Daume III (max-ent classifier)

   python package for disentangling [chat]_ : Micha Elsner

   [CPLEX]_ : ILOG Inc. (solves LP and ILP problems)

   [DSDP]_ : Benson, Ye and Zhang (solves SDPs)

   [ConicBundle]_ : Christoph Helmberg (solves SDPs) We used version
   0.2i of Conic Bundle.

Notes
=========

Program descriptions below only show the most useful features. To
get a usage message, try running the program without arguments, or with the
flag *-h*.

In general (but not always, thanks CPLEX!), programs write the
output of their computation to stdout and everything else (debugging
notes, iteration tickers, bonus information) to stderr. Most of the
time you will want to use the output redirect *>*.

Making data
=============

Making Artificial Points
------------------------

A dataset is essentially a vector of points and their cluster
labels, in the format::

		   [label] [feat val]*

*label* determines the cluster to which the point belongs, and is a
positive integer. *feat* can be any string, and *val* is a real
number.

Each datafile usually contains two sections, one for training and one
for testing, separated by a newline.

There are two primary ways to make artificial datasets: featureless data
(just cluster labels) and Gaussian data (points drawn from spherical
clusters). Both methods pick the number of clusters and cluster sizes
stochastically (some parameters of the distributions are controlled by
command line arguments).

To get featureless points, use:

::

	python script/genFeatureless.py -n [number of points]

The generative model for featureless points is:

::

	n (number of points) ~ command line
	k (number of clusters) ~ command line | log_1.5(n)
	alpha (prior parameter) ~ command line | 1
	p (k-dimensional vector of cluster probabilities) ~ Dirichlet(alpha) | uniform if command line parameter "balanced"
	z (number of points in each cluster) ~ Multinomial(p, n)

To get Gaussian points, use:

::

	python script/genGaussians.py -n [number of points] -t [number of training points] -f [number of features]

The generative model for Gaussian points is:

::

	z (number of points in each cluster) ~ same as featureless
	f (number of features) ~ command line | 1
	vv (optional parameter controlling feature variances) ~ command line
	sigma (k*f matrix of variances of each cluster) ~ command line | InverseGamma(vv, 1)
	mu (k*f matrix of means of each cluster) ~ Gaussian(0, sigma)
	F[i,j] (feature j of point i in cluster k_i) ~ Gaussian(mu[k_i,j], sigma[k_i,j])


To make a clustering with all clusters of size 1, use

::

	python script/genSingletons.py -n [number of points]


Making Artificial Weight Matrices
---------------------------------

To make a weight matrix, you first need to generate the points as described in the previous section. Then you need to run a classifier.

We provide three classifiers for general experimentation. To find the
names of the classifiers to use as *-c* arguments, run the script with
*-a*; it will print a list.

Our first classifier makes random errors. This is sometimes a useful
algorithmic model (and has been studied theoretically: see the
paper). This is also the only classifier which does anything useful
with featureless synthetic points. You can set the parameters
*epsilon*, *A* and *B* via command line options (the flag
*--simple* creates a 0-1 instance where *A* and *B* are ignored).

+------------------+---------------------+------------+
| True state       |  Classifier decision             |
+==================+=====================+============+
| Same cluster     | with p(1 - epsilon) | Beta(A,B)  |
|                  +---------------------+------------+
|                  | with p(epsilon)     | Beta(B,A)  |
+------------------+---------------------+------------+
| Diff cluster     | with p(1 - epsilon) | Beta(B,A)  |
|                  +---------------------+------------+
|                  | with p(epsilon)     | Beta(A,B)  |
+------------------+---------------------+------------+

The other two classifiers actually learn models from the data, and
require points with features. The first (Naive Bayes) assumes all features are
independent samples from Gaussians, and learns Gaussian distributions
on the differences of features for within-class and cross-class
instances. The second (Max Ent) learns a logistic regression on
feature differences instead.

::

	python script/classify.py -c [classifier] [data file]

The stderr stream output will look like this:

::

	test classifier performance
	P: 85.71 R: 100.00 F: 92.31 Acc: 90 (9/10)

Here P is precision of *same cluster* class (number of correct *same
cluster* decisions / number of *same cluster* decisions), R is recall
(number of correct *same cluster* decisions / number of true *same
cluster* edges), F is F-score (geometric mean of P and R), Acc is
accuracy (number correct / number of edges). In particular, watch out
for values like F = 0, Acc = .9: this means that the classifier is
useless-- it will always say *different cluster*, but since that's
usually the correct decision, accuracy is misleadingly high.

Making the Twenty Newsgroups Weight Matrices
--------------------------------------------

The best way to get the newsgroup weight matrices we used is directly
off the web, at *cs.brown.edu/~melsner*. We are providing our
newsgroup processing code in order to make our work replicable, not
because it is particularly general, elegant or effective.

If you are insistent on actually running the newsgroup code, first you
have to actually get the [mini_newsgroups]_ dataset from the UCI
machine learning repository.

Next, edit *script/newsgroup.py*, setting the path to your newsgroup
directory and a filename for the term frequency dump file which the
script will create. Run the script::

	   python script/newsgroup.py

Now edit *script/newsgroupToDataset.py* and set the same path to the
dump file. Also set a path to an empty directory where the program
will dump the data files. This program will transform the newsgroup
data into data files with integer cluster labels and term/count
features, and create five training/testing splits of the data.

::

	   python script/newsgroup.py

Now comes the really ugly part; we did the Latent Semantic Analysis decomposition by hand, in Matlab. There are plenty of ways to reimplement this using any linear algebra package, though.

Here's what we did: use *script/writeSparseMat.py* to write each data
matrix into a Matlab sparse matrix file. Now use Matlab to import the
file. Run the following Matlab commands::

	  [u, s, v] = svds(mat, 200)
	  save '<lsa-filename>' -ascii u

Now add these features back to the original dataset using::

	python script/lsa.py [original dataset] [lsa-filename] > [augmented file]

Finally, you can run the newsgroup classifier::

	python script/classify.py -c mxnews [lsa-augmented file]

Making the Chat Weight Matrices
--------------------------------------------

Make sure you have a copy of the chat disentanglement package
([chat]_). Follow the instructions to create a set of predictions for
your dataset (eg, using *classifierTest*). Now use::

	 python script/correlationClusteringData.py [chat file] [predictions file] [keys file] [output true labels] [output weight matrix]

Note that the "true labels" file will contain a fake training section
with a single point.

Viewing Weight Matrices
--------------------------------------------
To view a weight matrix such as one output by classify.py, use

::

	python script/colorMat.py [matrix file]

The matrix is color-coded. Red indicates 1, i.e. the classifier
believes the points belong in the same cluster with 100 percent
probability. Blue indicates 0, i.e. the points definitely belong in
different clusters.

This utility can also be used to view the output of the "print"
command in chainedSolvers.

You can also view the weight matrix and a clustering solution all at
once, by running::

	python script/colorMatAndTruth.py [matrix file] [solution file]

The solution file can be the original datafile from which the weight
matrix was produced, or the clustering provided by one of our
solvers. The display puts the weights in the upper triangle of the
matrix, and the solution in the lower triangle.

Finding Clusterings
=====================

Building
----------

Our Makefile has been tested on our x86 Debian GNU/Linux systems
only. Use on other platforms may require changes to the Makefile. The
Makefile automatically identifies whether the machine compiled on is
32 or 64 bit.

Our code supports two SDP solvers, DSDP and Conic Bundle, and one LP
solver, CPLEX. Without those solvers our code will still compile and
run, but LP and SDP based lower bounds will not be available. If you
have installed one or more of these solvers and wish to use them, set
the appropriate directories as documented in the Makefile.

Before building create a *bin32* directory and a *lib32* directory (or
64-bit equivalents). For convenience, it's nice to symlink *bin* to
*bin32*.

To make a binary such as chainedSolvers, go to the main directory and
"make chainedSolvers". The binary will be placed in the bin32 or bin64
directory as appropriate.

All these tools basically take a matrix file as input and write a
vector of cluster indices to stdout.

Using Our Recommended Solver
----------------------------

You can run the heuristic we recommend in the paper in the following way::

	bin/chainedSolvers log vote boem [matrix] > [solution]

In our experiments, we do this 100 times, check the objective values,
and take the best (lowest) objective.

Chained Solvers
---------------

You can run a long sequence of solvers (for instance, to solve an SDP
and then round the solution to integrality) using the *chainedSolvers*
program. Most of the solvers treat the output of the previous solver
as if it were the input matrix. The local search solvers BOEM and
simulated annealing act differently, treating the output of the
previous solver as an initial clustering to improve. The local search
solvers use one huge cluster as the initial clustering if run as the
first solver.

::

	bin/chainedSolvers [solver_1..solver_n] [matrix]

To preprocess the edge weights by taking logarithms, add "log" as the
first argument to chainedSolvers. (Actually you can hide "log" in the
middle of the solver list if you want to accomplish the same thing
with extra confusion.)

For instance, to run SDP, round with voting, and
improve the solution with best one-element, use:

::

	bin/chainedSolvers log sdp2 vote boem [matrix]

There are three additional special "solvers". The "stats" solver
prints information about the current solution and input. The "print"
solver prints the current solution to a file in the current directory
suitable for viewing with the colorMat.py script. It is often useful
to put "print" in the middle of a list of solver, e.g. to output the
SDP solution before rounding. The "read" solver restores the state as
of the last "print" solver. This is useful for reusing the painfully
slow SDP solutions. For example, first run:

::

	bin/chainedSolvers log sdp2 print [matrix]

to write the SDP solution to a file, currently hard-coded to "clustering.mat". Then round it, e.g.

::

	bin/chainedSolvers log read vote boem [matrix]

For a net result equivalent to:

::

	bin/chainedSolvers log sdp2 vote boem [matrix]

The advantage of using read and print is you can use the same SDP
solution multiple times.

To get a full list of solver names the program will accept, or tweak
construction arguments to any of the solvers, you'll need to edit the code.

Evaluation
============

The chainedSolvers application automatically prints objective function
values after each solution step for ease of debugging. To evaluate a
solution against the ground truth, run:

::

	python script/evaluate.py [data] [matrix] [output_1 .. output_n]

(The output files should be vectors of node indices. Make sure these
files don't contain log statements from CPLEX or something.)

The output will look something like this:

::

	True clustering has 2 clusters
	Objective value of truth: 2.44759404055
	Best Rand:
	File: data/featureless1/gpivot6
	Clusters: 2
	Objective: 2.45
	Objective (log): -18

	Some edge-counting metrics:
	Rand index (max 1): 1
	Jaccard index (max 1): 1
	Mirkin metric (min 0): 0
	(Same cluster) Prec: 1 Rec: 1 F: 1

	Some node-counting metrics:
	One-to-one match (max 1): 1
	Many-to-one match (max 1): 1
	Variation of information (min 0, max 2.32): 0
	Normalized mutual information (0-1): 1

For definitions of the metrics used, see the pydoc for the
ClusterMetrics package. Most of the metrics are defined in [Meila99]_.

References
----------

.. [ElsnerSchudy09] Micha Elsner and Warren Schudy. "Bounding and Comparing Methods for Correlation Clustering Beyond ILP". ILP-NLP '09.

.. [waterworks] http://www.cs.brown.edu/~dmcc/software/waterworks/

.. [mini_newsgroups] http://archive.ics.uci.edu/ml/databases/20newsgroups/20newsgroups.html

.. [Meila99] Marina Meila. "Comparing Clusterings". UW Statistics Technical Reports, COLT '03. http://www.stat.washington.edu/mmp/www.stat.washington.edu/mmp/Papers/compare-colt.pdf

.. [megam] Hal Daume III. Paper at http://pub.hal3.name#daume04cfg-bfgs.pdf, program at http://hal3.name/megam

.. [matplotlib] http://matplotlib.sourceforge.net/

.. [chat] Paper: Micha Elsner and Eugene Charniak. "You Talking To Me? A Corpus and Algorithm for Conversation Disentanglement". ACL '08. Software: http://cs.brown.edu/~melsner/chat-distr.tgz

.. [CPLEX] http://www.ilog.com/products/cplex/

.. [DSDP] Steven J. Benson, Yinyu Ye and Xiong Zhang. "DSDP5: Software For Semidefinite Programming". Tech report, 2005. Software: http://www.mcs.anl.gov/hs/software/DSDP/
 
.. [ConicBundle] Christoph Helmberg. "Semidefinite programming for combinatorial optimization". Tech report, 2000. http://www-user.tu-chemnitz.de/~helmberg/
