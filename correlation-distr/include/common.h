#ifndef COMMON_H
#define COMMON_H

//includes
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include <math.h>
#include <assert.h>

#include <tr1/unordered_map>
#include <tr1/unordered_set>

using std::cin;
using std::cout;
using std::cerr;

using std::string;

using std::istream;
using std::ifstream;
using std::ostream;

using std::istringstream;
using std::ostringstream;

//the naming conventions here should be obvious, right?
typedef std::vector<int> ints;

typedef std::tr1::unordered_set<int> intSet;
typedef std::tr1::unordered_map<int, intSet> intToIntSet;

typedef std::tr1::unordered_map<int, double> intToDouble;
typedef std::tr1::unordered_map<int, intToDouble> intToIntToDouble;

typedef std::vector<int*> intPtrs;
typedef std::vector<double*> doublePtrs;

typedef std::vector<string> strings;

//MJ's loop macros
// foreach is a simple loop construct
//
//   STORE should be an STL container
//   TYPE is the typename of STORE
//   VAR will be defined as a local variable of type TYPE::iterator
//
#define foreach(TYPE, VAR, STORE) \
   for (TYPE::iterator VAR = (STORE).begin(); VAR != (STORE).end(); ++VAR)

// cforeach is just like foreach, except that VAR is a const_iterator
//
//   STORE should be an STL container
//   TYPE is the typename of STORE
//   VAR will be defined as a local variable of type TYPE::const_iterator
//
#define cforeach(TYPE, VAR, STORE) \
   for (TYPE::const_iterator VAR = (STORE).begin(); VAR != (STORE).end(); ++VAR)

//get a random permutation of indices 0..n
bool isIdentity(ints& perm);
void randPermutation(ints& iperm, int size);

//turn an int into a string
string intToString(int x);

#define approxEq(p, q) (fabs((double)(p) - (double)(q)) < 1e-5)

#endif
