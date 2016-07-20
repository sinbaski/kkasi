#include <vector>
#include <algorithm>
#include <armadillo>

using namespace std;
using namespace arma;

#typedef vec (*proposal_density)(const mat &X);
#typedef vec (*target_density)(const mat &X);

/**
   @X: Each row of the matrix is a set of arguments at which the density is to be evaluated.
   @ret: The i-th element of @ret corresponds to the density evaluated at the i-th row of @X
 */
vec garch21_h_func(cosnt mat &X)
{
    return vec(0);
}



void metropolisHastings(proposal_density prop, target_density target, vector< vector<double> > &samples)
{
    unsigned long n = samples[0].size();
    unsigned long m = ceil(n*0.1);
    
}

