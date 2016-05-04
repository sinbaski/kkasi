#include <time.h>
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <vector>
#include "GarchPQ.hpp"

using namespace std;
using namespace arma;

/**
   alpha: a vector of size 2+ that contain the coefficients of past squared returns.
   beta: a vector of size 1+ that contain the coefficients of past variances.
   This notation is consistent with Mikosch and Davis'es paper. on GARCH(p,q)
 */
mat& gen_rand_matrices(const vector<double> &alpha,
		       const vector<double> &beta,
		       int measure_index,
		       mat& A)
{
    static chi_squared_distribution<double> dist;
    static random_device gen;
    
    double z = dist(gen);
    A.row(0) = rowvec({alpha[0] * z + beta[0], beta[1], alpha[1], alpha[2]});
    A.row(1) = rowvec({1, 0, 0, 0});
    A.row(2) = rowvec({z, 0, 0, 0});
    A.row(3) = rowvec({0, 0, 1, 0});
    return A;
}

int main(int argc, char*argv[])
{
    vector<double> alpha({0.6, 0.1, 0.1});
    vector<double> beta({0.1, 0.05});
    
}
