#include <time.h>
#include <cmath>
#include <algorithm>
// #include <armadillo>
#include "SmallNumber.hpp"

using namespace std;
using namespace arma;

chi_squared_distribution<double> dist;
random_device gen;

XMatrix & gen_rand_matrix(const vector<double> &alpha,
			  const vector<double> &beta,
			  int measure_index, XMatrix &M)
{
    unsigned n = 2;
    double z = dist(gen);
    M.entry.resize(n);
    for_each(M.entry.begin(), M.entry.end(),
	     [n](vector<SmallNumber> &x) {
		 x.resize(n);
	     });
    M(0, 0) = alpha[1] * z + beta[0];
    M(0, 1) = alpha[2];
    M(1, 0) = z;
    M(1, 1) = 0;
    return M;
}

int main(int argc, char*argv[])
{
    vector<double> alpha({1.0e-7, 0.6, 0.1});
    vector<double> beta({0.1, 0.05});
//    double power = stod(argv[2]);
    double norms[2];
    double Lambda[2];
    vector<XMatrix> A(400);
    XMatrix P;

    for_each(A.begin(), A.end(),
    	     [&](XMatrix &M) {
    		 gen_rand_matrix(alpha, beta, 0, M);
    	     });
    P = A[0];
    for (unsigned i = 1; i < A.size(); i++) {
	P *= A[i];
    }
    return 0;
}

/**
   alpha: a vector of size 2+ that contain the coefficients of past squared returns.
   beta: a vector of size 1+ that contain the coefficients of past variances.
   This notation is consistent with Mikosch and Davis'es paper. on GARCH(p,q)
*/
// mat& gen_rand_matrix(const vector<double> &alpha,
// 		     const vector<double> &beta,
// 		     int measure_index,
// 		     mat& A)
// {
//     static chi_squared_distribution<double> dist;
//     static random_device gen;
    
//     double z(dist(gen));
//     A(0, 0) = z * alpha[0] + beta[0];
//     A(0, 1) = alpha[1];
//     A(1, 0) = z;
//     A(1, 1) = 0;
//     return A;
// }

// int main(int argc, char*argv[])
// {
//     vector<double> alpha({1.0e-7, 0.6, 0.1});
//     vector<double> beta({0.1, 0.05});
//     cube A(2, 2, stol(argv[1]));
//     mat X(2, 2);
//     double power = stod(argv[2]);
//     double norms[2];
//     double Lambda[2];
    
//     gen_rand_matrix(alpha, beta, 0, X);
    
//     for (unsigned k = 0; k < A.n_slices; k++) {
// 	gen_rand_matrix(alpha, beta, 0, A.slice(k));
// 	X *= A.slice(k);
// 	for (unsigned i = 0; i < X.n_rows; i++) {
// 	    for (unsigned j = 0; j < X.n_cols; j++)
// 		if (X(i, j) == 0) {
// 		    printf("A(%d, %d, %d) = 0\n", k, i, j);
// 		}
// 	}
//     }
//     norms[0] = pow(norm(X, "inf"), power);
//     norms[1] = pow(norm(X, 2), power);
//     unsigned long n = A.n_slices + 1;

//     transform(norms, norms+2, Lambda,
// 	      [n] (double x) {
// 		  return log(x)/double(n);
// 	      });
    
//     printf("Norms,%e,%e\n",
// 	   norms[0], norms[1]);
//     printf("Lambda,%e,%e\n", Lambda[0], Lambda[1]);

//     return 0;
// }

