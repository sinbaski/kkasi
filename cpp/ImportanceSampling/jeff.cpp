#include <time.h>
#include <cmath>
#include <algorithm>
// #include <armadillo>
#include "SmallNumber.hpp"

using namespace std;
using namespace arma;

chi_squared_distribution<double> dist;
random_device gen;

extern void adjust_numbers(SmallNumber &V);

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

    /* adjust_numbers is good! */
    // vector<double> V = {
    // 	-14.509109561124117, -13.176534978491997, 0, 15.243856232806506, 15.720977487526168, 30.487712465613011	
    // };
    // SmallNumber X(V);
    // adjust_numbers(X);


    //This test case is passed.
    // XMatrix A(2, 2), B(2, 2), M;
    // A(0, 0) = 1.7533e+15;
    // A(0, 1) = 1.546e-29;
    // A(1, 0) = 2.0030e+14;
    // A(1, 1) = 0;

    // B(0, 0) = 1;
    // B(0, 1) = 3;
    // B(1, 0) = 2;
    // B(1, 1) = 4;

    // M = A * B;
    // M *= A;
    // M += B + A * A;

    // double data1[][4] = {
    // 	{19.403665, 4.077789e-01, 797.39273, 0.03255509},
    // 	{6.128070, 3.924188e-05,  27.35455, 1.36518070},
    // 	{7.012403, 6.138501e+00,  41.20729, 0.37090008}
    // };
    // double data2[][3] = {
    // 	{0.014663314,    0.187702,   1.3084776},
    // 	{0.002948789,    2.226695,   0.2997368},
    // 	{0.005873739, 1005.096251,   0.3484540},
    // 	{0.474892351,    1.205615, 124.4773247}
    // };
    // double data3[][3] = {
    // 	{0.0000418549, 0.7024121, 1.83917835},
    // 	{0.3144126534, 0.9380704, 0.07165514},
    // 	{4.6536517854, 0.5059325, 3.36968125},
    // };
    // double *p[] = {data1[0], data1[1], data1[2]};
    // double *q[] = {data2[0], data2[1], data2[2], data2[3]};
    // double *r[] = {data3[0], data3[1], data3[2]};
    // XMatrix A(p, 3, 4), B(q, 4, 3), C(r, 3, 3), M(3, 3);
    // M = A * B + C;
    // M += C * A * B;
    // M *= (C^2) + ((A * B)^3);
    

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

