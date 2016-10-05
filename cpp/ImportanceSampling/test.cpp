#include <stdio.h>
#include <time.h>
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <iostream>
#include <random>
#include "ExtremeNumber2.hpp"
#include "XMatrix.hpp"

using namespace std;
using namespace arma;

normal_distribution<double> dist;
// student_t_distribution<double> dist(3);
// default_random_engine gen;
random_device gen;

extern void adjust_numbers(ExtremeNumber &V);

XMatrix & gen_rand_matrix(const vector<double> &alpha,
			  const vector<double> &beta,
			  int measure_index, XMatrix &M)
{
    unsigned p = beta.size(), q = alpha.size() - 1;
    unsigned n = p + q - 1;
    double z = pow(dist(gen), 2);
    ExtremeNumber e(0);
    M.entry.resize(n);
    for_each(M.entry.begin(), M.entry.end(),
	     [&](vector<ExtremeNumber> &row) {
		 row.resize(n, e);
	     });
    M(0, 0) = alpha[1] * z + beta[0];
    for (unsigned i = 1; i < p; i++) {
	M(0, i) = beta[i];
    }
    for (unsigned i = p; i < n; i++) {
	M(0, i) = alpha[i - p + 2];
    }
    for (unsigned i = 1; i < p; i++) {
	M(i, i - 1) = 1;
    }
    if (n > 1) M(p, 0) = z;
    for (unsigned i = p + 1; i < n; i++) {
	M(i, i - 1) = 1;
    }
    return M;
}

ExtremeNumber estimateLambda(const vector<double>& alpha,
		      const vector<double>& beta,
		      double xi,
		      unsigned long n,
		      unsigned long K,
		      ExtremeNumber &sd)
{
    vector<ExtremeNumber> results(K);
#pragma omp parallel for
    for (unsigned k = 0; k < K; k++) {
	vector<XMatrix> A(n);
	for_each(A.begin(), A.end(),
		 [&](XMatrix &M) {
		     gen_rand_matrix(alpha, beta, 0, M);
		 });
	XMatrix P = A[0];
	for (unsigned i = 1; i < A.size(); i++) {
	    P *= A[i];
	}
	long power;
	double m;
	mat X = P.comptify(&power);
	m = norm(X, "inf");

	results[k] = ExtremeNumber(m);
	results[k] ^= xi;
	results[k].mylog += power * xi;
    }

    sort(results.begin(), results.end());
    /* Importance sampling using the empirical distribution */
    results.emplace(results.begin(), 0.0);
    double n_inv = 1/(double)n;
    unsigned long M = K;
    vector<ExtremeNumber> Q(K + 1);
    vector<ExtremeNumber> Y(M);
    partial_sum(
	results.begin(), results.end(),
	Q.begin(),
	[=](const ExtremeNumber& s, const ExtremeNumber& x) {
	    return s + (x^0.5);
	    // return s + (x^n_inv);
	});
    results.erase(results.begin());
    Q.erase(Q.begin());
    ExtremeNumber mean(0);
    sd = 0;
#pragma omp parallel for shared (gen)
    for (size_t i = 0; i < M; i++) {
	uniform_real_distribution<double> unif(0, static_cast<double>(Q.back()));
	double U;
#pragma omp critical	    
	U = unif(gen);
	size_t k = upper_bound(Q.begin(), Q.end(), U) - Q.begin();
	Y[i] = results[k];
	// mean += Y[i]/fabs(log10(Y[i])) / (double)M;
	mean += (Y[i]^0.5)/(double)M;
    }

// 	vector<ExtremeNumber> J(K);
// #pragma omp parallel for shared (gen)
// 	for (size_t j = 0; j < K; j++) {
// 	    uniform_real_distribution<double> unif(0, static_cast<double>(Q.back()));
// 	    double U;
// #pragma omp critical	    
// 	    U = unif(gen);
// 	    size_t k = upper_bound(Q.begin(), Q.end(), U) - Q.begin();
// 	    J[j] = results[k];
// 	}
// 	Y[i] = accumulate(J.begin(), J.end(), ExtremeNumber(0.0),
// 			   [=](const ExtremeNumber& s, const ExtremeNumber& x)
// 			   {
// 			       return x/fabs(log10(x))/(double)K + s;
// 			   });
//     mean += Y[i] / (double)M;
// }
    // double Elogx = accumulate(results.begin(), results.end(), 0.0,
    // 			      [=](double s, const ExtremeNumber& x) {
    // 				  return s + fabs(log10(x))/(double)K;
    // 			      });
    // mean *= Elogx;
    // ExtremeNumber mean =
    // 	accumulate(.cbegin(), results.cend(),
    // 		   ExtremeNumber(0.0))/q;
    sd = accumulate(
    	Y.cbegin(), Y.cend(), ExtremeNumber(0),
    	[M, &mean](const ExtremeNumber &y, const ExtremeNumber &x)
    	{
    	    return y + (abs(x - mean)^2) / (double)M;
    	});
    
    sd ^= 0.5;
    return mean;
//    return mean() * pow(10, mean.power());
}

void test_cases(void)
{
   /* adjust_numbers is good! */
    // vector<double> V = {
    // 	-14.509109561124117, -13.176534978491997, 0, 15.243856232806506, 15.720977487526168, 30.487712465613011	
    // };
    // ExtremeNumber X(V);
    // adjust_numbers(X);


    //This test case is passed.
    // XMatrix A(2, 2), B(2, 2), M;
    // mat X;
    // int p;
    // A(0, 0) = 1.7533e+15;
    // A(0, 1) = 1.546e-29;
    // A(1, 0) = 2.0030e+14;
    // A(1, 1) = 0;

    // B(0, 0) = 1;
    // B(0, 1) = 3;
    // B(1, 0) = 2;
    // B(1, 1) = 4;

    // M = A * B;
    // X = M.comptify(&p);
    // X.print();
    // M *= A;
    // X = M.comptify(&p);
    // X.print();
    // M += B + A * A;
    // X = M.comptify(&p);
    // X.print();

    // Test case passed.
    // XMatrix A(2, 3), B(3, 2);
    // mat V;
    // int p;
    // A(0, 0) = 1.2;
    // A(0, 1) = 0.9;
    // A(0, 2) = 0.01;
    // A(1, 0) = 1.5;
    // A(1, 1) = 2.2;
    // A(1, 2) = 0;

    // B(0, 0) = 3.9;
    // B(0, 1) = 0.3;
    // B(1, 0) = 10.19;
    // B(1, 1) = 1.22;
    // B(2, 0) = 1.2e-3;
    // B(2, 1) = 6.5;
    
    // XMatrix C = A * B;
    // V = C.comptify(&p);
    // V.print();

    // XMatrix M = (C^10);
    // V = M.comptify(&p);
    // V.print();

    // test case Not run.
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
}

int main(int argc, char*argv[])
{
    vector<double> alpha({9.2e-6, 0.08835834, 0.09685783});
    // vector<double> alpha({1.0e-7, stod(argv[3])});
    vector<double> beta({0.6543018});
    ExtremeNumber expected;
    ExtremeNumber sd;
    
    cout << "alpha[0]= "  << alpha[0] << ", alpha[1]=" <<
    	alpha[1] << ", alpha[2]=" << alpha[2] <<
    	", beta[1]=" << beta[0] << endl;
    cout << "n = " << argv[3] << endl;
    cout << argv[4] << " iterations" << endl;

    // double nu = stod(argv[1]);
    // Lambda = estimateLambda(
    // 	alpha, beta, nu, stoul(argv[2]),
    // 	stoul(argv[3]), bounds);
    // printf("%e    %e\n", nu, Lambda);
    double n = stod(argv[3]);
    for (double nu = stod(argv[1]); nu < stod(argv[2]); nu += 0.1) {
    	expected = estimateLambda(
    	    alpha, beta, nu, stoul(argv[3]), stoul(argv[4]), sd);
	double rel_err = static_cast<double>(sd / expected);
    	printf("%e    %e    %.4f\n", nu, log(expected)/n, rel_err);
    }
    

    // printf("Lambda(%s) = %.4e\n", argv[1], Lambda);
    // printf("lambda(xi)^n = %.4fE%+ld (%.4fE+%ld, %.4fE+%ld)\n\n",
    // 	   bounds[1].first, bounds[1].second,
    // 	   bounds[0].first, bounds[0].second,
    // 	   bounds[2].first, bounds[2].second
    // 	);
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

