#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <math.h> 
#include <algorithm>
#include <armadillo>
#include <iostream>
#include <random>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_erf.h>
#include "ExtremeNumber2.hpp"
#include "XMatrix.hpp"

using namespace std;
using namespace arma;

random_device gen;

#define SHIFT_PARAM 0.2

// struct param_structure
// {
//     const vector<double>& a;
//     const vector<double>& b;
// };


// double shift_param_def_func(double x, void *params)
// {
//     struct param_structure * p = (struct param_structure *)params;
//     const vector<double>& a = p -> a;
//     const vector<double>& b = p -> b;

//     double w = (a[2] + b[0])/(1 - a[1]);
//     double t1 = sqrt(0.2e1);
//     double t3 = exp(x * a[2]);
//     double t6 = exp(x * b[0]);
//     double t10 = sqrt(-0.4e1 * x * a[1] + 0.2e1);
//     double t11 = sqrt(w);
//     double t14 = gsl_sf_erf(t11 * t10 / 0.2e1);
//     double t21 = sqrt(-0.4e1 * x + 0.2e1);
//     double t24 = gsl_sf_erf(t11 * t21 / 0.2e1);
//     double t26 = 0.1e1 / t21;
//     double t29 = 0.1e1 / t10 * t14 * t6 * t3 * t1 - t26 * t24 * t1 + t26 * t1;
//     return t29 - 1;
// }

// double find_shift_param(const vector<double>& a, const vector<double>& b)
// {
//     gsl_root_fsolver *solver;
//     gsl_function F;
//     int iter = 0;
//     int status = 0;
//     int max_iter = 100;
//     double lb, ub;
//     double shift_par;
//     struct param_structure params = {a, b};
//     double bounds[2];
    

//     F.function = shift_param_def_func;
//     F.params = &params;
//     solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
//     gsl_root_fsolver_set(solver, &F, bounds[0], bounds[1]);
//     do {
// 	iter++;
// 	status = gsl_root_fsolver_iterate (solver);
// 	shift_par = gsl_root_fsolver_root (solver);
// 	lb = gsl_root_fsolver_x_lower(solver);
// 	ub = gsl_root_fsolver_x_upper(solver);
// 	status = gsl_root_test_interval(lb, ub, 0.0, 1.0e-4);
//     } while (status == GSL_CONTINUE && iter < max_iter);
//     gsl_root_fsolver_free(solver);
//     assert(status == GSL_SUCCESS);
//     return shift_par;
// }

double norm_fun(double x, const vector<double>& a, const vector<double>& b)
{
    double t5 = x * x;
    double t8 = (t5 >= t5 * a[1] + a[2] + b[0] ? t5 : t5 * a[1] + a[2] + b[0]);
    return t8;
}

double shifted_dist_semi_pdf(double x,
			     const vector<double>& a,
			     const vector<double>& b)
{
    double t1 = sqrt(2 * M_PI);
    double t12 = exp(norm_fun(x, a, b) * SHIFT_PARAM - x * x / 2);
    double t14 = t12 / t1;
    return t14;
}

double draw_from_shifted_dist(const vector<double>& a,
			      const vector<double>& b)
{
    double sigma = sqrt(0.5/(0.5 - SHIFT_PARAM));
    uniform_real_distribution<double> unif(0, 1);
    normal_distribution<double> dist(0, sigma);
    double rv = 0;
    bool accepted = false;
    double C = exp(SHIFT_PARAM * (a[2] + b[0])) * sigma;
    while (! accepted) {
	double X = dist(gen);
	double t1 = shifted_dist_semi_pdf(X, a, b);
	double t2 = gsl_ran_gaussian_pdf(X, sigma);
	
	double prob = t1 / (t2 * C);
	double U = unif(gen);
	if (U <= prob) {
	    rv = X;
	    accepted = true;
	}
    }
    assert(rv != 0);
    return rv;
}

XMatrix& gen_rand_matrix(const vector<double> &alpha,
		     const vector<double> &beta,
		     double z2, XMatrix &M)
{
    unsigned p = beta.size(), q = alpha.size() - 1;
    unsigned n = p + q - 1;

    M.set_size(n, n);
    M(0, 0) = alpha[1] * z2 + beta[0];
    for (unsigned i = 1; i < p; i++) {
	M(0, i) = beta[i];
    }
    for (unsigned i = p; i < n; i++) {
	M(0, i) = alpha[i - p + 2];
    }
    for (unsigned i = 1; i < p; i++) {
	M(i, i - 1) = 1;
    }
    if (n > 1) M(p, 0) = z2;
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

#pragma omp parallel for shared(gen)
    for (unsigned k = 0; k < K; k++) {
	vector<double> Z(n);
	vector<XMatrix> A(n);
	for (unsigned i = 0; i < n; i++) {
#pragma omp critical
	    Z[i] = draw_from_shifted_dist(alpha, beta);
	    gen_rand_matrix(alpha, beta, Z[i] * Z[i], A[i]);
	}
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

	ExtremeNumber t;
	t = accumulate(Z.begin(), Z.end(), t=0,
		       [&](const ExtremeNumber& s, const ExtremeNumber& a ) {
			   return s + exp(SHIFT_PARAM * norm_fun(a, alpha, beta))/(double)n;
		       });
	for (unsigned i = 0; i < n; i++) {
	    results[k] *= exp(-SHIFT_PARAM * norm_fun(Z[i], alpha, beta));
	}
    }

    ExtremeNumber mean;
    mean = accumulate(results.cbegin(), results.cend(), mean=ExtremeNumber(0),
		      [K](const ExtremeNumber& a, const ExtremeNumber& b) {
			  return a + b / (double)K;
		      });
    sd = accumulate(results.cbegin(), results.cend(), sd=0,
		    [&](const ExtremeNumber& a, const ExtremeNumber& b) {
			return a + (abs(b - mean)^2) / (double)K;
		    });
    sd ^= 0.5;
    return mean;
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

