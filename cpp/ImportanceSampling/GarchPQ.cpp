#include <stdio.h>
#include <time.h>
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <iostream>
#include <random>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

using namespace std;
using namespace arma;

// random_device gen;
unsigned long SEED;

mat& gen_rand_matrix(const vector<double> &alpha,
		     const vector<double> &beta,
		     double z2, mat &M)
{
    unsigned p = beta.size(), q = alpha.size() - 1;
    unsigned n = p + q - 1;

    M.zeros(n, n);
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

double estimateLambda(const vector<double>& a,
		      const vector<double>& b,
		      double theta,
		      unsigned long N,
		      unsigned long K,
		      double &sd)
{

    vector<vec> E(K);
    vector<vec> E_prev(K);
    vector<double> alpha(K, 1);
#pragma omp parallel for
    for (unsigned i = 0; i < K; i++) {
	uniform_real_distribution<double> unif;
	mt19937 gen(i + K + SEED);
	E_prev[i].set_size(a.size() + b.size() - 2);
	for_each(E_prev[i].begin(), E_prev[i].end(),
		 [&](double &x)
		 {
		     x = unif(gen);
		 });
	auto n = norm(E_prev[i], "Inf");
	E_prev[i] /= n;
    }
    double Lambda = 0;
    // According to me
    // vector<double> Lambda_samples(K, 1);
    sd = 0;
    for (unsigned j = 0; j < N; j++) {
	uniform_real_distribution<double> unif;
	chi_squared_distribution<double> chi2;
	mt19937 gen(j + SEED);

	vector<mat> A(K);
	vector<double> Q(K);
	partial_sum(alpha.begin(), alpha.end(), Q.begin());

#pragma omp parallel for schedule(dynamic) shared(gen, chi2)
	for (unsigned k = 0; k < K; k++) {
	    double z2;
#pragma omp critical
	    z2 = chi2(gen);
	    gen_rand_matrix(a, b, z2, A[k]);
	}

#pragma omp parallel for shared(gen, unif)
	for (unsigned k = 0; k < K; k++) {
	    double U;
#pragma omp critical
	    U = unif(gen) * Q.back();
	    unsigned l = upper_bound(Q.begin(), Q.end(), U) - Q.begin();
	    E[k] = A[k] * E_prev[l];
	    alpha[k] = norm(E[k], "inf");
	    E[k] /= alpha[k];
	    alpha[k] = pow(alpha[k], theta);

	    // according to me
	    // Lambda_samples[k] *= log(alpha[k]);
	}
	copy(E.begin(), E.end(), E_prev.begin());
	// According to Anand
	// double lbt = log(Q.back()/K);
	Lambda += log(Q.back()/K) / N;
	sd += pow(lbt, 2) / N;
    }
    // According to Anand
    double lbt = log(accumulate(alpha.begin(), alpha.end(), 0.0)/K);
    Lambda += lbt/N;
    sd += pow(lbt, 2)/N;
    sd = sqrt(sd - pow(Lambda, 2));
    return Lambda;

    // according to me
    // Lambda = 
    // accumulate(Lambda_samples.begin(), Lambda_samples.end(), 0.0,
    // 	       [=](double s, double x) {
    // 		   return s + x/K;
    // 	       });
    // sd =
    // accumulate(Lambda_samples.begin(), Lambda_samples.end(), 0.0,
    // 	       [=](double s, double x) {
    // 		   double t = exp(x) - Lambda;
    // 		   return s + t * t/K;
    // 	       });
    // sd = sqrt(sd);
    // return log(Lambda)/N;
}

struct func_par
{
    const vector<double>& a;
    const vector<double>& b;
    unsigned long N;
    unsigned long K;
};

double func(double theta, void *p)
{
    struct func_par* par = (struct func_par*) p;
    double sd;
    return estimateLambda(par->a, par->b, theta, par->N, par->K, sd);
}

double find_root(const vector<double>& a,
		 const vector<double>& b,
		 unsigned long N,
		 unsigned long K,
		 double bounds[2])
{
    gsl_root_fsolver *solver;
    gsl_function F;
    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double lb, ub;
    double xi;
    struct func_par par = {
	a, b, N, K
    };
    F.function = func;
    F.params = &par;
    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(solver, &F, bounds[0], bounds[1]);
    do {
	iter++;
	status = gsl_root_fsolver_iterate (solver);
	xi = gsl_root_fsolver_root (solver);
	lb = gsl_root_fsolver_x_lower(solver);
	ub = gsl_root_fsolver_x_upper(solver);
	status = gsl_root_test_interval(lb, ub, 0.0, 1.0e-4);
	if (status == GSL_SUCCESS)
	    cout << "Tail index found: xi = " << xi << endl;
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(solver);
    if (status != GSL_SUCCESS) {
	cout << "The Brent algorithm did not converge after " << max_iter
	     << " iterations." << endl;
	xi = -1;
    }
    return xi;
}

int main(int argc, char*argv[])
{
    // // madeup
    // vector<double> alpha({1.0e-7, 0.11, 0});
    // vector<double> beta({0.88});
    // DAX
    // vector<double> alpha({1.0e-7, 0.02749864, 0.04228535});
    // vector<double> beta({0.8968533});
    // FTSE100
    // vector<double> alpha({1.0e-7, 0.10623464, 0.02904907});
    // vector<double> beta({0.7829784});
    // // FTSE100 GARCH(1,1)
    // vector<double> alpha({1.0e-7, 1.268747e-01});
    // vector<double> beta({8.008847e-01});

    // SP500
    vector<double> alpha({9.225747e-06, 8.835834e-02, 9.685783e-02});
    vector<double> beta({6.543018e-01});
    
    double Lambda;
    
    // cout << "alpha[0]= "  << alpha[0] << ", alpha[1]=" <<
    // 	alpha[1] << ", alpha[2]=" << alpha[2] <<
    // 	", beta[1]=" << beta[0] << endl;
    cout << "N = " << argv[3] << endl;
    cout << "K = " << argv[4] << endl;

    unsigned long N = stoul(argv[3]), K = stoul(argv[4]);
    SEED = stoul(argv[5]);
    
    double bounds[2];
    bool flag = false;
    for (double nu = stod(argv[1]); nu <= stod(argv[2]); nu += 0.1) {
	double sd;
	Lambda = estimateLambda(alpha, beta, nu, N, K, sd);
	// according to Anand
	printf("%.4f    % .4f    %.4f    %.4f\n", nu, Lambda, sd, sd/abs(Lambda));

	// according to me
	// printf("%.4f    % .4f    %.4f    %.4f\n", nu, exp(Lambda), sd, sd/exp(Lambda));
	if (!flag && Lambda > 0) {
	    bounds[1] = nu;
	    bounds[0] = nu - 0.1;
	    flag = true;
	}
    }
    double xi = find_root(alpha, beta, N, K, bounds);
    cout << "Lambda(" << xi << ") = 0" << endl;
    return 0;
}


