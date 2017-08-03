#include <stdio.h>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <armadillo>
#include <iostream>
#include <random>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include "GarchPQ.hpp"

using namespace std;
using namespace arma;

GarchPQ::GarchPQ(const vector<double> &alpha, const vector<double> &beta)
    :alpha(alpha), beta(beta), gen(time(NULL))
{
}


mat& GarchPQ::gen_rand_matrix(double z2, mat &M) const
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

/**
 * DO NOT parallelize the code because parallel execution
 * leads to different estimate for the same theta
 */
double GarchPQ::estimateLambda
(double theta, double &sd, unsigned long N, unsigned long K) const
{
    gen.seed(floor(theta));
    vector<vec> E(K);
    vector<vec> E_prev(K);
    vector<double> alpha(K, 1);
// #pragma omp parallel for
    for (unsigned i = 0; i < K; i++) {
	uniform_real_distribution<double> unif;
	E_prev[i].set_size(this->alpha.size() + this->beta.size() - 2);
	for_each(E_prev[i].begin(), E_prev[i].end(),
		 [&](double &x)
		 {
		     x = unif(gen);
		 });
	auto n = norm(E_prev[i], "Inf");
	E_prev[i] /= n;
    }
    double Lambda = 0;
    sd = 0;
    for (unsigned j = 0; j < N; j++) {
	uniform_real_distribution<double> unif;
	chi_squared_distribution<double> chi2;
	vector<mat> A(K);
	vector<double> Q(K);
	partial_sum(alpha.begin(), alpha.end(), Q.begin());

// #pragma omp parallel for schedule(dynamic) shared(chi2)
	for (unsigned k = 0; k < K; k++) {
	    double z2;
	    z2 = chi2(gen);
	    gen_rand_matrix(z2, A[k]);
	}

// #pragma omp parallel for shared(unif)
	for (unsigned k = 0; k < K; k++) {
	    double U;
	    U = unif(gen) * Q.back();
	    unsigned l = upper_bound(Q.begin(), Q.end(), U) - Q.begin();
	    E[k] = A[k] * E_prev[l];
	    alpha[k] = norm(E[k], "inf");
	    E[k] /= alpha[k];
	    alpha[k] = pow(alpha[k], theta);

	}
	copy(E.begin(), E.end(), E_prev.begin());
	// According to Anand
	double lbt = log(Q.back()/K);
	Lambda += lbt / N;
	sd += pow(lbt, 2) / N;
    }
    // According to Anand
    double lbt = log(accumulate(alpha.begin(), alpha.end(), 0.0)/K);
    Lambda += lbt/N;
    sd += pow(lbt, 2)/N;
    sd = sqrt(sd - pow(Lambda, 2));
    return Lambda;

}

struct func_par
{
    // const vector<double>& a;
    // const vector<double>& b;
    const GarchPQ *garch;
    unsigned long N;
    unsigned long K;
};

double func(double theta, void *p)
{
    struct func_par* par = (struct func_par*) p;
    const GarchPQ *garch = par->garch;
    double sd;
    return garch->estimateLambda(theta, sd, par->N, par->K);
}

/*
 * We add more and more estimation points until their 
 */
double GarchPQ::find_tail_index
(double bounds[2], unsigned long N, unsigned long K) const
{
    gsl_root_fsolver *solver;
    gsl_function F;
    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double lb, ub;
    double xi;
    struct func_par par = {
	this, N, K
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
	status = gsl_root_test_interval(lb, ub, 0.0, 1.0e-2);
	// if (status == GSL_SUCCESS)
	//     cout << "Tail index found: xi = " << xi << endl;
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(solver);
    if (status != GSL_SUCCESS) {
	cout << "The Brent algorithm did not converge after " << max_iter
	     << " iterations." << endl;
	xi = -1;
    }
    return xi;
}

#ifndef MAIN_ALGO
int main(int argc, char*argv[])
{
    // // madeup
    // vector<double> alpha({1.0e-7, 0.11, 0});
    // vector<double> beta({0.88});
    // vector<double> alpha({1.0e-7, 0.2, 0.1});
    // vector<double> beta({0.2});
    // DAX
    // vector<double> alpha({3.374294e-06, 2.074732e-02, 4.104949e-02});
    // vector<double> beta({9.102376e-01});
    // FTSE100
    // vector<double> alpha({1.0e-7, 0.10623464, 0.02904907});
    // vector<double> beta({0.7829784});
    // // FTSE100 GARCH(1,1)
    // vector<double> alpha({1.0e-7, 1.268747e-01});
    // vector<double> beta({8.008847e-01});

    // DJIA
    vector<double> alpha({3.374294e-06, 0.061577928, 0.12795424});
    vector<double> beta({0.6610499});

    // SP500
    // vector<double> alpha({9.376992e-06, 7.949678e-02, 8.765884e-02});
    // vector<double> beta({6.683833e-01});
    
    GarchPQ garch(alpha, beta);
    double Lambda;
    
    // cout << "alpha[0]= "  << alpha[0] << ", alpha[1]=" <<
    // 	alpha[1] << ", alpha[2]=" << alpha[2] <<
    // 	", beta[1]=" << beta[0] << endl;
    cout << "N = " << argv[3] << endl;
    cout << "K = " << argv[4] << endl;

    unsigned long N = stoul(argv[3]), K = stoul(argv[4]);
    // SEED = stoul(argv[5]);
    
    double bounds[2];
    bool flag = false;
    for (double nu = stod(argv[1]); nu <= stod(argv[2]); nu += 0.1) {
	double sd;
	Lambda = garch.estimateLambda(nu, sd, N, K);
	// according to Anand
	printf("%.4f    % .4f    %.4f    %.4f\n",
	       nu, Lambda, sd, sd/abs(Lambda));

	if (!flag && Lambda > 0) {
	    bounds[1] = nu;
	    bounds[0] = nu - 0.1;
	    flag = true;
	}
    }
    double xi = garch.find_tail_index(bounds, N, K);
    cout << "Lambda(" << xi << ") = 0" << endl;
    return 0;
}
#endif

