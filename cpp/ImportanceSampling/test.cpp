#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <math.h> 
#include <algorithm>
#include <armadillo>
#include <iostream>
#include <random>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_min.h>
#include "ExtremeNumber2.hpp"
#include "XMatrix.hpp"

using namespace std;
using namespace arma;

#define URNG mt19937

struct param_set1
{
    const vector<double>& alpha;
    const vector<double>& beta;
    double xi;
    unsigned long n;
    unsigned long K;
    ExtremeNumber value;
};

double norm_const(double x,
		  const vector<double>& a,
		  const vector<double>& b)
{
    double w = (a[2] + b[0])/(1 - a[1]);
    double t1 = sqrt(0.2e1);
    double t3 = exp(x * a[2]);
    double t6 = exp(x * b[0]);
    double t10 = sqrt(-0.4e1 * x * a[1] + 0.2e1);
    double t11 = sqrt(w);
    double t14 = gsl_sf_erf(t11 * t10 / 0.2e1);
    double t21 = sqrt(-0.4e1 * x + 0.2e1);
    double t24 = gsl_sf_erf(t11 * t21 / 0.2e1);
    double t26 = 0.1e1 / t21;
    double t29 = 0.1e1 / t10 * t14 * t6 * t3 * t1 - t26 * t24 * t1 + t26 * t1;
    return t29;
}

double norm_fun(double x, const vector<double>& a, const vector<double>& b)
{
    double t5 = x * x;
    double t8 = (t5 >= t5 * a[1] + a[2] + b[0] ? t5 : t5 * a[1] + a[2] + b[0]);
    return t8;
}

double shifted_dist_semi_pdf(double x, double shift_param,
			     const vector<double>& a,
			     const vector<double>& b)
{
    double t1 = sqrt(2 * M_PI);
    double t12 = exp(norm_fun(x, a, b) * shift_param - x * x / 2);
    double t14 = t12 / t1;
    return t14;
}

double draw_from_shifted_dist(const vector<double>& a,
			      const vector<double>& b,
			      double shift_param,
			      URNG& gen)
{
    double sigma = sqrt(0.5/(0.5 - shift_param));
    uniform_real_distribution<double> unif(0, 1);
    normal_distribution<double> dist(0, sigma);
    double rv = 0;
    bool accepted = false;
    double C = exp(shift_param * (a[2] + b[0])) * sigma;
    while (! accepted) {
	double X = dist(gen);
	double t1 = shifted_dist_semi_pdf(X, shift_param, a, b);
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
			     double shift_param,
			     double xi,
			     unsigned long n,
			     unsigned long K,
			     ExtremeNumber &sd)
{
    vector<ExtremeNumber> results(K);

#pragma omp parallel for
    for (unsigned k = 0; k < K; k++) {
	mt19937 gen(k);
	vector<double> Z(n);
	vector<XMatrix> A(n);
	for (unsigned i = 0; i < n; i++) {
	    Z[i] = draw_from_shifted_dist(alpha, beta, shift_param, gen);
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

	ExtremeNumber C = norm_const(shift_param, alpha, beta);

	for (unsigned i = 0; i < n; i++) {
	    results[k] *= C * exp(-shift_param * norm_fun(Z[i], alpha, beta));
	}
    }

    ExtremeNumber mean;
    mean = accumulate(results.cbegin(), results.cend(), mean=0,
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

double fun1(double shift_param, void *param)
{
    struct param_set1* par = (struct param_set1 *)param;
    ExtremeNumber sd;
    par->value = estimateLambda(par->alpha, par->beta, shift_param, par->xi, par->n, par->K, sd);
    return sd.mylog - par->value.mylog;
}

ExtremeNumber refined_estimate(const vector<double>& alpha,
			       const vector<double>& beta,
			       double xi,
			       unsigned long n,
			       unsigned long K,
			       ExtremeNumber &sd)
{
    int status;
    int iter = 0, max_iter = 100;
    gsl_min_fminimizer *s;
    double m = 0.5 - 1/(xi + 2);
    double a = 0, b = 0.45;
    gsl_function F;

    struct param_set1 params = {
    	alpha, beta, xi, n, K, ExtremeNumber()
    };
    F.function = &fun1;
    F.params = &params;
    
    vector< array<double, 2> > fs;
    fs.push_back({0, fun1(0, &params)});
    while (fs.back()[1] >= fs.front()[1]) {
    	a += 0.1 - fs.back()[0] / 4.9;
    	fs.push_back({a, fun1(a, &params)});
    }
    fs.erase(fs.begin(), fs.end()-1);
    b = a;
    while (fs.back()[1] <= fs.front()[1]) {
    	b += 0.1 - fs.back()[0] / 4.9;
    	fs.push_back({b, fun1(b, &params)});
    }
    s = gsl_min_fminimizer_alloc(gsl_min_fminimizer_brent);
    gsl_min_fminimizer_set(s, &F, a, 0, b);

    do {
    	iter++;
    	status = gsl_min_fminimizer_iterate (s);

    	m = gsl_min_fminimizer_x_minimum (s);
    	a = gsl_min_fminimizer_x_lower (s);
    	b = gsl_min_fminimizer_x_upper (s);

    	status = gsl_min_test_interval(a, b, 0.001, 0.0);

    } while (status == GSL_CONTINUE && iter < max_iter);
    assert(status == GSL_SUCCESS);
    sd.mylog = m;
    return params.value;
}

double func(double theta, void *p)
{
    struct param_set1* par = (struct param_set1*) p;
    ExtremeNumber sd;
    ExtremeNumber expected = refined_estimate(par->alpha, par->beta, theta, par->n, par->K, sd);
    par->value = sd;
    return log(expected)/par->n;
}

double find_root(const vector<double>& a,
		 const vector<double>& b,
		 unsigned long N,
		 unsigned long K,
		 double bounds[2],
		 ExtremeNumber& sd)
{
    gsl_root_fsolver *solver;
    gsl_function F;
    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double lb, ub;
    double xi;
    struct param_set1 par = {
	a, b, 0, N, K, ExtremeNumber()
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
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(solver);
    assert (status == GSL_SUCCESS);
    sd = par.value;
    return xi;
}

int main(int argc, char*argv[])
{
    vector<double> alpha({9.2e-6, 0.08835834, 0.09685783});
    // vector<double> alpha({1.0e-7, stod(argv[3])});
    vector<double> beta({0.6543018});
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
    unsigned long n = stoul(argv[3]);
    unsigned long K = stoul(argv[4]);
    double bounds[2];
    bool flag = false;
    for (double nu = stod(argv[1]); nu < stod(argv[2]); nu += 0.1) {
    	ExtremeNumber expected = refined_estimate(
    	    alpha, beta, nu, n, K, sd);
	double Lambda = log(expected)/(double)n;
	double rel_err = static_cast<double>(sd / expected);
	if (!flag && Lambda > 0) {
	    bounds[1] = nu;
	    bounds[0] = nu - 0.1;
	    flag = true;
	}
    	printf("%e    %e    %.4f\n", nu, Lambda, rel_err);
    }
    
    double xi = find_root(alpha, beta, n, K, bounds, sd);
    cout << "Lambda(" << xi << ") = 0" << endl;

    return 0;
}

