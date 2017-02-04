#include <stdio.h>
#include <assert.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <string>
#include <armadillo>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_eigen.h>
#include <nlopt.hpp>

#define ARRAY_SIZE(X) sizeof(X)/sizeof((X)[0])

using namespace std;
using namespace arma;

double t_density_function(double x, void *params)
{
    double nu = *(double *)params;
    return gsl_ran_tdist_pdf(x, nu);
}

double norm_density_function(double x, void *params)
{
    return gsl_ran_gaussian_pdf(x, 1);
}

/**
   params:
   0th: tail index
   1st: prob. density func.
   2nd: a vector of ARCH and GARCH parameters
 */
double uni_moment_integrand_function (double t, void *params)
{
    vector<void *> *par = (vector<void *> *)params;
    double exponent = *(double *)(par->at(0));
    gsl_function *pdf = (gsl_function *)(par->at(1));
    vector<double> *mdl_par = (vector<double> *)(par->at(2));
    double y;
    if (mdl_par->size() == 1) {
	y = pow(t * t, exponent) *
	    pdf->function(t, pdf->params);
    } else {
	y = pow(t * t + mdl_par->at(1) / mdl_par->at(0), exponent) *
	    pdf->function(t, pdf->params);
    }
    return y;
}

/**
   params:
   0th: a vector of ARCH and GARCH parameters
   1st: prob. density func.
 */
double uni_moment_function(double index, void *params)
{
    vector<void *> *par = (vector<void *> *)params;
    // double a = ((vector<double> *)par->at(0))->at(0);
    gsl_integration_workspace *w 
	= gsl_integration_workspace_alloc (1000);
  
    double result, error;
    gsl_function G;
    vector<void *> G_params;
    vector<double> *garch = (vector<double> *)par->at(0);
    
    G_params.push_back(&index);
    G_params.push_back(par->at(1));
    G_params.push_back(par->at(0));
    
    G.function = &uni_moment_integrand_function;
    G.params = &G_params;
    gsl_integration_qagi(&G, 0, 1.0e-3, 1000, w, &result, &error);
    return pow(garch->at(0), index) * result - 1;
}

double uni_tail_index(vector<double>* mdl_par, gsl_function *pdf)
{
    gsl_root_fsolver *solver;
    gsl_function F;
    vector<void *> F_params;

    int iter = 0;
    int status = 0;
    int max_iter = 100;
    double lb, ub;
    double xi = 2;

    F.function = &uni_moment_function;
    F_params.push_back(mdl_par);
    F_params.push_back(pdf);
    F.params = &F_params;

    // double test[10];
    // xi = 0.1;
    // for (size_t i = 0; i < ARRAY_SIZE(test); i++) {
    // 	test[i] = uni_moment_function(xi += i, &F_params);
    // 	if (i > 0 && ((test[i-1] < 0 && test[i] > 0) ||
    // 		      (test[i-1] > 0 && test[i] < 0))) {
    // 	    lb = xi - 1;
    // 	    ub = xi;
    // 	    break;
    // 	}
    // }
    
    double t = uni_moment_function(xi, &F_params);
    if (t < 0) {
    	do {
    	    xi *= 2;
    	    t = uni_moment_function(xi, &F_params);
    	} while (t < 0);
    	ub = xi;
    	lb = xi / 2;
    } else {
    	do {
    	    xi /= 2;
    	    t = uni_moment_function(xi, &F_params);
    	} while (t > 0);
    	ub = xi * 2;
    	lb = xi;
    }
    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(solver, &F, lb, ub);
    do {
	iter++;
	status = gsl_root_fsolver_iterate (solver);
	xi = gsl_root_fsolver_root (solver);
	lb = gsl_root_fsolver_x_lower(solver);
	ub = gsl_root_fsolver_x_upper(solver);
	status = gsl_root_test_interval(lb, ub, 0.0, 1.0e-4);
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(solver);
    assert(status == GSL_SUCCESS);
    return xi;
}

// arg 0: estimated values of alpha1, beta1
// arg 1: inverse of the covariance matrix of alpha1, beta1
// arg 2: Hill estimate
double objective_func(const vector<double> &X, vector<double> &grad,
		      void * params)
{
    vector<void *> *par = (vector<void *> *)params;
    const vec *estimated = (const vec *)par->at(0);
    const mat *cov_inv = (const mat *)par->at(1);
    vec arg(X);
    double density = 0;

    // compute the multivariate normal density function
    arg = (arg - *estimated);
    
    density = exp(-as_scalar(arg.t() * *cov_inv * arg)/2);
    density *= sqrt(abs(det(*cov_inv)));
    density /= 2 * M_PI;

    return density;
}

double constraint_func(const vector<double> &X,
		      vector<double> &grad,
		      void * params)
{
    double index = *(double *)params;
    vector<double> arg(X);
    gsl_function F;
    F.function = norm_density_function;
    F.params = NULL;
    double alpha = uni_tail_index(&arg, &F);
    return alpha - index;
}

// params[0]: alpha
// params[1]: Hill estimate
double beta_function(double beta, void *params)
{
    vector<double> *par = (vector<double> *)params;
    vector<double> mdl_par(2);
    mdl_par[0] = par->at(0);
    mdl_par[1] = beta;

    gsl_function pdf;
    pdf.function = norm_density_function;
    pdf.params = NULL;

    double index = uni_tail_index(&mdl_par, &pdf);
    return index - par->at(1);
}

void initial_compliant_params(vector<double> &garch,
			      double hill)
{
    gsl_root_fsolver *solver;
    gsl_function F;
    vector<double> params(2);
    params[0] = garch[0];
    params[1] = hill;
    
    F.function = &beta_function;
    F.params = &params;

    double beta = garch[1];
    double lb, ub;
    double t = beta_function(beta, &params);
    if (t > 0) {
	do {
	    beta += (1 - beta) / 2;
	    t = beta_function(beta, &params);
	} while (t > 0);
	lb = 2 * beta - 1;
	ub = beta;
    } else {
	do {
	    beta /= 2;
	    t = beta_function(beta, &params);	    
	} while (t < 0);
	lb = beta;
	ub = beta * 2;
    }
    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(solver, &F, lb, ub);
    int iter = 0;
    int status;
    int max_iter = 100;
    do {
	iter++;
	status = gsl_root_fsolver_iterate (solver);
	beta = gsl_root_fsolver_root (solver);
	lb = gsl_root_fsolver_x_lower(solver);
	ub = gsl_root_fsolver_x_upper(solver);
	status = gsl_root_test_interval(lb, ub, 0.0, 1.0e-4);
    } while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free(solver);
    assert(status == GSL_SUCCESS);
    garch[1] = beta;
}

// maximize the prob. density while making the implied tail index
// equal to the Hill estimate
vector<double> garch_refine_params(vector<double> params,
				   vector<double> params_cov,
				   double tail_index)
{
    size_t iter = 0;
    int status;

    gsl_multimin_fminimizer *s =
	gsl_multimin_fminimizer_alloc(
	    gsl_multimin_fminimizer_nmsimplex2, 3);
    vector<double> X(params);
    initial_compliant_params(X, tail_index);

    nlopt::opt opt(nlopt::LN_COBYLA, 2);

    // return X;


// alpha1, beta1, lambda: Lagrange Multiplier
    // vector<double> X(2);
    // gsl_vector_set(X, 0, X[0]);
    // gsl_vector_set(X, 1, X[1]);
    // gsl_vector_set(X, 2, 1);

    vec estimated(params);
    mat cov_inv(params_cov);
    cov_inv.reshape(2, 2);
    cov_inv = cov_inv.i();

    // gsl_function density_func;
    // density_func.function = &norm_density_function;
    // density_func.params = NULL;
    
    vector<void *> objective_params(2);
    objective_params[0] = &estimated;
    objective_params[1] = &cov_inv;

    opt.set_max_objective(objective_func, &objective_params);
    opt.add_equality_constraint(constraint_func, &tail_index, 0.5e-3);
    opt.set_xtol_rel(0.5e-3);

    double density;
    nlopt::result result = opt.optimize(X, density);
    // gsl_multimin_function target_func;
    // target_func.n = 3;
    // target_func.f = objective_func;
    // target_func.params = &objective_params;

    // gsl_vector *ss = gsl_vector_alloc(3);
    // gsl_vector_set(ss, 0, 0.5 * sqrt(params_cov[0]));
    // gsl_vector_set(ss, 1, 0.5 * sqrt(params_cov[3]));
    // gsl_vector_set(ss, 2, 0.5);
    
    // gsl_multimin_fminimizer_set(s, &target_func, X, ss);
    // double size;
    // do {
    // 	iter++;
    // 	status = gsl_multimin_fminimizer_iterate(s);
    // 	if (status) break;

    // 	size = gsl_multimin_fminimizer_size (s);
    // 	status = gsl_multimin_test_size (size, 1e-2);

    // 	if (status == GSL_SUCCESS) {
    // 	    printf ("converged to minimum at\n");
    //     }
    // } while (status == GSL_CONTINUE && iter < 100);
    // gsl_vector *new_estimate = gsl_multimin_fminimizer_x(s);
    // vector<double> refined(2);
    // refined[0] = gsl_vector_get(new_estimate, 0);
    // refined[1] = gsl_vector_get(new_estimate, 1);
    
    // gsl_vector_free(ss);
    // gsl_vector_free(new_estimate);
    // gsl_multimin_fminimizer_free(s);
    return X;
}

int main(void)
{

    //a1, b1, Hill Estimate of the squares
    double mdl_par[][3] = {
	{0.04029375, 0.9437317, 2.141817},
	{0.04747335, 0.9353835, 2.128157},
	{0.05731124, 0.9252467, 2.261255},
	{0.04494217, 0.9334640, 1.914272},
	{0.08561078, 0.8715603, 2.438216},
	{0.03269157, 0.9554981, 2.094851},
	{0.03125818, 0.9320087, 1.899987},
	{0.04317779, 0.9360484, 2.131553},
	{0.03197846, 0.9575640, 2.028524}
    };

    // covariance matrix of a1, b1
    double param_cov[][4] = {
	{1.579851e-04, -2.358054e-04, -2.358054e-04, 4.088979e-04},
	{1.171753e-04, -1.526732e-04, -1.526732e-04, 2.679860e-04},
	{1.924188e-04, -2.484854e-04, -2.484854e-04, 4.034176e-04},
	{9.585965e-05, -1.144427e-04, -1.144427e-04, 2.186690e-04},
	{2.640917e-04, -3.590647e-04, -3.590647e-04, 7.224615e-04},
	{5.208773e-05, -6.319880e-05, -6.319880e-05, 1.053359e-04},
	{1.248856e-04, -2.809559e-04, -2.809559e-04, 8.419657e-04},
	{1.185908e-04, -1.548410e-04, -1.548410e-04, 2.848309e-04},
	{4.970730e-05, -5.603032e-05, -5.603032e-05, 9.138853e-05}
    };
    
    // Compute the implied tail index
    // gsl_function pdf;
    // vector<double> garch;
    // double kappa;
    // pdf.function = &norm_density_function;
    // for (size_t i = 0; i < ARRAY_SIZE(mdl_par); i++) {
    // 	garch.assign(mdl_par[i], mdl_par[i] + ARRAY_SIZE(mdl_par[i]) - 1);
    // 	pdf.params = mdl_par[i] + ARRAY_SIZE(mdl_par[i]) - 1;
    // 	kappa = uni_tail_index(&garch, &pdf);
    // 	printf("%e\n", kappa);
    // }

    for (size_t i = 0; i < ARRAY_SIZE(mdl_par); i++) {
    	vector<double> estimated(mdl_par[i], mdl_par[i] + 2);
    	vector<double> cov_mat(param_cov[i], param_cov[i] + 4);
    	double hill_estimate = mdl_par[i][2];
    	vector<double> refined = garch_refine_params(
    	    estimated, cov_mat, hill_estimate);
    	printf("%e\t%e\t%8.3f\n", refined[0], refined[1],
	       refined[0] + refined[1]);
    }

    // vector<double> temp(param_cov[0], param_cov[1]);
    // vec V(temp);
    // mat M = V * V.t();
    // M.print("M");
    return 0;
}
